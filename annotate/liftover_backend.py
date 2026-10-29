"""One backend selector/argument builder for Picard liftover and dictionary tools."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import os
from pathlib import Path
import shlex
import shutil
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
MODES = ("auto", "gatk", "picard-jar", "picard-cli")
TOOLS = ("CreateSequenceDictionary", "LiftoverVcf")
DEFAULT_JAVA_OPTS = "-Xms1g -Xmx8g"


class BackendError(ValueError):
    pass


def java_options(raw):
    try:
        values = shlex.split(raw)
    except ValueError as exc:
        raise BackendError(f"invalid JVM option quoting: {exc}") from exc
    launch_options = {"-jar", "-cp", "-classpath", "--class-path", "-m", "--module",
                      "-p", "--module-path", "--source"}
    for value in values:
        if (not value.startswith("-") or value.split("=", 1)[0] in launch_options or
                any(ord(c) < 32 or ord(c) == 127 for c in value)):
            raise BackendError(f"expected self-contained JVM option, got {value!r}")
    return tuple(values)


def tool_options_value(options):
    """Quote tokens for the JVM's environment parser, not a shell parser."""
    quoted = []
    for value in options:
        if '"' not in value:
            quoted.append('"' + value + '"')
        elif "'" not in value:
            quoted.append("'" + value + "'")
        else:
            raise BackendError("mixed quote characters in a JVM option cannot be safely passed "
                               "through JAVA_TOOL_OPTIONS; use --backend picard-jar")
    return " ".join(quoted)


@dataclass(frozen=True)
class Backend:
    kind: str
    executable: str
    jar: str | None
    java_options: tuple[str, ...]

    def metadata(self):
        data = asdict(self)
        data["java_options"] = list(self.java_options)
        return data

    def invocation(self, tool, arguments, environ=None):
        if tool not in TOOLS:
            raise BackendError(f"unsupported Picard tool: {tool}")
        env = os.environ if environ is None else environ
        updates = {}
        if self.kind == "picard-jar":
            argv = [self.executable, *self.java_options, "-jar", self.jar, tool]
        elif self.kind == "gatk":
            # The upstream launcher uses str.split(), not shlex.split(). Refuse
            # unrepresentable options rather than silently breaking quoted paths.
            if any(any(c.isspace() for c in value) for value in self.java_options):
                raise BackendError("GATK's launcher splits JVM options on whitespace; "
                                   "use --backend picard-jar for space-containing JVM values")
            argv = [self.executable, "--java-options", " ".join(self.java_options), tool]
        elif self.kind == "picard-cli":
            argv = [self.executable, tool]
            value = tool_options_value(self.java_options)
            if value:
                existing = env.get("JAVA_TOOL_OPTIONS", "")
                updates["JAVA_TOOL_OPTIONS"] = (existing + " " + value).strip()
        else:
            raise BackendError(f"unknown backend: {self.kind}")
        return argv + [str(arg) for arg in arguments], updates

    @classmethod
    def from_metadata(cls, data):
        if (not isinstance(data, dict) or set(data) != {"kind", "executable", "jar", "java_options"}
                or data["kind"] not in MODES[1:] or not isinstance(data["executable"], str)
                or not Path(data["executable"]).is_absolute()
                or not isinstance(data["java_options"], list)
                or not all(isinstance(v, str) for v in data["java_options"])):
            raise BackendError("invalid pinned backend specification")
        options = java_options(shlex.join(data["java_options"]))
        executable = Path(data["executable"])
        if not executable.is_file() or not os.access(executable, os.X_OK):
            raise BackendError(f"backend executable is missing/not executable: {executable}")
        jar = data["jar"]
        if data["kind"] == "picard-jar":
            jar = str(readable_jar(jar)) if isinstance(jar, str) else None
            if jar is None:
                raise BackendError("pinned picard-jar backend requires a jar")
        elif jar is not None:
            raise BackendError("non-jar backend must not contain a jar")
        backend = cls(data["kind"], str(executable), jar, options)
        backend.invocation(TOOLS[0], [])  # Validate transport before starting a tool.
        return backend


def readable_jar(value):
    path = Path(value).expanduser().resolve()
    if not path.is_file() or not os.access(path, os.R_OK) or not path.stat().st_size:
        raise BackendError(f"not a readable, nonempty Picard jar: {path}")
    return path


def discover_jar(repo_root=None):
    root = SCRIPT_DIR.parent if repo_root is None else Path(repo_root)
    for name in ("mapping/tst/picard.jar", "util/bin/picard.jar", "mapping/picard.jar"):
        path = root / name
        if path.exists() or path.is_symlink():
            return readable_jar(path)
    return None


def select_backend(mode="auto", jar=None, options=DEFAULT_JAVA_OPTS, *, environ=None, repo_root=None):
    env = os.environ if environ is None else environ
    if mode not in MODES:
        raise BackendError(f"unknown backend: {mode}")
    parsed = java_options(options)
    selected_jar = None
    if mode == "auto":
        if jar:
            mode = "picard-jar"
        elif shutil.which("gatk", path=env.get("PATH", os.defpath)):
            mode = "gatk"
        else:
            selected_jar = discover_jar(repo_root)
            mode = "picard-jar" if selected_jar else "picard-cli"
    if mode == "picard-jar":
        selected_jar = readable_jar(jar) if jar else selected_jar or discover_jar(repo_root)
        if selected_jar is None:
            raise BackendError("set --picard-jar or PICARD_JAR; no Picard jar was discovered")
    command = {"picard-jar": "java", "gatk": "gatk", "picard-cli": "picard"}[mode]
    executable = shutil.which(command, path=env.get("PATH", os.defpath))
    if not executable:
        raise BackendError(f"missing backend command: {command}; install GATK/Picard or set --picard-jar")
    backend = Backend(mode, str(Path(executable).absolute()),
                      str(selected_jar) if selected_jar else None, parsed)
    backend.invocation(TOOLS[0], [], env)
    return backend


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend", choices=MODES, default=os.environ.get("LIFTOVER_BACKEND", "auto"))
    parser.add_argument("--picard-jar", default=os.environ.get("PICARD_JAR"))
    parser.add_argument("--java-opts", default=os.environ.get("PICARD_JAVA_OPTS", DEFAULT_JAVA_OPTS))
    parser.add_argument("--spec", help="trusted pinned backend JSON from --describe")
    parser.add_argument("--describe", action="store_true", help="select/validate without running a tool")
    parser.add_argument("--check-options", action="store_true", help="check mode/JVM syntax without requiring installed tools")
    parser.add_argument("tool", nargs="?", choices=TOOLS)
    parser.add_argument("arguments", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    try:
        if args.check_options:
            java_options(args.java_opts)
            return 0
        backend = (Backend.from_metadata(json.loads(args.spec)) if args.spec else
                   select_backend(args.backend, args.picard_jar, args.java_opts))
        if args.describe:
            print(json.dumps(backend.metadata()))
            return 0
        if args.tool is None:
            parser.error("a tool or --describe is required")
        arguments = args.arguments[1:] if args.arguments[:1] == ["--"] else args.arguments
        argv, updates = backend.invocation(args.tool, arguments)
        print(f"Backend: {backend.kind}; command: {shlex.join(argv)}", file=sys.stderr, flush=True)
        # Replace the shim so terminal signals and batch process-group cleanup
        # reach the backend. Standalone parent shells do not forward PID-only
        # signals automatically. Never eval options or split tool paths.
        os.execve(argv[0], argv, {**os.environ, **updates})
    except (BackendError, OSError, ValueError) as exc:
        parser.exit(2, f"ERROR: {exc}\n")


if __name__ == "__main__":
    main()
