#!/bin/sh

# Shared helpers for liftover-oriented shell scripts.
# Source from other scripts:
#   . "$(dirname "$0")/lib_liftover.sh"

load_liftover_config() {
  # Load before defaults/argument parsing: explicit CLI options must win. This
  # is trusted executable shell code, not a data-only configuration format.
  liftover_config="${BIO_TOOLS_CFG:-./bio-tools.cfg}"
  case "$liftover_config" in
    */*) ;;
    *) liftover_config="./$liftover_config" ;;
  esac
  if [ -n "${BIO_TOOLS_CFG:-}" ] && [ ! -f "$liftover_config" ]; then
    echo "ERROR: BIO_TOOLS_CFG is not a file: $liftover_config" >&2
    return 2
  fi
  [ -f "$liftover_config" ] || return 0
  case "$-" in *a*) liftover_restore_export=0 ;; *) liftover_restore_export=1 ;; esac
  set -a
  # shellcheck disable=SC1090
  . "$liftover_config"
  if [ "$liftover_restore_export" -eq 1 ]; then set +a; fi
}

find_liftover_ref() {
  liftover_allow_extended="$4"
  set -- --assembly "$1" --reference "$2" --ref-root "$3"
  [ "$liftover_allow_extended" -eq 0 ] || set -- "$@" --allow-extended-ref
  python3 "${SCRIPT_DIR}/liftover_resources.py" "$@"
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "ERROR: missing command: $1" >&2
    exit 1
  }
}

select_liftover_backend() {
  python3 "${SCRIPT_DIR}/liftover_backend.py" --describe \
    --backend "$1" --picard-jar "$2" --java-opts="$3"
}

run_liftover_tool() {
  liftover_backend_spec="$1"
  shift
  python3 "${SCRIPT_DIR}/liftover_backend.py" --spec "$liftover_backend_spec" "$@"
}

dict_path_for_ref() {
  python3 "${SCRIPT_DIR}/check_reference.py" "$1" --print-dict-path
}

link_dict_alias() {
  src="$1"
  dst="$2"
  [ -f "$src" ] || return 0
  if [ -e "$dst" ] || [ -L "$dst" ]; then
    [ -f "$dst" ] || { echo "ERROR: invalid dictionary alias: $dst" >&2; return 2; }
    return 0
  fi
  ln -s "$(basename "$src")" "$dst"
}

ensure_ref_index() {
  ref="$1"
  python3 "${SCRIPT_DIR}/check_reference.py" "$ref" --format-only >/dev/null || return
  python3 -c 'import pysam' >/dev/null 2>&1 || {
    echo "ERROR: reference validation requires pysam (python3 -m pip install pysam)" >&2
    return 1
  }
  if [ ! -e "${ref}.fai" ] && [ ! -L "${ref}.fai" ]; then
    if [ -e "${ref}.gzi" ] || [ -L "${ref}.gzi" ]; then
      echo "ERROR: existing .gzi without .fai for $ref; rebuild the reference indexes explicitly." >&2
      return 2
    fi
    samtools faidx "$ref" || return
  fi
  python3 "${SCRIPT_DIR}/check_reference.py" "$ref" --index-only >/dev/null || return
}

ensure_ref_sidecars() {
  ref="$1"
  picard_jar="$2"
  java_opts="$3"
  backend_spec="${4:-}"
  ensure_ref_index "$ref" || return

  dict_primary="$(dict_path_for_ref "$ref")" || return
  dict_alt="${ref%.*}.dict"
  if [ ! -e "$dict_primary" ] && [ ! -L "$dict_primary" ] && [ ! -e "$dict_alt" ] && [ ! -L "$dict_alt" ]; then
    if [ -z "$backend_spec" ]; then
      backend_spec="$(select_liftover_backend "${LIFTOVER_BACKEND:-auto}" "$picard_jar" "$java_opts")" || return
    fi
    run_liftover_tool "$backend_spec" CreateSequenceDictionary \
      --REFERENCE "$ref" --OUTPUT "$dict_primary" || return
  fi
  python3 "${SCRIPT_DIR}/check_reference.py" "$ref" >/dev/null || return
  link_dict_alias "$dict_primary" "$dict_alt" || return
  link_dict_alias "$dict_alt" "$dict_primary" || return
}
