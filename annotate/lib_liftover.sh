#!/bin/sh

# Shared helpers for liftover-oriented shell scripts.
# Source from other scripts:
#   . "$(dirname "$0")/lib_liftover.sh"

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "ERROR: missing command: $1" >&2
    exit 1
  }
}

detect_picard_jar() {
  script_dir="$1"
  for p in \
    "${script_dir}/../mapping/tst/picard.jar" \
    "${script_dir}/../util/bin/picard.jar" \
    "${script_dir}/../mapping/picard.jar"
  do
    if [ -f "$p" ]; then
      echo "$p"
      return 0
    fi
  done
  return 1
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
  ensure_ref_index "$ref" || return

  dict_primary="$(dict_path_for_ref "$ref")" || return
  dict_alt="${ref%.*}.dict"
  if [ ! -e "$dict_primary" ] && [ ! -L "$dict_primary" ] && [ ! -e "$dict_alt" ] && [ ! -L "$dict_alt" ]; then
    if command -v gatk >/dev/null 2>&1; then
      gatk CreateSequenceDictionary -R "$ref" -O "$dict_primary" || return
    elif [ -n "$picard_jar" ] && [ -f "$picard_jar" ]; then
      java $java_opts -jar "$picard_jar" CreateSequenceDictionary R="$ref" O="$dict_primary" || return
    elif command -v picard >/dev/null 2>&1; then
      picard CreateSequenceDictionary R="$ref" O="$dict_primary" || return
    else
      echo "ERROR: missing sequence dictionary for $ref and no tool available to create it." >&2
      echo "Install gatk/picard or set --picard-jar." >&2
      exit 2
    fi
  fi
  python3 "${SCRIPT_DIR}/check_reference.py" "$ref" >/dev/null || return
  link_dict_alias "$dict_primary" "$dict_alt" || return
  link_dict_alias "$dict_alt" "$dict_primary" || return
}
