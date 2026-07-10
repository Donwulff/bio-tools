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
  ref="$1"
  case "$ref" in
    *.fa.gz) echo "${ref%.fa.gz}.dict" ;;
    *.fna.gz) echo "${ref%.fna.gz}.dict" ;;
    *.fasta.gz) echo "${ref%.fasta.gz}.dict" ;;
    *.fas.gz) echo "${ref%.fas.gz}.dict" ;;
    *.fa) echo "${ref%.fa}.dict" ;;
    *.fna) echo "${ref%.fna}.dict" ;;
    *.fasta) echo "${ref%.fasta}.dict" ;;
    *.fas) echo "${ref%.fas}.dict" ;;
    *) echo "${ref%.*}.dict" ;;
  esac
}

link_dict_alias() {
  src="$1"
  dst="$2"
  [ -f "$src" ] || return 0
  [ -f "$dst" ] && return 0
  ln -sfn "$(basename "$src")" "$dst"
}

ensure_ref_sidecars() {
  ref="$1"
  picard_jar="$2"
  java_opts="$3"

  [ -f "${ref}.fai" ] || samtools faidx "$ref"

  dict_primary="$(dict_path_for_ref "$ref")"
  dict_alt="${ref%.*}.dict"
  if [ ! -f "$dict_primary" ] && [ ! -f "$dict_alt" ]; then
    if command -v gatk >/dev/null 2>&1; then
      gatk CreateSequenceDictionary -R "$ref" -O "$dict_primary"
    elif [ -n "$picard_jar" ] && [ -f "$picard_jar" ]; then
      java $java_opts -jar "$picard_jar" CreateSequenceDictionary R="$ref" O="$dict_primary"
    elif command -v picard >/dev/null 2>&1; then
      picard CreateSequenceDictionary R="$ref" O="$dict_primary"
    else
      echo "ERROR: missing sequence dictionary for $ref and no tool available to create it." >&2
      echo "Install gatk/picard or set --picard-jar." >&2
      exit 2
    fi
  fi
  # Some tools look for <ref>.dict while others expect <ref-with-fasta-suffix-removed>.dict.
  link_dict_alias "$dict_primary" "$dict_alt"
  link_dict_alias "$dict_alt" "$dict_primary"
}
