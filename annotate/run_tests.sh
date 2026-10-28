#!/bin/sh
set -eu

SCRIPT_DIR=$(CDPATH='' cd -- "$(dirname -- "$0")" && pwd)
python3 "${SCRIPT_DIR}/test_compare_marker_status_sets.py"
python3 "${SCRIPT_DIR}/test_reference_utils.py"
python3 "${SCRIPT_DIR}/test_y_haplo_from_markers.py"
python3 "${SCRIPT_DIR}/test_ytree.py"
