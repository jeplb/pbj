#!/usr/bin/env bash
# compare pbj merge behavior against BBmerge and PEAR on a real paired
# FASTQ dataset. fails if pbj's merge rate diverges from BBmerge by more
# than the tolerance.
#
# usage: compare_bbmerge_pear.sh R1.fq[.gz] R2.fq[.gz] [tolerance]

set -u

R1="${1:?R1.fq required}"
R2="${2:?R2.fq required}"
TOL="${3:-0.02}"

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$( cd "$SCRIPT_DIR/../.." && pwd )"
PBJ="${PBJ_BIN:-$ROOT_DIR/pbj}"

if ! command -v bbmerge.sh >/dev/null 2>&1; then
    echo "skipping: bbmerge.sh not in PATH"
    exit 77
fi
if ! command -v pear >/dev/null 2>&1; then
    echo "skipping: pear not in PATH"
    exit 77
fi

OUT="$ROOT_DIR/tests/out/compare"
rm -rf "$OUT"; mkdir -p "$OUT"

echo "running BBmerge"
bbmerge.sh in1="$R1" in2="$R2" out="$OUT/bbmerge.fq" 2>"$OUT/bbmerge.log"

echo "running PEAR"
pear -f "$R1" -r "$R2" -o "$OUT/pear" >"$OUT/pear.log" 2>&1

echo "running pbj"
"$PBJ" -f "$R1" -r "$R2" -o "$OUT/pbj" -t 4 2>"$OUT/pbj.log"

python3 "$SCRIPT_DIR/compare_merge_rates.py" \
    --bbmerge "$OUT/bbmerge.fq" \
    --pear    "$OUT/pear.assembled.fastq" \
    --pbj     "$OUT/pbj.assembled.fastq" \
    --tolerance "$TOL"
