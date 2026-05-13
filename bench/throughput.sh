#!/usr/bin/env bash
# microbenchmark: generate a synthetic pair set, time pbj.
# usage: ./bench/throughput.sh [N_PAIRS] [READ_LEN] [FRAGMENT_LEN]

set -u
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"
PBJ="${PBJ_BIN:-$ROOT_DIR/pbj}"

N=${1:-100000}
L=${2:-150}
FRAG=${3:-250}

OUT="$SCRIPT_DIR/out"
mkdir -p "$OUT"

if [ ! -x "$PBJ" ]; then
    echo "error: pbj binary not at $PBJ" >&2
    exit 2
fi

python3 - <<PY
import random
random.seed(17)
B='ACGT'
RC=str.maketrans('ACGT','TGCA')
def rc(s): return s.translate(RC)[::-1]
N=$N; L=$L; F=$FRAG
err=0.01
with open("$OUT/R1.fq","w") as r1, open("$OUT/R2.fq","w") as r2:
    for i in range(N):
        frag=''.join(random.choice(B) for _ in range(F))
        a=list(frag[:L]); b=list(rc(frag[-L:]))
        qa=['?']*L; qb=['?']*L
        swap={'A':'T','C':'G','G':'C','T':'A'}
        for j in range(L):
            if random.random()<err: a[j]=swap[a[j]]; qa[j]=chr(33+10)
            if random.random()<err: b[j]=swap[b[j]]; qb[j]=chr(33+10)
        r1.write(f"@r{i}/1\n{''.join(a)}\n+\n{''.join(qa)}\n")
        r2.write(f"@r{i}/2\n{''.join(b)}\n+\n{''.join(qb)}\n")
PY

echo "N=$N L=$L fragment=$FRAG overlap=$((2*L-FRAG))"
echo "-- single-threaded --"
"$PBJ" -f "$OUT/R1.fq" -r "$OUT/R2.fq" -o "$OUT/single" -t 1 2>&1 | grep -E 'merged|throughput|elapsed'
for T in 2 4 8; do
    echo "-- threads=$T --"
    "$PBJ" -f "$OUT/R1.fq" -r "$OUT/R2.fq" -o "$OUT/t$T" -t $T 2>&1 | grep -E 'merged|throughput|elapsed'
done
