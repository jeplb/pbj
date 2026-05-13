#!/usr/bin/env bash
set -u

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"
PBJ_BIN="${PBJ_BIN:-$ROOT_DIR/pbj}"
FIX_DIR="$SCRIPT_DIR/fixtures"
OUT_DIR="$SCRIPT_DIR/out"

if [ ! -x "$PBJ_BIN" ]; then
    echo "error: pbj binary not found at $PBJ_BIN (run 'make' first)" >&2
    exit 2
fi

python3 "$SCRIPT_DIR/scripts/gen_fixtures.py" "$FIX_DIR" >/dev/null

rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"

PASS=0
FAIL=0
FAILED_NAMES=()

get_expected() {
    local file="$1"
    local key="$2"
    grep "^$key=" "$file" 2>/dev/null | head -1 | cut -d= -f2-
}

count_records() {
    local f="$1"
    if [ ! -s "$f" ]; then echo 0; return; fi
    local lines
    lines=$(wc -l < "$f" | tr -d ' ')
    echo $((lines / 4))
}

extract_seq() {
    awk 'NR==2{print; exit}' "$1"
}

extract_name() {
    # strip the leading @ and any trailing /1 /2 from the header line
    awk 'NR==1{sub(/^@/,""); print; exit}' "$1"
}

extract_qual_first_byte() {
    awk 'NR==4{printf "%d", substr($0,1,1)*1; exit}' "$1" 2>/dev/null
    # posix awk doesn't have ord() so use printf substitute:
}

ord_first() {
    python3 -c "import sys; line=open('$1').readlines()[3]; print(ord(line[0]))"
}

run_fixture() {
    local name="$1"
    local fix="$FIX_DIR/$name"
    local out="$OUT_DIR/$name"
    mkdir -p "$out"

    local exp="$fix/EXPECTED.txt"
    local extra_flags
    extra_flags=$(get_expected "$exp" pbj_flags)
    local input_ext
    input_ext=$(get_expected "$exp" input_ext)
    input_ext="${input_ext:-.fq}"
    local single_file
    single_file=$(get_expected "$exp" input_file)

    if [ -n "$single_file" ]; then
        "$PBJ_BIN" -f "$fix/$single_file" -o "$out/out" \
            ${extra_flags} 2>"$out/stderr" >"$out/stdout"
    else
        "$PBJ_BIN" -f "$fix/R1$input_ext" -r "$fix/R2$input_ext" -o "$out/out" \
            ${extra_flags} 2>"$out/stderr" >"$out/stdout"
    fi
    local rc=$?
    if [ $rc -ne 0 ]; then
        echo "  FAIL $name: pbj exit code $rc"
        sed 's/^/      /' "$out/stderr"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name"); return
    fi

    local got_merged got_unassembled
    got_merged=$(count_records "$out/out.assembled.fastq")
    got_unassembled=$(count_records "$out/out.unassembled.forward.fastq")

    local exp_merged exp_unassembled
    exp_merged=$(get_expected "$exp" merged)
    exp_unassembled=$(get_expected "$exp" unassembled)

    local ok=1
    if [ "$got_merged" != "$exp_merged" ]; then ok=0; fi
    if [ "$got_unassembled" != "$exp_unassembled" ]; then ok=0; fi

    # optional: exact merged length
    local exp_len len_ok=1
    exp_len=$(get_expected "$exp" merged_len)
    if [ -n "$exp_len" ] && [ "$got_merged" = "1" ]; then
        local got_len
        got_len=$(awk 'NR==2{print length($0); exit}' "$out/out.assembled.fastq")
        if [ "$got_len" != "$exp_len" ]; then
            len_ok=0
            echo "    length mismatch: expected $exp_len, got $got_len"
        fi
    fi

    # optional: exact merged sequence
    local exp_seq seq_ok=1
    exp_seq=$(get_expected "$exp" merged_seq)
    if [ -n "$exp_seq" ] && [ "$got_merged" = "1" ]; then
        local got_seq
        got_seq=$(extract_seq "$out/out.assembled.fastq")
        if [ "$got_seq" != "$exp_seq" ]; then
            seq_ok=0
            echo "    sequence mismatch:"
            echo "      expected: $exp_seq"
            echo "      got:      $got_seq"
        fi
    fi

    # optional: exact merged name
    local exp_name name_ok=1
    exp_name=$(get_expected "$exp" merged_name)
    if [ -n "$exp_name" ] && [ "$got_merged" = "1" ]; then
        local got_name
        got_name=$(extract_name "$out/out.assembled.fastq")
        if [ "$got_name" != "$exp_name" ]; then
            name_ok=0
            echo "    name mismatch:"
            echo "      expected: $exp_name"
            echo "      got:      $got_name"
        fi
    fi

    # optional: first quality byte as integer (validates phred+33 output)
    local exp_qual0 qual_ok=1
    exp_qual0=$(get_expected "$exp" merged_qual_first_byte)
    if [ -n "$exp_qual0" ] && [ "$got_merged" = "1" ]; then
        local got_qual0
        got_qual0=$(ord_first "$out/out.assembled.fastq")
        if [ "$got_qual0" != "$exp_qual0" ]; then
            qual_ok=0
            echo "    output quality byte mismatch: expected $exp_qual0, got $got_qual0"
        fi
    fi

    if [ "$ok" = "1" ] && [ "$seq_ok" = "1" ] && [ "$name_ok" = "1" ] && [ "$qual_ok" = "1" ] && [ "$len_ok" = "1" ]; then
        echo "  PASS $name (merged=$got_merged unassembled=$got_unassembled)"
        PASS=$((PASS+1))
    else
        echo "  FAIL $name (merged=$got_merged exp $exp_merged; unassembled=$got_unassembled exp $exp_unassembled)"
        sed 's/^/      /' "$out/stderr"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name")
    fi
}

run_negative_dovetail() {
    local name="06_dovetail_no_flag"
    local fix="$FIX_DIR/06_dovetail"
    local out="$OUT_DIR/$name"
    mkdir -p "$out"
    "$PBJ_BIN" -f "$fix/R1.fq" -r "$fix/R2.fq" -o "$out/out" 2>"$out/stderr" >"$out/stdout"
    local got_merged
    got_merged=$(count_records "$out/out.assembled.fastq")
    if [ "$got_merged" = "0" ]; then
        echo "  PASS $name (merged=0 without --dovetail)"
        PASS=$((PASS+1))
    else
        echo "  FAIL $name (merged=$got_merged, expected 0)"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name")
    fi
}

# pair-mismatch detection
run_pair_mismatch() {
    local name="pair_mismatch"
    local out="$OUT_DIR/$name"
    mkdir -p "$out"
    printf "@foo/1\nACGT\n+\n!!!!\n" > "$out/R1.fq"
    printf "@bar/2\nACGT\n+\n!!!!\n" > "$out/R2.fq"
    "$PBJ_BIN" -f "$out/R1.fq" -r "$out/R2.fq" -o "$out/o" 2>"$out/stderr" >"$out/stdout"
    local rc=$?
    if [ $rc -eq 2 ] && grep -q "read name mismatch" "$out/stderr"; then
        echo "  PASS $name (exit 2 with name mismatch)"
        PASS=$((PASS+1))
    else
        echo "  FAIL $name (exit=$rc)"
        sed 's/^/      /' "$out/stderr"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name")
    fi
}

# /1 + /1 (same explicit suffix on both sides) must NOT be treated as a pair
run_pair_same_suffix() {
    local name="pair_same_suffix"
    local out="$OUT_DIR/$name"
    mkdir -p "$out"
    printf "@foo/1\nACGT\n+\n!!!!\n" > "$out/R1.fq"
    printf "@foo/1\nACGT\n+\n!!!!\n" > "$out/R2.fq"
    "$PBJ_BIN" -f "$out/R1.fq" -r "$out/R2.fq" -o "$out/o" 2>"$out/stderr" >"$out/stdout"
    local rc=$?
    if [ $rc -eq 2 ] && grep -q "read name mismatch" "$out/stderr"; then
        echo "  PASS $name (rejected, both ended /1)"
        PASS=$((PASS+1))
    else
        echo "  FAIL $name (exit=$rc, expected 2)"
        sed 's/^/      /' "$out/stderr"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name")
    fi
}

# interleaved input with an odd record count must error
run_interleaved_odd() {
    local name="interleaved_odd"
    local out="$OUT_DIR/$name"
    mkdir -p "$out"
    printf "@r1/1\nACGT\n+\n!!!!\n@r1/2\nACGT\n+\n!!!!\n@r2/1\nACGT\n+\n!!!!\n" \
        > "$out/INTERLEAVED.fq"
    "$PBJ_BIN" -f "$out/INTERLEAVED.fq" --interleaved-in -o "$out/o" \
        2>"$out/stderr" >"$out/stdout"
    local rc=$?
    if [ $rc -eq 2 ] && grep -q "R2 ended before R1" "$out/stderr"; then
        echo "  PASS $name (odd record count detected)"
        PASS=$((PASS+1))
    else
        echo "  FAIL $name (exit=$rc)"
        sed 's/^/      /' "$out/stderr"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name")
    fi
}

# interleaved input in r2,r1 order with casava headers must be detected
run_interleaved_r2_first() {
    local name="interleaved_r2_first"
    local out="$OUT_DIR/$name"
    mkdir -p "$out"
    printf "@HISEQ:1 2:N:0:ATCG\nACGT\n+\n!!!!\n@HISEQ:1 1:N:0:ATCG\nACGT\n+\n!!!!\n" \
        > "$out/INTERLEAVED.fq"
    "$PBJ_BIN" -f "$out/INTERLEAVED.fq" --interleaved-in -o "$out/o" \
        2>"$out/stderr" >"$out/stdout"
    local rc=$?
    if [ $rc -eq 2 ] && grep -q "appears to be in (R2,R1) order" "$out/stderr"; then
        echo "  PASS $name (direction swap detected)"
        PASS=$((PASS+1))
    else
        echo "  FAIL $name (exit=$rc)"
        sed 's/^/      /' "$out/stderr"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name")
    fi
}

# adapter false positive guard: no real adapter present, only a coincidental
# 5-mer match in the middle of the read. with --adapter-min-match 5 and a
# 3'-anchored algorithm, this must NOT trim the read. before the round-2
# fix, a mid-read 5/5 match would truncate.
run_adapter_no_false_positive() {
    local name="adapter_no_false_positive"
    local out="$OUT_DIR/$name"
    mkdir -p "$out"
    # 100bp read with the adapter's first 5 bases at position 30; 3' end is
    # all-T (which doesn't match adapter prefix).
    local adapter="AGATCGGAAGAGCACACGTCT"
    local r1="CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAGATC""TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
    local r2="$(python3 -c "
import sys
s='${r1}'[:80]
tab=str.maketrans('ACGT','TGCA')
print('AAAAAAAAAAAAAAAAAAAA' + s.translate(tab)[::-1])
")"
    printf "@noadapter/1\n%s\n+\n%s\n" "$r1" "$(python3 -c "print('?'*${#r1})")" > "$out/R1.fq"
    printf "@noadapter/2\n%s\n+\n%s\n" "$r2" "$(python3 -c "print('?'*${#r2})")" > "$out/R2.fq"
    "$PBJ_BIN" -f "$out/R1.fq" -r "$out/R2.fq" -o "$out/o" \
        -a "$adapter" --adapter-min-match 5 2>"$out/stderr" >"$out/stdout"
    local rc=$?
    if [ $rc -ne 0 ]; then
        echo "  FAIL $name (pbj exit $rc)"; sed 's/^/      /' "$out/stderr"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name"); return
    fi
    local trim_r1
    trim_r1=$(grep -E '^  adapter trimmed R1' "$out/stderr" | awk '{print $4}')
    if [ "$trim_r1" = "0" ]; then
        echo "  PASS $name (no false-positive trim)"
        PASS=$((PASS+1))
    else
        echo "  FAIL $name (trimmed R1 unexpectedly: $trim_r1)"
        sed 's/^/      /' "$out/stderr"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name")
    fi
}

# g1: separate-file mode with casava 1.8+ headers passed in (r2, r1) order
run_separate_file_swap() {
    local name="separate_file_swap"
    local fix="$FIX_DIR/15_casava_header"
    local out="$OUT_DIR/$name"
    mkdir -p "$out"
    "$PBJ_BIN" -f "$fix/R2.fq" -r "$fix/R1.fq" -o "$out/o" \
        2>"$out/stderr" >"$out/stdout"
    local rc=$?
    if [ $rc -eq 2 ] && grep -q 'appears to be in (R2,R1) order' "$out/stderr"; then
        echo "  PASS $name (-f R2 -r R1 swap detected via Casava direction)"
        PASS=$((PASS+1))
    else
        echo "  FAIL $name (exit=$rc)"; sed 's/^/      /' "$out/stderr"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name")
    fi
}

# G2: --no-bonferroni smoke test
run_no_bonferroni_smoke() {
    local name="no_bonferroni_smoke"
    local fix="$FIX_DIR/01_perfect_overlap"
    local out="$OUT_DIR/$name"
    mkdir -p "$out"
    "$PBJ_BIN" -f "$fix/R1.fq" -r "$fix/R2.fq" -o "$out/o" --no-bonferroni \
        2>"$out/stderr" >"$out/stdout"
    local rc=$?
    local got_merged
    got_merged=$(count_records "$out/o.assembled.fastq")
    if [ $rc -eq 0 ] && [ "$got_merged" = "1" ]; then
        echo "  PASS $name (flag accepted, merged=1)"
        PASS=$((PASS+1))
    else
        echo "  FAIL $name (exit=$rc merged=$got_merged)"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name")
    fi
}

# G3: qtrim collapses merged read to zero with --min-length 0; must FAIL_LENGTH
run_qtrim_zero_minlen_zero() {
    local name="qtrim_zero_minlen_zero"
    local out="$OUT_DIR/$name"
    mkdir -p "$out"
    python3 - <<PY > /dev/null
import random
random.seed(101)
B='ACGT'; rc=str.maketrans('ACGT','TGCA')
frag=''.join(random.choice(B) for _ in range(60))
r1=frag[:40]; r2=frag[-40:].translate(rc)[::-1]
q=chr(33+5)*40
open('$out/R1.fq','w').write(f'@q0/1\n{r1}\n+\n{q}\n')
open('$out/R2.fq','w').write(f'@q0/2\n{r2}\n+\n{q}\n')
PY
    "$PBJ_BIN" -f "$out/R1.fq" -r "$out/R2.fq" -o "$out/o" \
        --qtrim 30 -l 0 -m 20 2>"$out/stderr" >"$out/stdout"
    local rc=$?
    local fail_len
    fail_len=$(grep '^  fail length' "$out/stderr" | awk '{print $3}')
    if [ $rc -eq 0 ] && [ "$fail_len" = "1" ]; then
        echo "  PASS $name (qtrim -> 0 forced FAIL_LENGTH with -l 0)"
        PASS=$((PASS+1))
    else
        echo "  FAIL $name (exit=$rc fail_length=$fail_len)"
        sed 's/^/      /' "$out/stderr"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name")
    fi
}

# G4: fixture-07 input WITHOUT --strict merges cleanly
run_fixture07_no_strict() {
    local name="fixture07_no_strict_merges"
    local fix="$FIX_DIR/07_strict_reject"
    local out="$OUT_DIR/$name"
    mkdir -p "$out"
    "$PBJ_BIN" -f "$fix/R1.fq" -r "$fix/R2.fq" -o "$out/o" \
        2>"$out/stderr" >"$out/stdout"
    local rc=$?
    local got_merged
    got_merged=$(count_records "$out/o.assembled.fastq")
    if [ $rc -eq 0 ] && [ "$got_merged" = "1" ]; then
        echo "  PASS $name (merges without --strict)"
        PASS=$((PASS+1))
    else
        echo "  FAIL $name (exit=$rc merged=$got_merged)"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name")
    fi
}

# G5: --adapter-min-match override below the default still 3'-anchors
run_adapter_min_match_override() {
    local name="adapter_min_match_5_override"
    local fix="$FIX_DIR/17_adapter_trim"
    local out="$OUT_DIR/$name"
    mkdir -p "$out"
    local adapter="AGATCGGAAGAGCACACGTCT"
    "$PBJ_BIN" -f "$fix/R1.fq" -r "$fix/R2.fq" -o "$out/o" \
        -a "$adapter" --adapter-min-match 5 2>"$out/stderr" >"$out/stdout"
    local rc=$?
    local got_merged
    got_merged=$(count_records "$out/o.assembled.fastq")
    if [ $rc -eq 0 ] && [ "$got_merged" = "1" ]; then
        echo "  PASS $name (adapter trim works at min-match=5)"
        PASS=$((PASS+1))
    else
        echo "  FAIL $name (exit=$rc merged=$got_merged)"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name")
    fi
}

# G6: stats counters sum to total_pairs on a multi-pair input
run_stats_sum_invariant() {
    local name="stats_sum_invariant"
    local out="$OUT_DIR/$name"
    mkdir -p "$out"
    python3 - <<PY > /dev/null
import random
random.seed(53)
B='ACGT'; rc=str.maketrans('ACGT','TGCA')
with open('$out/R1.fq','w') as a, open('$out/R2.fq','w') as b:
    for i in range(300):
        # mix of overlapping, non-overlapping, and low-quality pairs
        if i % 3 == 0:
            f=''.join(random.choice(B) for _ in range(160))
            a.write(f'@m{i}/1\n{f[:100]}\n+\n{"?"*100}\n')
            b.write(f'@m{i}/2\n{f[-100:].translate(rc)[::-1]}\n+\n{"?"*100}\n')
        elif i % 3 == 1:
            f1=''.join(random.choice(B) for _ in range(150))
            f2=''.join(random.choice(B) for _ in range(150))
            a.write(f'@n{i}/1\n{f1}\n+\n{"?"*150}\n')
            b.write(f'@n{i}/2\n{f2}\n+\n{"?"*150}\n')
        else:
            f=''.join(random.choice(B) for _ in range(160))
            a.write(f'@l{i}/1\n{f[:100]}\n+\n{chr(33+2)*100}\n')
            b.write(f'@l{i}/2\n{f[-100:].translate(rc)[::-1]}\n+\n{chr(33+2)*100}\n')
PY
    "$PBJ_BIN" -f "$out/R1.fq" -r "$out/R2.fq" -o "$out/o" \
        --max-ee 1.0 --max-mono-frac 0.95 -Q 5 2>"$out/stderr" >"$out/stdout"
    local rc=$?
    python3 - <<PY > "$out/check.txt"
import re
total = None; counters = {}
with open('$out/stderr') as f:
    for line in f:
        m = re.match(r'^  (total pairs|merged|no candidate|fail [\w -]+?)\s+(\d+)', line)
        if not m: continue
        k = m.group(1).strip(); v = int(m.group(2))
        if k == 'total pairs': total = v
        else: counters[k] = v
s = sum(counters.values())
print(f'total={total} sum={s}')
print('ok' if total == s else 'mismatch')
PY
    if [ $rc -eq 0 ] && grep -q '^ok$' "$out/check.txt"; then
        echo "  PASS $name ($(head -1 $out/check.txt))"
        PASS=$((PASS+1))
    else
        echo "  FAIL $name (exit=$rc, $(cat $out/check.txt | tr '\n' ' '))"
        sed 's/^/      /' "$out/stderr"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name")
    fi
}

# G7: concurrency stress with all filters on; -t 1 vs -t 8 byte-identical
run_concurrency_with_filters() {
    local name="concurrency_with_filters"
    local out="$OUT_DIR/$name"
    mkdir -p "$out"
    python3 - <<PY > /dev/null
import random
random.seed(7)
B='ACGT'; rc=str.maketrans('ACGT','TGCA')
swap={'A':'T','C':'G','G':'C','T':'A'}
with open('$out/R1.fq','w') as r1, open('$out/R2.fq','w') as r2:
    for i in range(2000):
        f=''.join(random.choice(B) for _ in range(220))
        a=list(f[:150]); b=list(f[-150:].translate(rc)[::-1])
        qa=['?']*150; qb=['?']*150
        for j in range(150):
            if random.random() < 0.01: a[j]=swap[a[j]]; qa[j]=chr(33+10)
            if random.random() < 0.01: b[j]=swap[b[j]]; qb[j]=chr(33+10)
        r1.write(f"@r{i}/1\n{''.join(a)}\n+\n{''.join(qa)}\n")
        r2.write(f"@r{i}/2\n{''.join(b)}\n+\n{''.join(qb)}\n")
PY
    local ad="AGATCGGAAGAGCACACGTCT"
    "$PBJ_BIN" -f "$out/R1.fq" -r "$out/R2.fq" -o "$out/t1" -t 1 \
        -a "$ad" --max-ee 1.0 --max-mono-frac 0.9 --qtrim 10 2>/dev/null
    "$PBJ_BIN" -f "$out/R1.fq" -r "$out/R2.fq" -o "$out/t8" -t 8 -b 64 \
        -a "$ad" --max-ee 1.0 --max-mono-frac 0.9 --qtrim 10 2>/dev/null
    if diff -q "$out/t1.assembled.fastq" "$out/t8.assembled.fastq" >/dev/null \
       && diff -q "$out/t1.unassembled.forward.fastq" "$out/t8.unassembled.forward.fastq" >/dev/null \
       && diff -q "$out/t1.unassembled.reverse.fastq" "$out/t8.unassembled.reverse.fastq" >/dev/null; then
        echo "  PASS $name (filters + threading deterministic across all 3 outputs)"
        PASS=$((PASS+1))
    else
        echo "  FAIL $name (outputs differ)"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name")
    fi
}

# G8: verbose-mode smoke test, pins per-pair stderr format
run_verbose_smoke() {
    local name="verbose_smoke"
    local fix="$FIX_DIR/01_perfect_overlap"
    local out="$OUT_DIR/$name"
    mkdir -p "$out"
    "$PBJ_BIN" -f "$fix/R1.fq" -r "$fix/R2.fq" -o "$out/o" -v \
        2>"$out/stderr" >"$out/stdout"
    if grep -q 'OK s=' "$out/stderr" && grep -q 'name=perfect' "$out/stderr"; then
        echo "  PASS $name (per-pair format intact)"
        PASS=$((PASS+1))
    else
        echo "  FAIL $name"
        sed 's/^/      /' "$out/stderr"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name")
    fi
}

# round-4 R1: lowercase adapter strings must also work (real bug regression)
run_adapter_lowercase() {
    local name="adapter_lowercase"
    local fix="$FIX_DIR/17_adapter_trim"
    local out="$OUT_DIR/$name"
    mkdir -p "$out"
    "$PBJ_BIN" -f "$fix/R1.fq" -r "$fix/R2.fq" -o "$out/o" \
        -a "agatcggaagagcacacgtct" 2>"$out/stderr" >"$out/stdout"
    local rc=$?
    local got_merged
    got_merged=$(count_records "$out/o.assembled.fastq")
    if [ $rc -eq 0 ] && [ "$got_merged" = "1" ]; then
        echo "  PASS $name (lowercase -a works)"
        PASS=$((PASS+1))
    else
        echo "  FAIL $name (exit=$rc merged=$got_merged)"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name")
    fi
}

# round-4 R2: stdin input via -f -
run_stdin_input() {
    local name="stdin_input"
    local fix="$FIX_DIR/01_perfect_overlap"
    local out="$OUT_DIR/$name"
    mkdir -p "$out"
    "$PBJ_BIN" -f - -r "$fix/R2.fq" -o "$out/o" < "$fix/R1.fq" \
        2>"$out/stderr" >"$out/stdout"
    local rc=$?
    local got_merged
    got_merged=$(count_records "$out/o.assembled.fastq")
    if [ $rc -eq 0 ] && [ "$got_merged" = "1" ]; then
        echo "  PASS $name (-f - reads R1 from stdin)"
        PASS=$((PASS+1))
    else
        echo "  FAIL $name (exit=$rc merged=$got_merged)"
        sed 's/^/      /' "$out/stderr"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name")
    fi
}

# round-4 r3: max-phred input (phred 93 = ascii 126 = '~')
run_max_phred_input() {
    local name="max_phred_input"
    local out="$OUT_DIR/$name"
    mkdir -p "$out"
    python3 - <<PY > /dev/null
import random
random.seed(91)
B='ACGT'; rc=str.maketrans('ACGT','TGCA')
frag=''.join(random.choice(B) for _ in range(150))
r1=frag[:100]; r2=frag[-100:].translate(rc)[::-1]
q='~'*100  # ASCII 126 = Phred 93
open('$out/R1.fq','w').write(f'@m93/1\n{r1}\n+\n{q}\n')
open('$out/R2.fq','w').write(f'@m93/2\n{r2}\n+\n{q}\n')
PY
    "$PBJ_BIN" -f "$out/R1.fq" -r "$out/R2.fq" -o "$out/o" \
        2>"$out/stderr" >"$out/stdout"
    local rc=$?
    local got_merged
    got_merged=$(count_records "$out/o.assembled.fastq")
    if [ $rc -eq 0 ] && [ "$got_merged" = "1" ]; then
        echo "  PASS $name (Phred 93 input merges, no overflow)"
        PASS=$((PASS+1))
    else
        echo "  FAIL $name (exit=$rc merged=$got_merged)"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name")
    fi
}

# round-4 R4: asymmetric read lengths; pins both merged_len and merged_seq
run_asymmetric_reads() {
    local name="asymmetric_reads"
    local out="$OUT_DIR/$name"
    mkdir -p "$out"
    python3 - <<PY > "$out/expected_seq"
import random
random.seed(80)
B='ACGT'; rc=str.maketrans('ACGT','TGCA')
frag=''.join(random.choice(B) for _ in range(200))
r1=frag[:180]; r2=frag[-50:].translate(rc)[::-1]
open('$out/R1.fq','w').write(f'@asym/1\n{r1}\n+\n{"?"*180}\n')
open('$out/R2.fq','w').write(f'@asym/2\n{r2}\n+\n{"?"*50}\n')
print(frag, end='')
PY
    "$PBJ_BIN" -f "$out/R1.fq" -r "$out/R2.fq" -o "$out/o" -m 20 \
        2>"$out/stderr" >"$out/stdout"
    local rc=$?
    local got_merged got_len got_seq exp_seq
    got_merged=$(count_records "$out/o.assembled.fastq")
    got_len=$(awk 'NR==2{print length($0); exit}' "$out/o.assembled.fastq")
    got_seq=$(awk 'NR==2{print; exit}' "$out/o.assembled.fastq")
    exp_seq=$(cat "$out/expected_seq")
    local seq_ok=no
    [ "$got_seq" = "$exp_seq" ] && seq_ok=yes
    if [ $rc -eq 0 ] && [ "$got_merged" = "1" ] && [ "$got_len" = "200" ] && [ "$seq_ok" = "yes" ]; then
        echo "  PASS $name (R1=180, R2=50, merged_len=200, seq matches)"
        PASS=$((PASS+1))
    else
        echo "  FAIL $name (exit=$rc merged=$got_merged len=$got_len seq=$seq_ok)"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name")
    fi
}

# round-5 R1: FAIL_PVALUE branch is otherwise untested across the whole
# suite. construct a pair whose only overlap is exactly k bases long
# (n_free = 0 after seed adjustment) so p_val = 1.0 > 0.01 and the
# pair is rejected by the p-value gate rather than no_candidate or length.
run_fail_pvalue_branch() {
    local name="fail_pvalue_branch"
    local out="$OUT_DIR/$name"
    mkdir -p "$out"
    python3 - <<PY > /dev/null
import random
random.seed(50)
B='ACGT'; rc=str.maketrans('ACGT','TGCA')
# 36bp fragment, 24bp reads -> overlap = 12 = k
frag = ''.join(random.choice(B) for _ in range(36))
r1 = frag[:24]
r2 = frag[-24:].translate(rc)[::-1]
open('$out/R1.fq','w').write(f'@pv/1\n{r1}\n+\n{"?"*24}\n')
open('$out/R2.fq','w').write(f'@pv/2\n{r2}\n+\n{"?"*24}\n')
PY
    "$PBJ_BIN" -f "$out/R1.fq" -r "$out/R2.fq" -o "$out/o" \
        -k 12 -m 12 -l 12 2>"$out/stderr" >"$out/stdout"
    local rc=$?
    local fail_pv
    fail_pv=$(grep '^  fail p-value' "$out/stderr" | awk '{print $3}')
    if [ $rc -eq 0 ] && [ "$fail_pv" = "1" ]; then
        echo "  PASS $name (overlap exactly = k -> FAIL_PVALUE)"
        PASS=$((PASS+1))
    else
        echo "  FAIL $name (exit=$rc fail_pvalue=$fail_pv)"
        sed 's/^/      /' "$out/stderr"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name")
    fi
}

# round-5 R2: mixed-case adapter exercises the per-byte up() helper
run_adapter_mixed_case() {
    local name="adapter_mixed_case"
    local fix="$FIX_DIR/17_adapter_trim"
    local out="$OUT_DIR/$name"
    mkdir -p "$out"
    "$PBJ_BIN" -f "$fix/R1.fq" -r "$fix/R2.fq" -o "$out/o" \
        -a "AgAtCgGaAgAgCaCaCgTcT" 2>"$out/stderr" >"$out/stdout"
    local rc=$?
    local got_merged
    got_merged=$(count_records "$out/o.assembled.fastq")
    if [ $rc -eq 0 ] && [ "$got_merged" = "1" ]; then
        echo "  PASS $name (mixed-case adapter trims correctly)"
        PASS=$((PASS+1))
    else
        echo "  FAIL $name (exit=$rc merged=$got_merged)"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name")
    fi
}

# round-4 R5: -M max-overlap restricts the candidate range
run_max_overlap_restrict() {
    local name="max_overlap_restrict"
    local fix="$FIX_DIR/01_perfect_overlap"
    local out="$OUT_DIR/$name"
    mkdir -p "$out"
    "$PBJ_BIN" -f "$fix/R1.fq" -r "$fix/R2.fq" -o "$out/o" -M 30 \
        2>"$out/stderr" >"$out/stdout"
    local rc=$?
    local got_merged
    got_merged=$(count_records "$out/o.assembled.fastq")
    if [ $rc -eq 0 ] && [ "$got_merged" = "0" ]; then
        echo "  PASS $name (-M 30 rules out the real 50bp overlap)"
        PASS=$((PASS+1))
    else
        echo "  FAIL $name (exit=$rc merged=$got_merged)"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name")
    fi
}

# round-4 R6+R7: --version / --help exact stdout shape
run_version_help_exact() {
    local name="version_help_exact_stdout"
    local out="$OUT_DIR/$name"
    mkdir -p "$out"
    "$PBJ_BIN" --version > "$out/version.out" 2>/dev/null
    "$PBJ_BIN" --help    > "$out/help.out"    2>/dev/null
    if grep -q '^pbj ' "$out/version.out" \
       && grep -q '^usage: pbj' "$out/help.out" \
       && grep -q 'forward reads' "$out/help.out"; then
        echo "  PASS $name"
        PASS=$((PASS+1))
    else
        echo "  FAIL $name"
        head -3 "$out/version.out" "$out/help.out" 2>/dev/null | sed 's/^/      /'
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name")
    fi
}

# round-4 R8: empty adapter -a "" must be a no-op, not a crash
run_empty_adapter() {
    local name="empty_adapter"
    local fix="$FIX_DIR/01_perfect_overlap"
    local out="$OUT_DIR/$name"
    mkdir -p "$out"
    "$PBJ_BIN" -f "$fix/R1.fq" -r "$fix/R2.fq" -o "$out/o" -a "" \
        2>"$out/stderr" >"$out/stdout"
    local rc=$?
    local got_merged
    got_merged=$(count_records "$out/o.assembled.fastq")
    if [ $rc -eq 0 ] && [ "$got_merged" = "1" ]; then
        echo "  PASS $name (empty adapter is a clean no-op)"
        PASS=$((PASS+1))
    else
        echo "  FAIL $name (exit=$rc merged=$got_merged)"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name")
    fi
}

# round-4 R9: all-IUPAC (no ACGT, no N) read should be unmergeable
run_all_iupac_no_n() {
    local name="all_iupac_no_n"
    local out="$OUT_DIR/$name"
    mkdir -p "$out"
    python3 - <<PY > /dev/null
seq = ('YRSWKM' * 17)[:100]
q = '?' * 100
open('$out/R1.fq','w').write(f'@iu/1\n{seq}\n+\n{q}\n')
open('$out/R2.fq','w').write(f'@iu/2\n{seq}\n+\n{q}\n')
PY
    "$PBJ_BIN" -f "$out/R1.fq" -r "$out/R2.fq" -o "$out/o" \
        2>"$out/stderr" >"$out/stdout"
    local rc=$?
    local got_merged got_unassembled
    got_merged=$(count_records "$out/o.assembled.fastq")
    got_unassembled=$(count_records "$out/o.unassembled.forward.fastq")
    if [ $rc -eq 0 ] && [ "$got_merged" = "0" ] && [ "$got_unassembled" = "1" ]; then
        echo "  PASS $name (no kmers from IUPAC-only reads -> NO_CANDIDATE)"
        PASS=$((PASS+1))
    else
        echo "  FAIL $name (exit=$rc merged=$got_merged un=$got_unassembled)"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name")
    fi
}

# round-4 R10: mixed-rejection stats-sum invariant. construct a varied
# input that hits every failure branch, then assert the stats sum equals
# total_pairs. extends G6 from a partial check to all branches.
run_mixed_rejection_invariant() {
    local name="mixed_rejection_invariant"
    local out="$OUT_DIR/$name"
    mkdir -p "$out"
    python3 - <<PY > /dev/null
import random
random.seed(17)
B='ACGT'; rc=str.maketrans('ACGT','TGCA')
def revcomp(s): return s.translate(rc)[::-1]
swap={'A':'T','C':'A','G':'C','T':'G'}
pairs=[]
# 30 clean merges (MERGED_OK)
for _ in range(30):
    f=''.join(random.choice(B) for _ in range(160))
    pairs.append((f[:100], revcomp(f[-100:]), '?'*100, '?'*100))
# 30 with no real overlap (NO_CANDIDATE)
for _ in range(30):
    f1=''.join(random.choice(B) for _ in range(150))
    f2=''.join(random.choice(B) for _ in range(150))
    pairs.append((f1, f2, '?'*150, '?'*150))
# 30 with very low quality (FAIL_MIN_QUAL with -Q 20)
for _ in range(30):
    f=''.join(random.choice(B) for _ in range(160))
    pairs.append((f[:100], revcomp(f[-100:]), chr(33+3)*100, chr(33+3)*100))
# 30 with high-Q mismatches in overlap (FAIL_AMBIGUOUS under --strict)
for _ in range(30):
    f=''.join(random.choice(B) for _ in range(160))
    r1=list(f[:100])
    for pos in (60, 70, 80):
        r1[pos]=swap[r1[pos]]
    pairs.append((''.join(r1), revcomp(f[-100:]), '?'*100, '?'*100))
# 30 short merges that fail -l 200 (FAIL_LENGTH)
for _ in range(30):
    f=''.join(random.choice(B) for _ in range(120))
    pairs.append((f[:80], revcomp(f[-80:]), '?'*80, '?'*80))
# 30 all-N (would also be NO_CANDIDATE because no kmers)
for _ in range(30):
    pairs.append(('N'*100, 'N'*100, '?'*100, '?'*100))
with open('$out/R1.fq','w') as a, open('$out/R2.fq','w') as b:
    for i,(r1,r2,q1,q2) in enumerate(pairs):
        a.write(f'@m{i}/1\n{r1}\n+\n{q1}\n')
        b.write(f'@m{i}/2\n{r2}\n+\n{q2}\n')
PY
    "$PBJ_BIN" -f "$out/R1.fq" -r "$out/R2.fq" -o "$out/o" \
        -Q 20 -l 200 -s --max-ee 5.0 --max-mono-frac 0.95 \
        2>"$out/stderr" >"$out/stdout"
    local rc=$?
    python3 - <<PY > "$out/check.txt"
import re
total=None; counters={}
with open('$out/stderr') as f:
    for line in f:
        m=re.match(r'^  (total pairs|merged|no candidate|fail [\w -]+?)\s+(\d+)', line)
        if not m: continue
        k=m.group(1).strip(); v=int(m.group(2))
        if k=='total pairs': total=v
        else: counters[k]=v
s=sum(counters.values())
print(f'total={total} sum={s}')
print('ok' if total == s else 'mismatch')
PY
    if [ $rc -eq 0 ] && grep -q '^ok$' "$out/check.txt"; then
        echo "  PASS $name ($(head -1 $out/check.txt))"
        PASS=$((PASS+1))
    else
        echo "  FAIL $name (exit=$rc $(cat $out/check.txt | tr '\n' ' '))"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name")
    fi
}

# G9: adapter sequence containing N acts as a wildcard match
run_adapter_with_n() {
    local name="adapter_with_n_wildcard"
    local out="$OUT_DIR/$name"
    mkdir -p "$out"
    python3 - <<PY > /dev/null
import random
random.seed(31)
B='ACGT'; rc=str.maketrans('ACGT','TGCA')
frag=''.join(random.choice(B) for _ in range(80))
real='AGATCGGAAGAGCACACGTCT'
r1 = frag + real
r2 = frag[::-1].translate(rc) + real
q = '?' * len(r1)
open('$out/R1.fq','w').write(f'@nwc/1\n{r1}\n+\n{q}\n')
open('$out/R2.fq','w').write(f'@nwc/2\n{r2}\n+\n{q}\n')
PY
    # adapter passed has N at position 2 (was 'A')
    local masked="AGNTCGGAAGAGCACACGTCT"
    "$PBJ_BIN" -f "$out/R1.fq" -r "$out/R2.fq" -o "$out/o" \
        -a "$masked" 2>"$out/stderr" >"$out/stdout"
    local rc=$?
    local trim_r1
    trim_r1=$(grep -E '^  adapter trimmed R1' "$out/stderr" | awk '{print $4}')
    if [ $rc -eq 0 ] && [ "$trim_r1" = "1" ]; then
        echo "  PASS $name (N in adapter matches wildcard, trim=1)"
        PASS=$((PASS+1))
    else
        echo "  FAIL $name (trim=$trim_r1, exit=$rc)"
        sed 's/^/      /' "$out/stderr"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name")
    fi
}

# gzip-level round-trip: -1 and -9 produce decompressible output that is
# byte-identical when decompressed.
run_gzip_level_roundtrip() {
    local name="gzip_level_roundtrip"
    local out="$OUT_DIR/$name"
    mkdir -p "$out"
    local r1="$FIX_DIR/01_perfect_overlap/R1.fq"
    local r2="$FIX_DIR/01_perfect_overlap/R2.fq"
    "$PBJ_BIN" -f "$r1" -r "$r2" -o "$out/lvl1" -z --gzip-level 1 2>/dev/null
    "$PBJ_BIN" -f "$r1" -r "$r2" -o "$out/lvl9" -z --gzip-level 9 2>/dev/null
    gunzip -c "$out/lvl1.assembled.fastq.gz" > "$out/lvl1.assembled.fastq"
    gunzip -c "$out/lvl9.assembled.fastq.gz" > "$out/lvl9.assembled.fastq"
    if diff -q "$out/lvl1.assembled.fastq" "$out/lvl9.assembled.fastq" >/dev/null; then
        echo "  PASS $name (-1 and -9 decompress to identical content)"
        PASS=$((PASS+1))
    else
        echo "  FAIL $name"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name")
    fi
}

# mismatched R1/R2 record counts
run_count_mismatch() {
    local name="count_mismatch"
    local out="$OUT_DIR/$name"
    mkdir -p "$out"
    printf "@r/1\nACGT\n+\n!!!!\n@s/1\nACGT\n+\n!!!!\n" > "$out/R1.fq"
    printf "@r/2\nACGT\n+\n!!!!\n" > "$out/R2.fq"
    "$PBJ_BIN" -f "$out/R1.fq" -r "$out/R2.fq" -o "$out/o" 2>"$out/stderr" >"$out/stdout"
    local rc=$?
    if [ $rc -eq 2 ] && grep -q "R2 ended before R1" "$out/stderr"; then
        echo "  PASS $name (R2 ended before R1)"
        PASS=$((PASS+1))
    else
        echo "  FAIL $name (exit=$rc)"
        sed 's/^/      /' "$out/stderr"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name")
    fi
}

# thread determinism: same output -t 1 vs -t 4
run_thread_determinism() {
    local name="thread_determinism"
    local out="$OUT_DIR/$name"
    mkdir -p "$out"
    python3 - <<PY > /dev/null
import random, os
random.seed(3)
B='ACGT'; rc=str.maketrans('ACGT','TGCA')
with open('$out/R1.fq','w') as a, open('$out/R2.fq','w') as b:
    for i in range(500):
        f=''.join(random.choice(B) for _ in range(160))
        a.write(f'@r{i}/1\n{f[:100]}\n+\n{"?"*100}\n')
        b.write(f'@r{i}/2\n{f[-100:].translate(rc)[::-1]}\n+\n{"?"*100}\n')
PY
    "$PBJ_BIN" -f "$out/R1.fq" -r "$out/R2.fq" -o "$out/t1" -t 1     2>/dev/null
    "$PBJ_BIN" -f "$out/R1.fq" -r "$out/R2.fq" -o "$out/t4" -t 4 -b 16 2>/dev/null
    if diff -q "$out/t1.assembled.fastq" "$out/t4.assembled.fastq" >/dev/null; then
        echo "  PASS $name (-t 1 and -t 4 byte-identical)"
        PASS=$((PASS+1))
    else
        echo "  FAIL $name (-t 1 vs -t 4 differ)"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name")
    fi
}

# CLI validation tests
expect_exit() {
    local name="$1"; local expect_rc="$2"; shift 2
    local out="$OUT_DIR/$name"
    mkdir -p "$out"
    "$PBJ_BIN" "$@" 2>"$out/stderr" >"$out/stdout"
    local rc=$?
    if [ "$rc" = "$expect_rc" ]; then
        echo "  PASS $name (exit $rc)"
        PASS=$((PASS+1))
    else
        echo "  FAIL $name (exit=$rc, expected $expect_rc)"
        sed 's/^/      /' "$out/stderr"
        FAIL=$((FAIL+1)); FAILED_NAMES+=("$name")
    fi
}

run_cli_tests() {
    local r1="$FIX_DIR/01_perfect_overlap/R1.fq"
    local r2="$FIX_DIR/01_perfect_overlap/R2.fq"
    expect_exit "cli_missing_out"     1 -f "$r1" -r "$r2"
    expect_exit "cli_missing_f"       1 -r "$r2" -o /tmp/pbj_cli
    expect_exit "cli_missing_file"    2 -f /no/such/file -r "$r2" -o /tmp/pbj_cli
    expect_exit "cli_bad_kmer_low"    1 -f "$r1" -r "$r2" -o /tmp/pbj_cli -k 3
    expect_exit "cli_bad_kmer_high"   1 -f "$r1" -r "$r2" -o /tmp/pbj_cli -k 32
    expect_exit "cli_bad_pvalue_neg"  1 -f "$r1" -r "$r2" -o /tmp/pbj_cli -p -0.1
    expect_exit "cli_bad_pvalue_big"  1 -f "$r1" -r "$r2" -o /tmp/pbj_cli -p 1.5
    expect_exit "cli_bad_pmatch_zero" 1 -f "$r1" -r "$r2" -o /tmp/pbj_cli --p-match 0
    expect_exit "cli_bad_gzlevel"     1 -f "$r1" -r "$r2" -o /tmp/pbj_cli --gzip-level 0
    expect_exit "cli_bad_maxee_neg"   1 -f "$r1" -r "$r2" -o /tmp/pbj_cli --max-ee -1
    expect_exit "cli_bad_mono"        1 -f "$r1" -r "$r2" -o /tmp/pbj_cli --max-mono-frac 2
    expect_exit "cli_interleaved_with_r2" 1 -f "$r1" -r "$r2" -o /tmp/pbj_cli --interleaved-in
    expect_exit "cli_version"         0 --version
    expect_exit "cli_help"            0 --help
}

echo "running fixture tests"
for d in "$FIX_DIR"/*/; do
    run_fixture "$(basename "$d")"
done
run_negative_dovetail
run_pair_mismatch
run_pair_same_suffix
run_count_mismatch
run_thread_determinism
run_interleaved_odd
run_interleaved_r2_first
run_adapter_no_false_positive
run_gzip_level_roundtrip
run_separate_file_swap
run_no_bonferroni_smoke
run_qtrim_zero_minlen_zero
run_fixture07_no_strict
run_adapter_min_match_override
run_stats_sum_invariant
run_concurrency_with_filters
run_verbose_smoke
run_adapter_with_n
run_adapter_lowercase
run_stdin_input
run_max_phred_input
run_asymmetric_reads
run_max_overlap_restrict
run_version_help_exact
run_empty_adapter
run_all_iupac_no_n
run_mixed_rejection_invariant
run_fail_pvalue_branch
run_adapter_mixed_case

echo "running CLI validation tests"
run_cli_tests

echo ""
echo "summary: $PASS passed, $FAIL failed"
if [ $FAIL -gt 0 ]; then
    echo "failed: ${FAILED_NAMES[*]}"
    exit 1
fi
exit 0
