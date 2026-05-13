#!/usr/bin/env python3
"""generate pbj test fixtures into tests/fixtures/.

each fixture has R1.fq (or R1.fq.gz), R2.fq, and EXPECTED.txt
(key=value lines). deterministic given the random seed below.
"""

from __future__ import annotations

import gzip
import os
import random
import sys

SEED: int = 42
random.seed(SEED)

BASES: str = "ACGT"
RC = str.maketrans("ACGTacgt", "TGCAtgca")
HIGH_Q: str = chr(33 + 30)  # phred 30
MID_Q:  str = chr(33 + 20)
LOW_Q:  str = chr(33 + 5)   # phred  5
VLOW_Q: str = chr(33 + 2)   # phred  2


def rand_seq(n: int) -> str:
    return "".join(random.choice(BASES) for _ in range(n))


def revcomp(s: str) -> str:
    return s.translate(RC)[::-1]


def write_fq(path: str, name: str, seq: str, qual: str) -> None:
    with open(path, "w") as f:
        f.write(f"@{name}\n{seq}\n+\n{qual}\n")


def write_fq_gz(path: str, name: str, seq: str, qual: str) -> None:
    with gzip.open(path, "wt") as f:
        f.write(f"@{name}\n{seq}\n+\n{qual}\n")


def write_expected(path: str, **kwargs: object) -> None:
    with open(path, "w") as f:
        for k, v in kwargs.items():
            f.write(f"{k}={v}\n")


def fixture(out_dir: str, name: str) -> str:
    d = os.path.join(out_dir, name)
    os.makedirs(d, exist_ok=True)
    return d


def make_paired(
    out_dir: str,
    fragment: str,
    r1_len: int,
    r2_len: int,
    name: str = "pair",
    q_r1: str | None = None,
    q_r2: str | None = None,
) -> None:
    r1 = fragment[:r1_len]
    r2 = revcomp(fragment[-r2_len:])
    q1 = q_r1 if q_r1 is not None else HIGH_Q * len(r1)
    q2 = q_r2 if q_r2 is not None else HIGH_Q * len(r2)
    write_fq(os.path.join(out_dir, "R1.fq"), name + "/1", r1, q1)
    write_fq(os.path.join(out_dir, "R2.fq"), name + "/2", r2, q2)


def main() -> None:
    out_base: str = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "fixtures"
    )
    out_base = os.path.abspath(out_base)
    os.makedirs(out_base, exist_ok=True)

    # 01: perfect overlap
    d = fixture(out_base, "01_perfect_overlap")
    frag = rand_seq(150)
    make_paired(d, frag, 100, 100, name="perfect")
    write_expected(os.path.join(d, "EXPECTED.txt"),
                   merged=1, unassembled=0, merged_len=150,
                   merged_seq=frag)

    # 02: no overlap
    d = fixture(out_base, "02_no_overlap")
    frag = rand_seq(400)
    make_paired(d, frag, 100, 100, name="nooverlap")
    write_expected(os.path.join(d, "EXPECTED.txt"),
                   merged=0, unassembled=1)

    # 03: short overlap
    d = fixture(out_base, "03_short_overlap")
    frag = rand_seq(175)
    make_paired(d, frag, 100, 100, name="short")
    write_expected(os.path.join(d, "EXPECTED.txt"),
                   merged=1, unassembled=0, merged_len=175,
                   merged_seq=frag)

    # 04: full read overlap
    d = fixture(out_base, "04_full_overlap")
    frag = rand_seq(100)
    make_paired(d, frag, 100, 100, name="full")
    write_expected(os.path.join(d, "EXPECTED.txt"),
                   merged=1, unassembled=0, merged_len=100,
                   merged_seq=frag)

    # 05: mismatches with asymmetric quality
    d = fixture(out_base, "05_mismatches")
    frag = rand_seq(150)
    r1 = list(frag[:100])
    q1 = list(HIGH_Q * 100)
    swap: dict[str, str] = {"A":"T","C":"A","G":"C","T":"G"}
    for pos in (60, 70, 80):
        r1[pos] = swap[r1[pos]]
        q1[pos] = LOW_Q
    r1_seq = "".join(r1)
    r1_qual = "".join(q1)
    r2_seq = revcomp(frag[-100:])
    r2_qual = HIGH_Q * 100
    write_fq(os.path.join(d, "R1.fq"), "mismatch/1", r1_seq, r1_qual)
    write_fq(os.path.join(d, "R2.fq"), "mismatch/2", r2_seq, r2_qual)
    write_expected(os.path.join(d, "EXPECTED.txt"),
                   merged=1, unassembled=0, merged_len=150,
                   merged_seq=frag)

    # 06: adapter readthrough / dovetail
    d = fixture(out_base, "06_dovetail")
    frag = rand_seq(80)
    r1_seq = frag + rand_seq(20)
    r2_seq = revcomp(frag) + rand_seq(20)
    write_fq(os.path.join(d, "R1.fq"), "dovetail/1", r1_seq, HIGH_Q * 100)
    write_fq(os.path.join(d, "R2.fq"), "dovetail/2", r2_seq, HIGH_Q * 100)
    write_expected(os.path.join(d, "EXPECTED.txt"),
                   merged=1, unassembled=0,
                   merged_len=80, merged_seq=frag,
                   pbj_flags="--dovetail")

    # 07: strict-mode reject. high-Q mismatches on BOTH sides -> ambiguous.
    # without --strict the pair merges (mismatches are picked by score, but
    # since both reads are high-Q on the disagreement, ambiguity is flagged).
    # we run this fixture with --strict; without --strict it would merge.
    d = fixture(out_base, "07_strict_reject")
    frag = rand_seq(150)
    r1 = list(frag[:100])
    r2_tmpl = list(frag[-100:])
    # introduce 3 mismatches in the overlap region with HIGH quality on both
    swap = {"A":"T","C":"A","G":"C","T":"G"}
    for pos in (60, 70, 80):
        r1[pos] = swap[r1[pos]]  # R1 disagrees from fragment
    r1_seq = "".join(r1)
    r2_seq = revcomp("".join(r2_tmpl))
    write_fq(os.path.join(d, "R1.fq"), "strict/1", r1_seq, HIGH_Q * 100)
    write_fq(os.path.join(d, "R2.fq"), "strict/2", r2_seq, HIGH_Q * 100)
    write_expected(os.path.join(d, "EXPECTED.txt"),
                   merged=0, unassembled=1,
                   pbj_flags="--strict")

    # 08: min-mean-qual reject. valid overlap but every base at phred 2.
    # without -Q the pair merges. with -Q 20 it fails the mean-qual check.
    d = fixture(out_base, "08_min_mean_qual")
    frag = rand_seq(150)
    make_paired(d, frag, 100, 100, name="lowq",
                q_r1=VLOW_Q * 100, q_r2=VLOW_Q * 100)
    write_expected(os.path.join(d, "EXPECTED.txt"),
                   merged=0, unassembled=1,
                   pbj_flags="-Q 20")

    # 09: min-length reject. valid merge of 150bp, rejected by -l 200.
    d = fixture(out_base, "09_min_length")
    frag = rand_seq(150)
    make_paired(d, frag, 100, 100, name="minlen")
    write_expected(os.path.join(d, "EXPECTED.txt"),
                   merged=0, unassembled=1,
                   pbj_flags="-l 200")

    # 10: max-length reject. valid merge of 150bp, rejected by -L 100.
    d = fixture(out_base, "10_max_length")
    frag = rand_seq(150)
    make_paired(d, frag, 100, 100, name="maxlen")
    write_expected(os.path.join(d, "EXPECTED.txt"),
                   merged=0, unassembled=1,
                   pbj_flags="-L 100")

    # 11: IUPAC ambiguity codes in input. one 'Y' in R1's overlap region.
    # kmers around Y still seed (we treat non-ACGT as invalid kmer chars but
    # kmers elsewhere in the read are fine). overlap scoring treats Y as a
    # real base (not N), so Y vs the correct base counts as a mismatch.
    d = fixture(out_base, "11_iupac")
    frag = rand_seq(150)
    r1 = list(frag[:100])
    r1[55] = "Y"  # ambiguity in overlap region
    r1_seq = "".join(r1)
    r2_seq = revcomp(frag[-100:])
    write_fq(os.path.join(d, "R1.fq"), "iupac/1", r1_seq, HIGH_Q * 100)
    write_fq(os.path.join(d, "R2.fq"), "iupac/2", r2_seq, HIGH_Q * 100)
    write_expected(os.path.join(d, "EXPECTED.txt"),
                   merged=1, unassembled=0)

    # 12: all-N reads. no valid kmers -> no candidates -> no merge.
    d = fixture(out_base, "12_all_n")
    write_fq(os.path.join(d, "R1.fq"), "alln/1", "N" * 100, HIGH_Q * 100)
    write_fq(os.path.join(d, "R2.fq"), "alln/2", "N" * 100, HIGH_Q * 100)
    write_expected(os.path.join(d, "EXPECTED.txt"),
                   merged=0, unassembled=1)

    # 13: reads shorter than k. n1 < k short-circuits the kmer pass.
    d = fixture(out_base, "13_short_reads")
    write_fq(os.path.join(d, "R1.fq"), "short/1", "ACGTACGT", HIGH_Q * 8)
    write_fq(os.path.join(d, "R2.fq"), "short/2", "TGCATGCA", HIGH_Q * 8)
    write_expected(os.path.join(d, "EXPECTED.txt"),
                   merged=0, unassembled=1)

    # 14: gzipped input. same content as fixture 01, but .fq.gz.
    d = fixture(out_base, "14_gzip_input")
    frag = rand_seq(150)
    r1 = frag[:100]
    r2 = revcomp(frag[-100:])
    write_fq_gz(os.path.join(d, "R1.fq.gz"), "gzipped/1", r1, HIGH_Q * 100)
    write_fq_gz(os.path.join(d, "R2.fq.gz"), "gzipped/2", r2, HIGH_Q * 100)
    write_expected(os.path.join(d, "EXPECTED.txt"),
                   merged=1, unassembled=0,
                   merged_len=150, merged_seq=frag,
                   input_ext=".fq.gz")

    # 15: casava 1.8+ header. tests that the merged read's name keeps only
    # the core id (everything before the space) and strips the directional
    # flag carried in the comment field.
    d = fixture(out_base, "15_casava_header")
    frag = rand_seq(150)
    core = "HISEQ:123:FC:1:1101:1234:5678"
    write_fq(os.path.join(d, "R1.fq"),
             f"{core} 1:N:0:ATCG", frag[:100], HIGH_Q * 100)
    write_fq(os.path.join(d, "R2.fq"),
             f"{core} 2:N:0:ATCG", revcomp(frag[-100:]), HIGH_Q * 100)
    write_expected(os.path.join(d, "EXPECTED.txt"),
                   merged=1, unassembled=0,
                   merged_name=core,
                   merged_seq=frag)

    # 16: quality-offset 64. same merge as fixture 01 but qualities are
    # encoded at phred+64. output should still be phred+33.
    d = fixture(out_base, "16_quality_offset_64")
    frag = rand_seq(150)
    q64 = chr(64 + 30)  # phred 30 under offset 64
    write_fq(os.path.join(d, "R1.fq"), "q64/1", frag[:100], q64 * 100)
    write_fq(os.path.join(d, "R2.fq"), "q64/2", revcomp(frag[-100:]), q64 * 100)
    write_expected(os.path.join(d, "EXPECTED.txt"),
                   merged=1, unassembled=0,
                   merged_len=150, merged_seq=frag,
                   pbj_flags="-q 64",
                   merged_qual_first_byte=str(33 + 30))  # 30 + 33 = ASCII 63 '?'

    # 17: adapter trim. R1 = fragment + adapter, R2 = revcomp(fragment) + adapter.
    # without -a the reads dovetail and refuse to merge (no --dovetail).
    # with -a the adapter is trimmed off both reads, then they merge cleanly.
    d = fixture(out_base, "17_adapter_trim")
    frag = rand_seq(80)
    adapter = "AGATCGGAAGAGCACACGTCT"
    r1_seq = frag + adapter
    r2_seq = revcomp(frag) + adapter
    write_fq(os.path.join(d, "R1.fq"), "adapter/1", r1_seq, HIGH_Q * len(r1_seq))
    write_fq(os.path.join(d, "R2.fq"), "adapter/2", r2_seq, HIGH_Q * len(r2_seq))
    write_expected(os.path.join(d, "EXPECTED.txt"),
                   merged=1, unassembled=0,
                   merged_len=80, merged_seq=frag,
                   pbj_flags=f"-a {adapter}")

    # 18: max-ee reject. valid merge but every base at phred 10 -> ee = 0.1 *
    # 150 = 15 expected errors; rejected by --max-ee 1.
    d = fixture(out_base, "18_max_ee")
    frag = rand_seq(150)
    q_ee = chr(33 + 10) * 100
    make_paired(d, frag, 100, 100, name="ee", q_r1=q_ee, q_r2=q_ee)
    write_expected(os.path.join(d, "EXPECTED.txt"),
                   merged=0, unassembled=1,
                   pbj_flags="--max-ee 1")

    # 19: low-complexity reject. fragment is ~90% A. merged passes the merge
    # logic (high score from many matches) but the mono-base fraction exceeds
    # --max-mono-frac 0.85, so the merge is rejected.
    d = fixture(out_base, "19_low_complexity")
    random.seed(99)
    frag = "".join("A" if random.random() < 0.9 else random.choice("CGT") for _ in range(150))
    random.seed(SEED)  # restore so later fixtures stay deterministic
    make_paired(d, frag, 100, 100, name="mono")
    write_expected(os.path.join(d, "EXPECTED.txt"),
                   merged=0, unassembled=1,
                   pbj_flags="--max-mono-frac 0.85")

    # 20: qtrim. last 30 bases of both reads at phred 5; --qtrim 20 trims the
    # low-Q tail off the merged output. merged_len drops from 150 to 120 (the
    # overlap region has high Q from both reads, only the suffix R2_rc tail
    # is low). actually: the merged record is r1_prefix(50) + overlap(50) +
    # r2_rc_suffix(50). only R2_rc's suffix tail is at low Q in the merged
    # output. so qtrim trims that 50bp tail -> merged_len = 100.
    d = fixture(out_base, "20_qtrim")
    frag = rand_seq(150)
    q_r1 = HIGH_Q * 100  # R1 all high
    # r2: first 50 high, last 50 low. after rc of r2, the low region becomes the
    # leading 50bp of r2_rc. so r2_rc's suffix (after overlap) is positions
    # 50..100 in r2_rc, originally r2[0..50] (high), so high.
    # r2 reads from 3' end of fragment, so r2 = revcomp(fragment[-100:]).
    # r2 bytes [0..50] correspond to fragment[100..150] reversed/complemented;
    # r2 bytes [50..100] correspond to fragment[50..100] reversed/complemented.
    # after rc: r2_rc bytes [0..50] = fragment[50..100], r2_rc [50..100] = fragment[100..150].
    # so r2_rc's tail (bytes 50..100) = the suffix region of the merged read.
    # we want low quality on r2_rc's tail, i.e., low q on r2's head (bytes 0..50).
    q_r2 = (chr(33 + 5) * 50) + (HIGH_Q * 50)
    make_paired(d, frag, 100, 100, name="qtrim", q_r1=q_r1, q_r2=q_r2)
    # after qtrim 20 the trailing 50bp (low Q) is dropped, leaving 100bp
    write_expected(os.path.join(d, "EXPECTED.txt"),
                   merged=1, unassembled=0,
                   merged_len=100,
                   pbj_flags="--qtrim 20")

    # 21: interleaved input. one file with R1, R2, R1, R2 interleaved.
    d = fixture(out_base, "21_interleaved")
    n_pairs = 3
    with open(os.path.join(d, "INTERLEAVED.fq"), "w") as f:
        for i in range(n_pairs):
            frag = rand_seq(150)
            r1 = frag[:100]
            r2 = revcomp(frag[-100:])
            f.write(f"@inter{i}/1\n{r1}\n+\n{HIGH_Q*100}\n")
            f.write(f"@inter{i}/2\n{r2}\n+\n{HIGH_Q*100}\n")
    write_expected(os.path.join(d, "EXPECTED.txt"),
                   merged=n_pairs, unassembled=0,
                   input_file="INTERLEAVED.fq",
                   pbj_flags="--interleaved-in")

    # 22: adapter non-3' position. real adapter sequence placed at position 30
    # of a 100bp read but NOT extending to the read end. with 3'-anchored
    # trimming the adapter should NOT be detected and the read should merge
    # normally. without anchoring (the old algorithm) this would over-trim.
    d = fixture(out_base, "22_adapter_non_3prime")
    adapter = "AGATCGGAAGAGCACACGTCT"
    # build R1 with adapter at pos 30 surrounded by non-adapter sequence
    # construct R1 such that R1's 3' end (last 21 bases) does NOT match the
    # adapter prefix so the 3'-anchored algorithm can't trim there either.
    frag_left = "TGCATGCATGCATGCATGCATGCATGCATGCA"   # 32bp non-adapter
    tail = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"  # 36bp polyC tail
    r1_seq = (frag_left[:30] + adapter + tail)[:100]
    # mirrored R2
    r2_seq = revcomp(r1_seq[:80])  # 80bp overlap with R1's first 80bp
    # need a 100bp R2 that overlaps R1[0:80] (i.e., R2_rc = R1[0:80] + suffix)
    suffix = "GGGGGGGGGGGGGGGGGGGG"  # 20bp polyG
    r2_seq = revcomp(r1_seq[:80] + suffix)
    write_fq(os.path.join(d, "R1.fq"), "noadapter/1", r1_seq, HIGH_Q * len(r1_seq))
    write_fq(os.path.join(d, "R2.fq"), "noadapter/2", r2_seq, HIGH_Q * len(r2_seq))
    write_expected(os.path.join(d, "EXPECTED.txt"),
                   merged=1, unassembled=0,
                   pbj_flags=f"-a {adapter}")

    # 23: qtrim drops merged below min_assembly_len. R2's head is low-Q over
    # its first 80 bases, which (after RC) lands at the tail of R2_rc and
    # therefore at the tail of the merged read. qtrim 20 strips that 70bp
    # tail down to merged_len=100. with -l 110 the post-trim check rejects.
    d = fixture(out_base, "23_qtrim_below_min")
    frag = rand_seq(170)
    q_r2 = (chr(33 + 5) * 80) + (HIGH_Q * 20)
    make_paired(d, frag, 100, 100, name="qshort",
                q_r1=HIGH_Q * 100, q_r2=q_r2)
    write_expected(os.path.join(d, "EXPECTED.txt"),
                   merged=0, unassembled=1,
                   pbj_flags="--qtrim 20 -l 110")

    # 24: max-ee acceptance. fixture 01 has phred 30 throughout merged 150bp,
    # so EE = 150 * 10^-3 = 0.15. --max-ee 1.0 accepts, --max-ee 0.1 rejects
    # (we test the acceptance side here; rejection is fixture 18).
    d = fixture(out_base, "24_max_ee_accept")
    frag = rand_seq(150)
    make_paired(d, frag, 100, 100, name="eeok")
    write_expected(os.path.join(d, "EXPECTED.txt"),
                   merged=1, unassembled=0,
                   pbj_flags="--max-ee 1.0")

    # 25: empty input. zero-byte R1 and R2.
    d = fixture(out_base, "25_empty_input")
    open(os.path.join(d, "R1.fq"), "w").close()
    open(os.path.join(d, "R2.fq"), "w").close()
    write_expected(os.path.join(d, "EXPECTED.txt"),
                   merged=0, unassembled=0)

    # 26: FASTA input (no quality lines). kseq treats this as FASTA and
    # pbj_fastq_read synthesizes phred 30 for every base.
    d = fixture(out_base, "26_fasta_input")
    frag = rand_seq(150)
    with open(os.path.join(d, "R1.fq"), "w") as f:
        f.write(f">fasta/1\n{frag[:100]}\n")
    with open(os.path.join(d, "R2.fq"), "w") as f:
        f.write(f">fasta/2\n{revcomp(frag[-100:])}\n")
    write_expected(os.path.join(d, "EXPECTED.txt"),
                   merged=1, unassembled=0,
                   merged_len=150, merged_seq=frag)

    # 27: long reads. 5kbp R1 and R2 with 2kbp overlap exercise the kmer
    # table grow path and the candidate bitset grow path in overlap.c.
    d = fixture(out_base, "27_long_reads")
    frag = rand_seq(8000)
    make_paired(d, frag, 5000, 5000, name="long")
    write_expected(os.path.join(d, "EXPECTED.txt"),
                   merged=1, unassembled=0,
                   merged_len=8000, merged_seq=frag)

    fixtures_count = 27
    print(f"wrote {fixtures_count} fixtures to {out_base}")


if __name__ == "__main__":
    main()
