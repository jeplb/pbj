# pbj

A fast paired-end FASTQ read merger in C. Drop-in alternative to PEAR and
BBmerge with a single static binary and no runtime dependencies beyond
`zlib`, `pthread`, and libc. MIT licensed.

The name stands for **P**EAR + **B**Bmerge **J**oined.

## What it does

Takes Illumina-style paired-end FASTQ inputs (R1 and R2) and combines
pairs whose inserts are shorter than `len(R1) + len(R2)` into a single
longer consensus read. The non-overlapping pairs pass through unchanged
into separate "unassembled" outputs.

Typical merge rates on real Illumina 2x150 with ~250bp inserts run
99%+ at default settings. The merger is fork-join parallel; throughput
scales near-linearly to ~8 threads.

## Install

```
git clone <repo-url>
cd pbj
make
```

Build deps are `cc` (gcc or clang), `make`, `zlib`, `pthread`. Optional:
`python3` for regenerating test fixtures.

Available make targets:

```
make            # release build (-O3)
make native     # release + -march=native
make debug      # AddressSanitizer + UBSan
make test       # run the full test suite (70 tests)
make bench      # synthetic-throughput microbenchmark
```

The build produces a single `./pbj` binary.

## Quick start

```
pbj -f R1.fq.gz -r R2.fq.gz -o sample
```

Produces three files under the `-o` prefix:

- `sample.assembled.fastq` &mdash; merged reads
- `sample.unassembled.forward.fastq` &mdash; R1 of pairs that didn't merge
- `sample.unassembled.reverse.fastq` &mdash; R2 of those pairs

Add `-z` to gzip-compress all three outputs. Input is read transparently
whether `.fq` or `.fq.gz`. Output qualities are always written in
Phred+33 regardless of input offset.

## Common workflows

**Standard merge.** Random insert distribution, default settings.

```
pbj -f R1.fq.gz -r R2.fq.gz -o sample -t 4
```

**Trim adapter readthrough.** When inserts are shorter than read length,
sequencing reads past the fragment into the adapter on the 3' end. Pass
the adapter sequence to trim before merging.

```
pbj -f R1.fq.gz -r R2.fq.gz -o sample \
    -a AGATCGGAAGAGCACACGTCT \
    -A AGATCGGAAGAGCGTCGTGTA
```

If you pass only `-a`, `-A` defaults to the same sequence (fine for
TruSeq libraries where the first 13bp are identical).

**Dovetail / adapter readthrough without `-a`.** If the insert is
shorter than the read, the reads overshoot each other. `--dovetail`
lets pbj merge those pairs and drop the adapter tails.

```
pbj -f R1.fq.gz -r R2.fq.gz -o sample --dovetail
```

**Quality filtering.** Apply post-merge gates.

```
pbj -f R1.fq.gz -r R2.fq.gz -o sample \
    --max-ee 1.0          # USEARCH-style expected-error cap
    --max-mono-frac 0.85  # reject low-complexity merges
    --qtrim 20            # trim 3' bases below Phred 20
    -Q 25                 # reject if mean merged quality < 25
```

`--max-ee 1.0` is the USEARCH default for 16S amplicons. Leave off for
DNA-storage workflows that handle errors with downstream decoders.

**Interleaved input.** One file with alternating R1/R2 records.

```
pbj -f interleaved.fq.gz --interleaved-in -o sample
```

**Stdin.** Use `-` as the path for either `-f` or `-r`.

```
zcat R1.fq.gz | pbj -f - -r R2.fq.gz -o sample
```

**Strict mode.** Reject any merge with a high-quality disagreement
between R1 and R2.

```
pbj -f R1.fq.gz -r R2.fq.gz -o sample -s --strict-q 25
```

## All options

| flag | default | description |
| --- | --- | --- |
| `-f, --forward FILE` | required | forward reads (R1), `-` for stdin |
| `-r, --reverse FILE` | required | reverse reads (R2), `-` for stdin |
| `-o, --output PREFIX` | required | output path prefix |
| `-z, --gzip` | off | gzip-compress output |
| `--gzip-level N` | 6 | gzip compression level (1..9) |
| `--interleaved-in` | off | read alternating R1/R2 from `-f` |
| `-m, --min-overlap N` | 20 | minimum overlap length |
| `-M, --max-overlap N` | 0 | maximum overlap length (0 = unbounded) |
| `-l, --min-length N` | 50 | minimum merged length |
| `-L, --max-length N` | 0 | maximum merged length (0 = unbounded) |
| `-k, --kmer-size N` | 12 | seed kmer size (4..31) |
| `-p, --p-value F` | 0.01 | p-value cutoff for accepting a merge |
| `--p-match F` | 0.25 | null match probability under random alignment |
| `--no-bonferroni` | off | skip Bonferroni correction over candidates |
| `-q, --quality-offset N` | 33 | input Phred offset (33 or 64); output is always +33 |
| `--match REWARD` | 1 | match reward in scoring |
| `--mismatch PENALTY` | 1 | mismatch penalty in scoring |
| `-s, --strict` | off | reject merges with high-Q disagreements |
| `--strict-q N` | 20 | Phred above which a disagreement counts as "high-Q" |
| `-d, --dovetail` | off | allow dovetail / adapter readthrough alignments |
| `-a, --adapter-forward SEQ` | none | adapter expected at R1 3' end (trimmed before merge) |
| `-A, --adapter-reverse SEQ` | none | adapter expected at R2 3' end; defaults to `-a` |
| `--adapter-min-match N` | 8 | min bp of adapter overlap to trigger trim |
| `--adapter-match-frac F` | 0.9 | min match fraction over the adapter overlap |
| `-Q, --min-mean-qual N` | 0 (off) | reject if mean merged Phred below N |
| `--max-ee F` | 0 (off) | reject if expected errors `Σ 10^(-q/10)` exceeds F |
| `--max-mono-frac F` | 0 (off) | reject if any single base exceeds fraction F |
| `--qtrim N` | 0 (off) | trim merged 3' bases below Phred N |
| `-t, --threads N` | 1 | worker threads |
| `-b, --batch-size N` | 1024 | pairs per work batch |
| `-v, --verbose` | off | per-pair status to stderr |
| `-h, --help` | | this help |
| `--version` | | print version |

## Output

For each pair, pbj writes to one of three files:

- successful merge &rarr; `*.assembled.fastq[.gz]` with the consensus sequence
- failed merge (no overlap, low quality, filter rejection) &rarr;
  `*.unassembled.forward.fastq[.gz]` (R1) and
  `*.unassembled.reverse.fastq[.gz]` (R2)

Output preserves input order. Multi-threaded runs produce byte-identical
output to single-threaded runs.

Merged-read headers use R1's "core id" (everything before the first
space, with `/1` or `/2` stripped). Casava 1.8+ direction tokens
(`1:N:0:ATCG`) are dropped from merged records to avoid confusing
downstream re-pairing tools. Unassembled records keep their full
original headers.

## How it works

Four stages per pair:

1. **Kmer seed.** R1's `k`-mers are indexed in an open-addressing hash.
   `k`-mers from `reverse_complement(R2)` probe the index. Each match
   gives a candidate overlap offset.

2. **Quality-weighted scoring.** Each candidate offset is scored by
   summing `+reward * min(q1,q2)` on matches and `-penalty * min(q1,q2)`
   on mismatches.

3. **Statistical test.** The best-scoring candidate must pass a
   seed-adjusted binomial test. Under the null hypothesis (random
   alignment given the kmer seed), the `overlap_len - k` non-seed
   positions are i.i.d. Bernoulli(`p_match`). p-value is multiplied by
   the number of candidates (Bonferroni). Reject if the corrected
   p-value exceeds `--p-value`.

4. **Consensus.** Per overlap position, agreement gives summed quality
   (capped at Phred 93). Disagreement gives the higher-Q base with a
   Bayesian-derived consensus quality under an independent-error model
   with uniform base prior.

### Adapter trimming

When `-a` or `-A` is provided, each read is scanned 3'-anchored for a
partial-prefix match to the adapter. Match positions must extend to the
read end (`pos + matched_length == read_len`), which corresponds to the
biology of adapter readthrough. If at least `--adapter-min-match` bases
match with at least `--adapter-match-frac` fraction, the read is
truncated at that position. `N` in either the read or the adapter is
treated as a wildcard.

The leftmost qualifying position is taken, which trims the longest
acceptable match.

### Post-merge filters

Three filters can be enabled independently, applied to the merged read
in this order:

1. `--qtrim N` strips 3' bases with Phred below N. If trimming brings
   the merged length below `--min-length`, the merge is rejected.
2. `--max-ee F` rejects if `Σ 10^(-q/10)` over the trimmed merged read
   exceeds F (USEARCH/vsearch convention).
3. `--max-mono-frac F` rejects if any of A, C, G, T, or N exceeds
   fraction F (catches homopolymers and N-rich reads).

## Defaults rationale

The defaults `min_overlap = 20` and `kmer_size = 12` are deliberately
stricter than PEAR's (PEAR uses 10 / no kmer seed). The seed-adjusted
p-value test requires the overlap to extend `k` bases beyond the seed
for the binomial to be sensitive, so `min_overlap` must exceed `k` with
margin for sequencing errors.

If your insert distribution sits at or below 20bp overlap, increase
read length, accept some loss, or lower both `-k` and `-m` together
(smaller k inflates the spurious-merge rate slightly).

`p_match = 0.25` assumes uniform A/C/G/T composition. For GC-biased
pools (DNA-storage oligos, primer-heavy amplicons), set `--p-match` to
`Σ f_i²` over the four base frequencies in your data.

## Performance

Single-threaded throughput on a recent x86 desktop runs ~130K pairs/sec
on 150bp reads and ~95K pairs/sec on 300bp reads. Multi-threaded
fork-join scales near-linearly to 8 threads (~470K pairs/sec on 150bp).
The benchmark script is at `bench/throughput.sh`.

Memory is bounded by per-thread workspace plus a per-batch buffer.
Default `--batch-size 1024 -t 4` runs in well under 100 MiB.

## Limits

- Per-record length cap: 256 MiB. Realistic paired-end records are far
  below this.
- Maximum Phred output: 93 (FASTQ convention).
- Input pairs must be in matching order: R1[i] is the mate of R2[i].
  pbj checks Casava 1.8+ direction tokens on the first pair and aborts
  with exit 2 if R2 was passed in place of R1.
- Composition defaults to uniform unless `--p-match` is set.
- pbj is paired-end specific. For long-read merging (PacBio, ONT), use a
  different tool.

## Exit codes

```
0  success (even with zero merged pairs)
1  argument error (invalid flag, missing required input)
2  I/O error (cannot open file, write failure, mismatched pairs)
3  memory error
```

## License

MIT. Vendored `src/kseq.h` is also MIT (Heng Li / GRL).
