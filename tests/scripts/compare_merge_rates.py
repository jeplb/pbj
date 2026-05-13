#!/usr/bin/env python3
"""compare merge rates and sequence agreement between bbmerge, pear, pbj."""

from __future__ import annotations

import argparse
import gzip
import sys
from typing import IO


def open_maybe_gz(path: str) -> IO[str]:
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def read_fastq(path: str) -> dict[str, str]:
    out: dict[str, str] = {}
    with open_maybe_gz(path) as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().rstrip("\n")
            _   = f.readline()
            _   = f.readline()
            name = header[1:].split()[0].rstrip("/12")
            out[name] = seq
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bbmerge", required=True)
    ap.add_argument("--pear", required=True)
    ap.add_argument("--pbj", required=True)
    ap.add_argument("--tolerance", type=float, default=0.02)
    args = ap.parse_args()

    bb:  dict[str, str] = read_fastq(args.bbmerge)
    pr:  dict[str, str] = read_fastq(args.pear)
    pbj: dict[str, str] = read_fastq(args.pbj)

    n_bb  = len(bb)
    n_pr  = len(pr)
    n_pbj = len(pbj)
    print(f"bbmerge merged: {n_bb}")
    print(f"pear    merged: {n_pr}")
    print(f"pbj     merged: {n_pbj}")

    if n_bb == 0:
        print("error: bbmerge produced no merged reads", file=sys.stderr)
        sys.exit(2)

    rate_diff: float = abs(n_pbj - n_bb) / n_bb
    print(f"|pbj - bbmerge| / bbmerge = {rate_diff:.4f}  (tolerance {args.tolerance})")
    fail: bool = rate_diff > args.tolerance

    common: set[str] = set(bb) & set(pr) & set(pbj)
    eq_all = sum(1 for n in common if bb[n] == pr[n] == pbj[n])
    print(f"three-way exact agreement on {len(common)} common reads: {eq_all}")
    if common:
        eq_frac = eq_all / len(common)
        print(f"three-way exact agreement fraction: {eq_frac:.4f}")
        if eq_frac < 0.99:
            print("warning: three-way agreement below 99%", file=sys.stderr)

    if fail:
        print(f"FAIL: pbj merge rate diverges from bbmerge by > {args.tolerance}",
              file=sys.stderr)
        sys.exit(1)
    print("OK")


if __name__ == "__main__":
    main()
