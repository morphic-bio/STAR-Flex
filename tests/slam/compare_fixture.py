#!/usr/bin/env python3
"""
Compare STAR-Slam output against GRAND-SLAM reference fixture.

Expected STAR-Slam columns (tab-separated, header required):
  Gene (or GeneID), ReadCount, Conversions, Coverage, NTR (or MAP)
"""

import argparse
import csv
import gzip
import math
import sys


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def find_col(header, names):
    lower = {name.lower(): i for i, name in enumerate(header)}
    for name in names:
        idx = lower.get(name.lower())
        if idx is not None:
            return idx
    return None


def find_suffix(header, suffixes):
    for suffix in suffixes:
        for i, name in enumerate(header):
            if name.endswith(" " + suffix) or name.endswith(suffix):
                return i
    return None


def parse_reference(path):
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError("Reference file is empty")

        gene_idx = find_col(header, ["Gene", "GeneID"])
        read_idx = find_suffix(header, ["Readcount"])
        conv_idx = find_suffix(header, ["Conversions"])
        cov_idx = find_suffix(header, ["Coverage"])
        map_idx = find_suffix(header, ["MAP"])

        missing = [name for name, idx in [
            ("Gene", gene_idx),
            ("Readcount", read_idx),
            ("Conversions", conv_idx),
            ("Coverage", cov_idx),
            ("MAP", map_idx),
        ] if idx is None]
        if missing:
            raise ValueError("Reference missing columns: " + ", ".join(missing))

        data = {}
        for row in reader:
            if not row or len(row) <= map_idx:
                continue
            gene = row[gene_idx]
            if not gene:
                continue
            data[gene] = {
                "readcount": float(row[read_idx]),
                "conversions": float(row[conv_idx]),
                "coverage": float(row[cov_idx]),
                "ntr": float(row[map_idx]),
            }
        return data


def parse_test(path):
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError("Test file is empty")

        gene_idx = find_col(header, ["Gene", "GeneID"])
        read_idx = find_col(header, ["ReadCount", "Readcount"])
        conv_idx = find_col(header, ["Conversions", "ConversionCount", "MismatchCount"])
        cov_idx = find_col(header, ["Coverage", "TCount"])
        ntr_idx = find_col(header, ["NTR", "NTR_MAP", "MAP"])

        missing = [name for name, idx in [
            ("Gene", gene_idx),
            ("ReadCount", read_idx),
            ("Conversions", conv_idx),
            ("Coverage", cov_idx),
            ("NTR", ntr_idx),
        ] if idx is None]
        if missing:
            raise ValueError("Test output missing columns: " + ", ".join(missing))

        data = {}
        for row in reader:
            if not row or len(row) <= ntr_idx:
                continue
            gene = row[gene_idx]
            if not gene:
                continue
            data[gene] = {
                "readcount": float(row[read_idx]),
                "conversions": float(row[conv_idx]),
                "coverage": float(row[cov_idx]),
                "ntr": float(row[ntr_idx]),
            }
        return data


def pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    num = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    den_x = sum((x - mean_x) ** 2 for x in xs)
    den_y = sum((y - mean_y) ** 2 for y in ys)
    den = math.sqrt(den_x * den_y)
    if den == 0.0:
        return float("nan")
    return num / den


def main():
    parser = argparse.ArgumentParser(description="Compare STAR-Slam output to GRAND-SLAM fixture.")
    parser.add_argument("--reference", required=True, help="GRAND-SLAM fixture TSV(.gz)")
    parser.add_argument("--test", required=True, help="STAR-Slam output TSV")
    parser.add_argument("--min-read-count", type=float, default=50.0,
                        help="Minimum read count for NTR correlation")
    parser.add_argument("--count-tol", type=float, default=0.0,
                        help="Tolerance for read/conversion/coverage deltas")
    parser.add_argument("--ntr-abs-max", type=float, default=1e-3,
                        help="Max absolute NTR difference (MAP vs NTR)")
    parser.add_argument("--corr-min", type=float, default=0.999,
                        help="Minimum Pearson correlation for NTR")
    parser.add_argument("--max-report", type=int, default=10,
                        help="Max mismatches to print")
    args = parser.parse_args()

    ref = parse_reference(args.reference)
    test = parse_test(args.test)

    shared = sorted(set(ref.keys()) & set(test.keys()))
    if not shared:
        print("FAIL: no overlapping genes between reference and test")
        return 1

    def collect_deltas(field):
        deltas = []
        for gene in shared:
            delta = abs(ref[gene][field] - test[gene][field])
            if delta > args.count_tol:
                deltas.append((delta, gene, ref[gene][field], test[gene][field]))
        deltas.sort(reverse=True)
        return deltas

    read_deltas = collect_deltas("readcount")
    conv_deltas = collect_deltas("conversions")
    cov_deltas = collect_deltas("coverage")

    ntr_pairs = []
    ntr_bad = []
    for gene in shared:
        if ref[gene]["readcount"] < args.min_read_count:
            continue
        r = ref[gene]["ntr"]
        t = test[gene]["ntr"]
        ntr_pairs.append((r, t, gene))
        if abs(r - t) > args.ntr_abs_max:
            ntr_bad.append((abs(r - t), gene, r, t))
    ntr_bad.sort(reverse=True)

    xs = [r for r, _, _ in ntr_pairs]
    ys = [t for _, t, _ in ntr_pairs]
    corr = pearson(xs, ys) if ntr_pairs else float("nan")

    ok = True
    if read_deltas:
        ok = False
        print(f"FAIL: ReadCount mismatches: {len(read_deltas)}")
        for delta, gene, r, t in read_deltas[:args.max_report]:
            print(f"  {gene}\tref={r}\ttest={t}\tdelta={delta}")
    if conv_deltas:
        ok = False
        print(f"FAIL: Conversions mismatches: {len(conv_deltas)}")
        for delta, gene, r, t in conv_deltas[:args.max_report]:
            print(f"  {gene}\tref={r}\ttest={t}\tdelta={delta}")
    if cov_deltas:
        ok = False
        print(f"FAIL: Coverage mismatches: {len(cov_deltas)}")
        for delta, gene, r, t in cov_deltas[:args.max_report]:
            print(f"  {gene}\tref={r}\ttest={t}\tdelta={delta}")

    if math.isnan(corr) or corr < args.corr_min:
        ok = False
        print(f"FAIL: NTR correlation {corr:.6f} (min {args.corr_min})")
    if ntr_bad:
        ok = False
        print(f"FAIL: NTR abs diff > {args.ntr_abs_max}: {len(ntr_bad)}")
        for delta, gene, r, t in ntr_bad[:args.max_report]:
            print(f"  {gene}\tref={r}\ttest={t}\tdelta={delta}")

    if ok:
        print("PASS: STAR-Slam parity checks")
        print(f"  Genes compared: {len(shared)}")
        print(f"  NTR correlation: {corr:.6f} (min {args.corr_min})")
        if ntr_pairs:
            print(f"  NTR genes (readcount >= {args.min_read_count}): {len(ntr_pairs)}")
        return 0

    return 1


if __name__ == "__main__":
    sys.exit(main())
