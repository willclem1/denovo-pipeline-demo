"""Scan a trio VCF for candidate de novo variants.

A candidate de novo is a site where the child is heterozygous and both
parents are confidently homozygous reference, passing basic QC filters
(depth, genotype quality, allele balance).

This is a teaching implementation: pure Python, no dependencies, easy to
read. A production pipeline would add joint genotyping, local realignment,
parental mosaicism checks, and population frequency filters.

Usage:
    python denovo_scan.py --vcf example/synthetic_trio.vcf \\
        --trio child,father,mother --out example/denovo_candidates.tsv
"""
import argparse
import sys


def parse_gt(gt):
    """Return sorted allele tuple, or None for missing."""
    if gt in (".", "./.", ".|."):
        return None
    gt = gt.replace("|", "/")
    try:
        return tuple(sorted(int(a) for a in gt.split("/")))
    except ValueError:
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf", required=True, help="Input trio VCF")
    ap.add_argument("--trio", required=True,
                    help="Sample names as child,father,mother")
    ap.add_argument("--out", required=True, help="Output TSV path")
    ap.add_argument("--min-dp", type=int, default=10,
                    help="Minimum child depth (default 10)")
    ap.add_argument("--min-gq", type=int, default=30,
                    help="Minimum child genotype quality (default 30)")
    ap.add_argument("--min-parent-dp", type=int, default=10,
                    help="Minimum parental depth (default 10)")
    ap.add_argument("--min-parent-gq", type=int, default=30,
                    help="Minimum parental genotype quality (default 30)")
    ap.add_argument("--ab-low", type=float, default=0.3,
                    help="Minimum child alt allele balance (default 0.3)")
    ap.add_argument("--ab-high", type=float, default=0.7,
                    help="Maximum child alt allele balance (default 0.7)")
    args = ap.parse_args()
    child_name, father_name, mother_name = args.trio.split(",")

    stats = {"scanned": 0, "candidates": 0, "fail_parent_gt": 0,
             "fail_child_gt": 0, "fail_qc": 0}
    out_lines = ["#chrom\tpos\tref\talt\tchild_dp\tchild_gq\tchild_ab\t"
                 "father_dp\tfather_gq\tmother_dp\tmother_gq"]

    with open(args.vcf) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                idx = {s: i for i, s in enumerate(samples)}
                for s in (child_name, father_name, mother_name):
                    if s not in idx:
                        sys.exit(f"Sample '{s}' not in VCF header")
                ci, fi, mi = idx[child_name], idx[father_name], idx[mother_name]
                continue
            stats["scanned"] += 1
            cols = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt = cols[0], cols[1], cols[3], cols[4]
            fmt = cols[8].split(":")
            calls = [dict(zip(fmt, cols[9 + i].split(":"))) for i in (ci, fi, mi)]
            c_gt, f_gt, m_gt = (parse_gt(c["GT"]) for c in calls)

            # Both parents must be confidently hom-ref.
            if f_gt != (0, 0) or m_gt != (0, 0):
                stats["fail_parent_gt"] += 1
                continue
            # Child must be het.
            if c_gt != (0, 1):
                stats["fail_child_gt"] += 1
                continue

            def num(call, key, default=0):
                try:
                    return int(call.get(key, default))
                except (ValueError, TypeError):
                    return default

            c_dp = num(calls[0], "DP")
            c_gq = num(calls[0], "GQ")
            f_dp, f_gq = num(calls[1], "DP"), num(calls[1], "GQ")
            m_dp, m_gq = num(calls[2], "DP"), num(calls[2], "GQ")
            try:
                ad = [int(x) for x in calls[0]["AD"].split(",")]
                ab = ad[1] / sum(ad) if sum(ad) else 0.0
            except (KeyError, ValueError, ZeroDivisionError):
                ab = 0.0

            qc_ok = (
                c_dp >= args.min_dp and c_gq >= args.min_gq
                and f_dp >= args.min_parent_dp and f_gq >= args.min_parent_gq
                and m_dp >= args.min_parent_dp and m_gq >= args.min_parent_gq
                and args.ab_low <= ab <= args.ab_high
            )
            if not qc_ok:
                stats["fail_qc"] += 1
                continue

            stats["candidates"] += 1
            out_lines.append(
                f"{chrom}\t{pos}\t{ref}\t{alt}\t{c_dp}\t{c_gq}\t{ab:.2f}\t"
                f"{f_dp}\t{f_gq}\t{m_dp}\t{m_gq}"
            )

    with open(args.out, "w") as f:
        f.write("\n".join(out_lines) + "\n")

    print(f"Scanned {stats['scanned']} variants")
    print(f"Candidate de novos: {stats['candidates']}")
    print(f"Filtered: {stats['fail_parent_gt']} parental genotype, "
          f"{stats['fail_child_gt']} child genotype, "
          f"{stats['fail_qc']} QC")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
