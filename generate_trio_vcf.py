"""Generate a small synthetic trio VCF for the de novo demo.

Creates entirely fake data on a fictional chromosome (chrSYN): random
background SNVs with Mendelian-consistent genotypes, plus a handful of
spiked-in candidate de novo variants and a few decoy Mendelian violations
that good QC filters should reject.

Usage:
    python generate_trio_vcf.py --out example/synthetic_trio.vcf --seed 42
"""
import argparse
import random

CHROM = "chrSYN"
CHROM_LEN = 1_000_000
N_BACKGROUND = 300
BASES = ["A", "C", "G", "T"]


def gt_str(gt):
    a, b = gt
    return f"{a}/{b}"


def make_record(pos, ref, alt, genos):
    """genos: dict sample -> (gt_tuple, ad_ref, ad_alt, dp, gq)."""
    fields = [CHROM, str(pos), ".", ref, alt, "60", "PASS",
              "AF=0.001", "GT:AD:DP:GQ"]
    for s in ("child", "father", "mother"):
        gt, ad_r, ad_a, dp, gq = genos[s]
        fields.append(f"{gt_str(gt)}:{ad_r},{ad_a}:{dp}:{gq}")
    return "\t".join(fields)


def het_geno(dp=None):
    dp = dp or random.randint(25, 60)
    ad_a = int(dp * random.uniform(0.35, 0.65))
    return ((0, 1), dp - ad_a, ad_a, dp, random.randint(90, 99))


def homref_geno(dp=None):
    dp = dp or random.randint(25, 60)
    return ((0, 0), dp, 0, dp, random.randint(90, 99))


def background_variant(pos):
    """A random SNV with Mendelian-consistent trio genotypes."""
    ref = random.choice(BASES)
    alt = random.choice([b for b in BASES if b != ref])
    # parental genotypes: mostly hom-ref, sometimes het
    f_gt = (0, 1) if random.random() < 0.15 else (0, 0)
    m_gt = (0, 1) if random.random() < 0.15 else (0, 0)
    # child inherits one allele from each parent
    c = (random.choice(f_gt), random.choice(m_gt))
    c = (min(c), max(c))
    genos = {
        "child": ((c[0], c[1]),) + het_geno()[1:] if c == (0, 1)
                 else ((0, 0),) + homref_geno()[1:],
        "father": ((f_gt[0], f_gt[1]),) + (het_geno()[1:] if f_gt == (0, 1)
                                           else homref_geno()[1:]),
        "mother": ((m_gt[0], m_gt[1]),) + (het_geno()[1:] if m_gt == (0, 1)
                                           else homref_geno()[1:]),
    }
    return make_record(pos, ref, alt, genos)


def spiked_denovo(pos, ref="A", alt="G", dp=40, gq=99, ab=0.5):
    """A textbook de novo: het child (balanced), confident hom-ref parents."""
    ad_a = int(dp * ab)
    genos = {
        "child": ((0, 1), dp - ad_a, ad_a, dp, gq),
        "father": homref_geno(dp=random.randint(30, 60)),
        "mother": homref_geno(dp=random.randint(30, 60)),
    }
    return make_record(pos, ref, alt, genos)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    rng = random.Random(args.seed)

    header = [
        "##fileformat=VCFv4.2",
        "##source=generate_trio_vcf.py (synthetic demo data)",
        f"##contig=<ID={CHROM},length={CHROM_LEN}>",
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tchild\tfather\tmother",
    ]
    records = []
    positions = sorted(rng.sample(range(1000, CHROM_LEN), N_BACKGROUND))
    for p in positions:
        records.append(background_variant(p))

    # --- spiked-in true de novos (should be called) ---
    for i, p in enumerate([111_111, 222_222, 333_333, 444_444, 555_555]):
        records.append(spiked_denovo(p, ref="C", alt="T",
                                     dp=rng.randint(30, 60),
                                     ab=round(rng.uniform(0.4, 0.6), 2)))
    # --- decoys: Mendelian violations that QC should reject ---
    # low genotype quality in child
    records.append(spiked_denovo(666_666, gq=18))
    # low depth in child
    records.append(spiked_denovo(777_777, dp=6))
    # skewed allele balance (likely artifact)
    records.append(spiked_denovo(888_888, dp=40, ab=0.12))
    # weak parental evidence (father low depth)
    rec = spiked_denovo(999_999)
    parts = rec.split("\t")
    parts[10] = "0/0:8,0:8:99"  # father DP=8
    records.append("\t".join(parts))

    records.sort(key=lambda r: int(r.split("\t")[1]))
    with open(args.out, "w") as f:
        f.write("\n".join(header + records) + "\n")
    print(f"Wrote {len(records)} synthetic variants to {args.out}")


if __name__ == "__main__":
    main()
