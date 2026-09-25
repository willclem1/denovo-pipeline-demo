# Trio De Novo Variant Calling Demo

A small, readable implementation of the core logic behind de novo variant
calling in a parent-child trio: scan a VCF for sites where the child is
heterozygous, both parents are confidently homozygous reference, and basic
QC filters pass.

This is a **teaching demo**, not a production caller. A real pipeline would
add joint genotyping, local realignment around indels, parental mosaicism
checks, population frequency filtering (gnomAD), and consequence annotation
(VEP). The point here is to show the logic clearly in ~150 lines of
dependency-free Python.

## How it works

`denovo_scan.py` applies three gates to every variant:

1. **Genotype gate**: child must be het (`0/1`); both parents must be `0/0`.
   Anything else is a different inheritance pattern, not a de novo.
2. **Confidence gate**: minimum depth and genotype quality for the child
   *and* for each parent. A "hom-ref" parent with 3 reads proves nothing.
3. **Allele balance gate**: the child's alt-allele fraction must fall in a
   plausible heterozygous range (default 0.3-0.7). Strongly skewed ratios
   usually mean sequencing artifact, not biology.

## From demo to real data: the full filter cascade

The three gates above are the core logic, but a production-grade de novo
analysis needs several more filters. This section documents the full cascade
as applied to a real 30x whole-genome parent-child quad (two parents, two
sons; GRCh38), which reduced ~2.8M heterozygous sites per child to 123
high-confidence de novo SNPs. Aggregate counts only; no individual-level
data is included in this repo.

**0. Joint evidence model.** The demo scans one VCF; the real analysis
streamed all four VCFs together (47.65M records). At each position, every
person was resolved as: called genotype, confident hom-ref (a PASS reference
block covering the position), or missing. 373,584 multiallelic or
conflicting-REF positions were excluded outright ("complex").

**1. De novo candidacy (per child, autosomes).** Child het, PASS, biallelic;
GQ ≥ 20; depth 10–120; allele balance 0.25–0.75. Both parents confident
hom-ref (PASS, depth ≥ 10, negligible alt depth). A candidate carried by the
*other* child was discarded: a shared allele is inherited or parental
mosaic, not de novo. chrX was scanned separately with a mother-only rule
(sons inherit X from the mother). Result: ~16K candidates per child.

**2. Novelty.** A true germline de novo is essentially always new to
science, so any candidate with a dbSNP rsID was dropped. This removed ~60%
of what remained: recurrent artifacts love to masquerade as de novos.

**3. Repeat mask.** Candidates in UCSC segmental duplications or
RepeatMasker regions were removed. These regions misalign routinely and
manufacture fake variants. Result: ~700 candidates per child.

**4. Declustering.** Genuine de novos don't cluster, so candidates within
1 kb of each other were dropped.

**5. Sanity check against biology.** Expectation: ~70 de novo SNPs per
genome, roughly half lost to the repeat mask. The surviving counts (37 and
87 SNPs) landed in that range, which is the signal that the list is real
rather than residual noise. Indels went through the same cascade but stayed
artifact-heavy (600+ each vs. ~5–10 expected biologically), so they were set
aside as untrustworthy without realignment-based validation.

**6. Population frequency (gnomAD v4.1).** A genuine de novo should be absent
or vanishingly rare in ~1.2M population genomes. One candidate was too
common and was cut as probable parental allelic dropout. Final: **123
high-confidence de novo SNPs.**

**7. Consequence annotation (VEP).** Each survivor was annotated for
predicted molecular consequence (missense, UTR, intronic, …) to prioritize
follow-up: protein-altering and regulatory candidates first.

## Try it

The repo ships with a synthetic trio VCF (fictional chromosome `chrSYN`,
entirely made-up data) containing 300 background variants, 5 spiked-in de
novos, and 4 decoy Mendelian violations designed to fail QC:

```bash
# regenerate the synthetic data (optional; a copy is in example/)
python generate_trio_vcf.py --out example/synthetic_trio.vcf --seed 42

# run the scan
python denovo_scan.py --vcf example/synthetic_trio.vcf \
    --trio child,father,mother \
    --out example/denovo_candidates.tsv
```

Expected output:

```
Scanned 309 variants
Candidate de novos: 5
Filtered: 85 parental genotype, 215 child genotype, 4 QC
```

All 5 spiked-in de novos are recovered, and all 4 decoys (low GQ, low
depth, skewed allele balance, weak parental evidence) are filtered. See
`example/denovo_candidates.tsv` for the result table.

## Tuning

```
--min-dp 10          minimum child depth
--min-gq 30          minimum child genotype quality
--min-parent-dp 10   minimum depth in each parent
--min-parent-gq 30   minimum genotype quality in each parent
--ab-low 0.3         minimum child alt allele balance
--ab-high 0.7        maximum child alt allele balance
```

## Files

- `denovo_scan.py`: the caller (pure Python, no dependencies)
- `generate_trio_vcf.py`: synthetic trio VCF generator
- `example/synthetic_trio.vcf`: example input (fake data)
- `example/denovo_candidates.tsv`: example output
