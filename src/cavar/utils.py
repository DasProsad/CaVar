"""
utils.py

Utility script for cavar
"""

import sys
import re
import pysam
import pysam.bcftools
from cyvcf2 import VCF
from pathlib import Path

#VALIDAPMS = re.compile(r"(?:[ACTG][AG]|[ATCG]GA)").fullmatch

def is_indexed(vcffile: Path) -> bool:
    """
    Checks if the VCF is indexed
    """
    index = vcffile.with_suffix(vcffile.suffix + ".tbi")
    status = index.exists()
    return status

def parse_bed(bedfile: Path) -> tuple:
    """
    Yields a tuple of entries by parsing BED
    """
    with bedfile.open(mode = "r", encoding = "utf-8") as bh:
        for line in bh:
            if line.strip() == "" or line.startswith("#"):
                continue
            chrom, start, end, name, score, strand, *_ = line.strip().split("\t")
            yield chrom, int(start), int(end), name, int(score), strand

def compress_vcf(vcffile: Path) -> Path:
    """
    Compress a vcf
    """
    outname = vcffile.with_suffix(vcffile.suffix + ".gz")
    pysam.tabix_compress(vcffile, outname, force = True)
    return outname

def sort_vcf(vcffile: Path) -> Path:
    """
    Sort a compressed vcf file with bcftools
    """
    out = vcffile.with_name(vcffile.stem.split(".")[0] + "_sorted.vcf.gz")
    pysam.bcftools.sort("-O", "z", "-o", str(out), str(vcffile), catch_stdout = False)
    return out

def index_vcf(vcffile: Path) -> Path:
    """
    Index the vcf
    """
    vcfgzindex = vcffile.with_suffix(vcffile.suffix + ".tbi")
    pysam.tabix_index(str(vcffile), preset = "vcf", force = True)
    return vcfgzindex

def subset_vcf(vcffile: Path, bedfile: Path) -> Path:
    """
    Subsets the vcf ot bed positions
    """
    vcfout = vcffile.with_name(vcffile.stem.split(".")[0] + "_subset.vcf")
    vcf = pysam.VariantFile(vcffile)
    out = pysam.VariantFile(vcfout, mode = "w", header = vcf.header)
    for chrom, start, end, *_ in parse_bed(bedfile):
            start = start + 1
            if chrom not in vcf.header.contigs:
                continue
            for rec in vcf.fetch(chrom, start, end):
                if len(rec.ref) == 1 and all(len(a) == 1 for a in rec.alts or []):
                    out.write(rec)
    out.close()
    return vcfout

def revcomp(seq: str) -> str:
    """
    Reverse complement
    """
    revcomp = seq.translate(str.maketrans("ACTGNactgn", "TGACNtgacn"))[::-1]
    return revcomp

def pam_status(pam_re: re.Pattern, pam: str, pam_alt: str) -> str:
    """
    Identity: 0
    Preserved: 1
    Created: 2
    Destroyed: 3
    Invalid: 4
    """
    _is_valid_pam = lambda seq: bool(pam_re.fullmatch(seq))
    valid_pam = _is_valid_pam(pam.upper())
    valid_pam_alt = _is_valid_pam(pam_alt.upper())
    if valid_pam and valid_pam_alt:
        return 0 if valid_pam  == pam_alt else 1
    if not valid_pam and valid_pam_alt:
        return 2
    if valid_pam and not valid_pam_alt:
        return 3
    return 4

def hamm_dist(s1: str, s2: str) -> int:
    """
    Calculates hamming distance of two strs
    """
    str1 = s1.upper()
    str2 = s2.upper()
    if len(str1) != len(str2):
        raise ValueError("Strings must be of equal length!")
    return sum(c1 != c2 for c1, c2 in zip(str1, str2))

def return_var(vcfh: VCF, chrom: str, start: int, end: int) -> tuple:
    """
    Reads from a vcf handle and return tuple of vars
    inside the grna
    """
    wstart = start + 1
    wend = end + 1
    for var in vcfh(f"{chrom}:{wstart}-{wend}"):
        if var.FILTER is None or var.FILTER == "PASS":
            yield var, var.start

def install_all_vars(vcfh: VCF, chrom: str, start: int, end: int, ref: str) -> str:
    """
    Finds all SNVs within (start, end), install
    them in ref and return the final mut ref seq.

    It assumes:
        VCF is position sorted
        Only SNVs
        First ALT is considered
        Ref matches VCF
    """
    snvs: list[tuple[str, int]] = []
    for var, v_start in return_var(vcfh, chrom, start, end):
        if len(var.REF) == 1 and len(var.ALT[0]) == 1:
            snvs.append((var.ALT[0], v_start))
    snvs.sort(key = lambda x: x[1])
    mutref = ref.upper()
    for alt, v_pos in snvs:
        offset = v_pos - start
        if offset < 0 or offset >= len(mutref):
            raise ValueError("Variant offset out of bounds")
        mutref = mutref[:offset] + alt.lower() + mutref[offset + 1: ]
    return mutref

def crrna_status(crrna: str, crrna_alt: str, d: int) -> int:
    """
    If a there are less than or eq d mismatches, the
    crrna is preserved
    """
    distance = d
    if hamm_dist(crrna, crrna_alt) <= distance:
        return 1
    return 0
