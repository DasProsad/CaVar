"""
maincavar.py

Main method for cavar
"""

import sys
import re
from cyvcf2 import VCF
from cavar.logger import log_message as logger
from cavar.utils import revcomp, parse_bed, is_indexed
from cavar.utils import subset_vcf, compress_vcf, sort_vcf, index_vcf
from cavar.utils import install_all_vars, crrna_status, pam_status
from cavar.startup import show_startup

def run_cavar(args) -> None:
    """
    Function to run cavar. It will take
    parsed args and do all param validation
    if true, then run the program
    """
    grna = args.grna
    bedfile = args.bedfile
    vcffile = args.vcffile
    pam_re = args.pam_regex
    distance = args.distance
    outbed  = args.outfile

    if not bedfile.exists():
        print("BED file does not exist!")
        sys.exit(3)

    if not vcffile.exists():
        print("VCF file does not exist!")
        sys.exit(3)

    if not grna:
        print("gRNA sequence was not provided!")
        sys.exit(1)

    show_startup(args)

    if not is_indexed(vcffile):
        logger("warning", f"{vcffile} is not indexed, creating index")
        source_index = index_vcf(vcffile)
        if not source_index.exists():
            logger("error", f"Could not index {vcffile}")
            sys.exit(3)
    
    logger("info", f"Subsetting {vcffile} with {bedfile}")
    sub_vcf = subset_vcf(vcffile, bedfile)
    if not sub_vcf.exists():
        logger("error", f"Could not subset {sub_vcf}")
        sys.exit(3)

    logger("info", f"Compressing {sub_vcf}")
    com_sub_vcf = compress_vcf(sub_vcf)
    if not com_sub_vcf.exists():
        logger("error", f"Could not compress {sub_vcf}")
        sys.exit(3)

    logger("info", f"Sorting {com_sub_vcf}")
    com_sr_sub_vcf = sort_vcf(com_sub_vcf)
    if not com_sr_sub_vcf.exists():
        logger("error", f"Could not sort {com_sr_sub_vcf}")
        sys.exit(3)

    logger("info", f"Indexing {com_sub_vcf}")
    sub_vcf_index = index_vcf(com_sub_vcf)
    if not sub_vcf_index.exists():
        logger("error", "Index for subset VCF not found!")
        sys.exit(3)

    vcf = VCF(com_sub_vcf)
    valid_pam = re.compile(pam_re)
    grna_len = len(grna)

    logger("info", f"Writing output to {outbed}")
    with outbed.open(mode = "w", encoding = "utf-8") as oh:
        for chrom, start, end, pts, _, strand in parse_bed(bedfile):
            logger("info", f"Analysing gRNA {chrom}:{start}-{end}")
            # variants are always on + strand
            if strand == "-":
                pts = revcomp(pts)
            mutpts = install_all_vars(vcf, chrom, start, end, pts)
            if strand == "-":
                mutpts = revcomp(mutpts)
                pts = revcomp(pts)
            pam_s = pam_status(valid_pam, pts[-3:], mutpts[-3:])
            crrna_s = crrna_status(grna[-grna_len:], mutpts[-grna_len:], distance)
            score = crrna_s * 10 + pam_s
            if score:
                logger("info", f"Score computed: {score}")
                oh.write(f"{chrom}\t{start}\t{end}\t{pts}:{mutpts}\t{score}\t{strand}\n")
