"""
cli.py

cli for cavar
"""

import argparse
from pathlib import Path
from cavar import __version__
from cavar.maincavar import run_cavar

def build_parser():
    """
    Builds parser for cli
    """
    parser = argparse.ArgumentParser(
        description = "Check status of CRISPR-Cas9 off-targets from genomic variants",
        usage = "\tcavar [OPTIONS] <grna> <bed file> <vcf file>\n\tcavar --help",
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "-v",
        "--version",
        action = "version",
        version = f"%(prog)s {__version__}",
    )

    parser.add_argument(
        "grna",
        type = str,
        metavar = "grna (STR)",
        help = "grna without the pam sequence"
    )

    parser.add_argument(
        "bedfile",
        type = Path,
        metavar = "bedfile (PATH)",
        help = "path of bedfile"
    )

    parser.add_argument(
        "vcffile",
        type = Path,
        metavar = "vcffile (PATH)",
        help = "path of vcffile"
    )

    parser.add_argument(
        "-p",
        "--pam-regex",
        type = str,
        metavar = "STR",
        default = "[ATCG]GG",
        help = "pam as regex"
    )

    parser.add_argument(
        "-d",
        "--distance",
        type = int,
        metavar = "INT",
        default = 4,
        help = "maximum number of mismatches in crRNA"
    )

    parser.add_argument(
        "-o",
        "--outfile",
        type = Path,
        required = False,
        default = "outfile.bed",
        metavar = "PATH",
        help = "name of outfile"
    )

    return parser

def main():
    """
    CLI entrypoint
    """
    parser = build_parser()
    parser.set_defaults(func = run_cavar)
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
