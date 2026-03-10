<div align="center">
	<img src="https://raw.githubusercontent.com/DasProsad/CaVar/main/assets/logo.png" alt="cavar logo" width="180"/>
	<h1>CaVar</h1>
	<h3>Check status of CRISPR-Cas9 off-targets from genomic variants</h3>
	<a href="https://github.com/dasprosad/cavar"><img src="https://img.shields.io/github/v/release/dasprosad/cavar?logo=github$logoColor=white" alt="GitHub Release"></a>
	<a href="https://pypi.org/project/cavar"><img src="https://img.shields.io/pypi/v/cavar?logo=pypi&logoColor=yellow" alt="PyPI Release"></a>
	<a href="https://github.com/dasprosad/cavar/blob/master/LICENSE"><img src="https://img.shields.io/badge/license-GPLv3-orange.svg" alt="License"></a>
	<a href=""><img src="https://img.shields.io/badge/python-3.9%20|%203.10%20|%203.11%20|%203.12%20|%203.13%20|%203.14-blue" alt="Python versions"></a>
	<a href="https://pypi.org/project/cavar"><img src="https://img.shields.io/pypi/dm/cavar" alt="PyPI Downloads"></a>
	<a href="https://github.com/dasprosad/cavar"><img src="https://img.shields.io/github/repo-size/DasProsad/CaVar" alt="Repo Size"></a>
</div>

## Introduction

CaVar (pronounced as ka-vee-a) is a CLI tool to check status of CRISPR-Cas off-target sites (OTs) from genomic variants. It utilizes output from an upstream OT finding tool, for example, Cas-Offinder and a VCF file to infer status of the OT sites. The following scoring scheme ([Table 1](#table1)) is used to report the status of the PAM (pm). To calculate the protospacer score (ps), a hamming-distance-based metric is used. If the number of mismatches between the supplied gRNA and the variant protospacer sequence of an off-target site is less than the maximum allowed value, the protospacer score (ps) is set to unity; otherwise, it is set to zero. Finally, the off-target score is calculated as: **os := (ps × 10) + pm**.

<div align="center">

| old PAM | new PAM      | Score (pm) |
|---------|--------------|------------|
| valid   | valid        |   0        |
| valid   | valid (same) |   1        |
| invalid | valid        |   2        |
| valid   | invalid      |   3        |
| invalid | invalid      |   4        |

</div>

<a id="table1"></a>
<p align="center"><b>Table 1. Off-target PAM scoring scheme</b></p>

## Dependencies

Cavar requires Python >= 3.9 and depends on the following external libraries

- cyvcf2 == 0.31.4
- pysam == 0.23.3

## Installation

### From PyPI

To install with pip use

```sh
pip install cavar
```

### From source

For installation from the source use the following

```sh
git clone https://github.com/dasprosad/cavar.git
cd cavar
pip install .
```

## Usage and options

### Prerequisites

- **Usage of reference genome**

The reference genome **must** be the same which was used for variant calling. VCF specification states that CHROM field must not contain any white space so, the FASTA headers must be reformated to contain only the FASTA ID and not any other description.

- **Finding off-target sites with Cas-Offinder**

Cas-Offinder usage is listed on its [GitHub project page](https://github.com/snugel/cas-offinder). Use the header-reformatted reference genome for quering with Cas-Offinder. Set the maximum number tolerated mismatches according to your Cas system of choice leaving the PAM as "NNN". An example input file for Cas-Offinder with a maximum of 4 mismatches looks like the following.

```
/home/user/path/to/ref/hg38.fa
NNNNNNNNNNNNNNNNNNNNNN
ATGTTGATGATAGGATGATNNN 4
```

Run Cas-Offinder with

```sh
cas-offinder input.txt C casoffinder_out.tab
```

- **Converting Cas-Offinder result to BED**

There are many ways of converting it to a BED file, but I have used the following `awk` oneliner

First add a BED header (optional)

```sh
echo -e "#BED3\n#CHROM\tSTART\tEND\tNAME\tSCORE\tSTRAND" >casoffinder_out.bed
```

Then use `awk` to format to BED

```awk
awk '$1 !~ /^#/ {print $1, $4, $5}' OFS='\t' casoffinder_out.tab >>casoffinder_out.bed
```

### Cavar usage

- For help use `cavar --help`

### Options

```sh
usage:  cavar [OPTIONS] <grna> <bed file> <vcf file>
        cavar --help

Check status of CRISPR-Cas9 off-targets from genomic variants

positional arguments:
  grna (STR)           grna without the pam sequence
  bedfile (PATH)       path of bedfile
  vcffile (PATH)       path of vcffile

options:
  -h, --help           show this help message and exit
  -v, --version        show program's version number and exit
  -p, --pam-regex STR  pam as regex (default: [ATCG]GG)
  -d, --distance INT   maximum number of mismatches in crRNA (default: 4)
  -o, --outfile PATH   name of outfile (default: outfile.bed)
```

### Examples

1. Find off-target status of a gRNA with a given Cas9

- Find OTs of the gRNA with Cas9 across the reference genome

```sh
cas-offinder input.txt C casoffinder_ots.tab
```

- Convert the `casoffinder_ots.tab` to BED

- Check status of those OTs with genomic variants

```sh
cavar --outfile p1_ots.bed ATGTTGATGATAGGATGAT casoffinder_ots.bed p1.vcf.gz
```

## Release

- **1.0b1** First beta release.

## Issues

If you have found any bugs or would like to request any new features please report it on [CAVAR](https://github/dasprosad/cavar/issues).

## License

CaVar is distributed under the GNU Public License v3 (GPL3). You should have received a copy of the license with CaVar.
