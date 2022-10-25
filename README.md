# variant-rc: Variant Read Classification

Tool to classify aligned sequencing reads into bins based upon distinguishing variants (SNPs and Indels) found in parental sequences.

Requirements for classification of sequencing reads:

* VCF file containing variants between two or more references
* aligned reads to be sorted in BAM format

## Installation

Requires: Python >=3.8, pip or pipx

Prefered method of installation into an isolated python enviroment using pipx:
(<https://pypa.github.io/pipx/>)

```
pip install pipx  # install pipx if unavailable

# then use pipx to install directly from this github repo
pipx install git+https://github.com/rnabioco/variant-rc.git
```

but can also be installed to your active python environemtn with vanilla pip:

```
pip install git+https://github.com/rnabioco/variant-rc.git
```

Installation creates a systemwide `variant-rc` command.

## Usage
Command Arguments:
`variant-rc --help`
![img](screenshot1.png)

Read Classification:
`variant-rc <bamfile> <vcf-file>`
![img](screenshot2.png)