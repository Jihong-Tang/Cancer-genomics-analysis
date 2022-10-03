[TOC levels=1-3]: ##

# Table of Contents
- [Table of Contents](#table-of-contents)
- [Name](#name)
- [Purpose](#purpose)
- [Part0 Introduction and Preparation](#part0-introduction-and-preparation)
  - [0.1 Introduction](#01-introduction)
  - [0.2 Getting Setup](#02-getting-setup)
    - [0.2.1 Installing Conda Environment](#021-installing-conda-environment)
    - [0.2.2 Installing Homebrew (Alternative choice)](#022-installing-homebrew-alternative-choice)
    - [0.2.3 Setting Up the Folder Structure](#023-setting-up-the-folder-structure)
- [Part1 Sequencing Data Quality Control](#part1-sequencing-data-quality-control)
- [Part2 Genome Alignment](#part2-genome-alignment)
- [Part3 Somatic Mutation Calling](#part3-somatic-mutation-calling)
- [Part4 Copy Number Variation Analysis](#part4-copy-number-variation-analysis)
- [Author](#author)
- [Reference](#reference)


# Name
Cancer Genomics Analysis Pipeline(CGAP): Notes for the computational pipeline focusing on the Cancer Genomics DNA-seq data

# Purpose
- Make record for the in-house DNA-seq data analysis pipeline including several modules for future usage:
  - Sequencing data Quality Control
  - Genome Alignment
  - Somatic Mutation Calling 
  - Copy Number Variation Analysis
- Take notes about cancer genomics data analysis for future review.

# Part0 Introduction and Preparation
## 0.1 Introduction
**UPDATING**
## 0.2 Getting Setup
### 0.2.1 Installing Conda Environment
Anaconda/Miniconda is one powerful and user-friendly package manager for Python. Miniconda is the smaller version of Anaconda with basic functions.

Using conda environment management tools, we can not only install and manage Python packages with different verisons, but also create virual environments and install various bioinformatics tools through access to large bioinformatics repositories(eg. [Bioconda](https://bioconda.github.io/))

The methods to download and install conda can be found [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). For example, we can use the following commands to install miniconda on MacOS machines.
```bash
# download the Miniconda3 installer to the home directory, only for MacOS
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O ~/miniconda.sh

# Run the miniconda installation in slient mode
bash ~/miniconda.sh -b -p $HOME/miniconda

# Choose YES for the question 'Do you wish the installer to initialize Miniconda3 by running conda init?', 
# then the bash profile will be changed automatically.

# Add bioinformatic channels for downloading required packages
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda # important channel to download and manage useful bioinformatics tools
```

### 0.2.2 Installing Homebrew (Alternative choice)
Homebrew is another powerful package manager that can install *the stuff you need* that Apple (or your Linux system) didn't have. The most well-known example is *wget*.

The methods to download and install brew can be found [here](https://brew.sh/). For example, we can use the following commands to install brew on MacOS machines.

```bash
# download the Homebrew installer and install brew
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

In addition, we could install the brew on Linux machines as following, referred from [wang-q Ubuntu](https://github.com/wang-q/ubuntu#install-linuxbrew).
```bash
echo "==> Install linuxbrew, copy the next *ONE* line to terminal"
bash -c "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/install/master/install.sh)"

test -d ~/.linuxbrew && PATH="$HOME/.linuxbrew/bin:$HOME/.linuxbrew/sbin:$PATH"
test -d /home/linuxbrew/.linuxbrew && PATH="/home/linuxbrew/.linuxbrew/bin:/home/linuxbrew/.linuxbrew/sbin:$PATH"

if grep -q -i linuxbrew $HOME/.bashrc; then
    echo "==> .bashrc already contains linuxbrew"
else
    echo "==> Update .bashrc"

    echo >> $HOME/.bashrc
    echo '# Linuxbrew' >> $HOME/.bashrc
    echo "export PATH='$(brew --prefix)/bin:$(brew --prefix)/sbin'":'"$PATH"' >> $HOME/.bashrc
    echo "export MANPATH='$(brew --prefix)/share/man'":'"$MANPATH"' >> $HOME/.bashrc
    echo "export INFOPATH='$(brew --prefix)/share/info'":'"$INFOPATH"' >> $HOME/.bashrc
    echo "export HOMEBREW_NO_ANALYTICS=1" >> $HOME/.bashrc
    echo "export HOMEBREW_NO_AUTO_UPDATE=1" >> $HOME/.bashrc
    echo >> $HOME/.bashrc
fi

source $HOME/.bashrc
```

Most of the used softwares and tools can be downloaded through both package managers, and it is highly recommended to use either of them to manage the related softwares and bioinformatics tools.

### 0.2.3 Setting Up the Folder Structure
Organizing is essential to proper reproducible research. During the processing and analysis steps, we will get many output files. A simple and clear folder structure is preferred not only for your own future usage, but also for other researchers' repeats.

The folder structure used in this pipeline will be summarized here. 

```bash
── ${project_WES}/
  │   └── reference/                <- Reference genome file (.FASTA)
  │       ├── hg19bwa/              <- hg19 version reference genome files
  │       ├── TCGA_GRCh38_bwa_ref/  <- hg38 version reference genome files
  │ 
  │   └── all_fq/                   <- Working Directory of Sequencing Data QC 
  │       ├── log/                  <- Log from Sequencing Data QC module
  │       ├── qc_res/               <- QC result from Sequencing Data QC module
  │ 
  │   └── all_mapping_hg19/         <- Working Directory of Genome Alignment: hg19 version 
  │       ├── log/                  <- Log from Genome Alignment module
  │       ├── md_metrix/            <- MarkDuplicate metrix from alignment duplicates marking step
  │       ├── md_tmp/               <- Location for temporary files from alignment duplicates marking step
  │       ├── tmp/                  <- Location for temporary files from Genome Alignment module
  │ 
  │   └── all_mapping_hg38/         <- Working Directory of Genome Alignment: hg38 version 
  │       ├── log/                  <- Log from Genome Alignment module
  │       ├── md_metrix/            <- MarkDuplicate metrix from alignment duplicates marking step
  │       ├── md_tmp/               <- Location for temporary files from alignment duplicates marking step
  │       ├── tmp/                  <- Location for temporary files from Genome Alignment module
  │ 
  │   └── savi_hg38/                <- Working Directory of Somatic Mutation Calling: hg38 version 
  │       ├── log/                  <- Log from savi_hg38 module
  │       ├── .../                  <- Other detailed subfolders containing working files, will be shown in future section
  │
  │   └── cnvkit_hg19/              <- Working Directory of CNV Analysis
  │       ├── log/                  <- Log from CNV Analysis module
  │       ├── .../                  <- Other detailed subfolders containing working files, will be shown in future section
  │ 
```

**Important Notes**
- The reference folder contains the reference genome used, two different versions of reference genome file will be covered in this pipeline (hg19 and hg38). For WangLab HPC3 user, the folder can be linked as following: 
  ```bash
  ln -s /scratch/PI/jgwang/jtangbd/reference ${project_WES}/reference/
  ```
  For other users, the reference genomes can be downloaded from the following pages: [hg19](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/), [hg38](https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files)
- Two different versions of reference genome are included in this pipeline. They are both used to generate different versions of genome alignment result files(BAMs). Afterwards, the bams are used differently for following mutation calling and CNV analysis module based on the reference genome version requirements for different tools.
- For the Somatic Mutation calling and CNV analysis modules, detailed folder arrangements will be explained in [Part3](#part3-somatic-mutation-calling) and [Part4](#part4-copy-number-variation-analysis).
# Part1 Sequencing Data Quality Control

# Part2 Genome Alignment

# Part3 Somatic Mutation Calling 

# Part4 Copy Number Variation Analysis


# Author
Jihong Tang &lt;jtangbd@connect.ust.hk&gt;instructed by @[Prof. Jiguang Wang](https://github.com/JiguangWang), @[Dr. Quanhua Mu](https://github.com/qhmu) and @[Dong Song](https://github.com/ForceField17)
# Reference 

