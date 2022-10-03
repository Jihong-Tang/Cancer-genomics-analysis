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
    - [0.2.4 Creating Cohort Index Files](#024-creating-cohort-index-files)
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
The in house cancer genomics analysis pipeline starts from the quanlity control of DNA-seq sequencing data. After aligning the sequencing data to different versions of reference genomes, the pipeline identifies the somatic mutations and copy number variations based on the comparisons between normal/blood and tumor tissues. 

The whole analysis pipeline is implemented across four main procedures based on different bioinformatics tools:
- Sequencing Data Quality Control - [fastp](https://github.com/OpenGene/fastp)
- Genome Alignment - [bwa](https://github.com/lh3/bwa) and [Picard](https://broadinstitute.github.io/picard/)
- Somatic Mutations Calling - [SAVI](https://github.com/WangLabHKUST/SAVI)
- CNV Analysis - [CNVkit](https://cnvkit.readthedocs.io/en/stable/)

The pipeline is initially arranged for the whole exome sequencing(WES) data, but with minor updates on the scripts, it can also be used for the whole genome sequencing(WGS) data. 

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

### 0.2.4 Creating Cohort Index Files
In order to facilitate the cohort data management and parallel processing, we need to create two index files for the WES/WGS sequencing data in this whole pipeline. Both index files are simple .txt files containing the basic naming information for all the sequencing data. 

**Index File 1st**

The first file is used for Sequencing Data QC and Genome Aligenment modules, and it contains three columns, representing row number(0-index), raw name prefix of fastq files, and changed name prefix of bam files, respectively.

Below is the example format, it is expected to contain index rows for both Blood(B) and Tumor(T) samples in this index file. It is recommended to label the blood and tumor samples in specific symbols, such as B and T after the patient number.

The first column is used for easier arrangement of parallel computing, and detailed information will be provided in the next subsection. The second column is the original name prefix of the fastq files, 
and it is highly recommended to rename your fastq files in your own style but noted them well in this index file. 
The third column is the changed name prefix you would like to name your bam files. It is highly recommended to contain DNA-seq methods(WES/WGS) and reference genome version(hg19/hg38) in the prefix for further confirmation.

The demo index file could be found [here](scripts/list_SP_fq_bwa.txt). It is used for the demo SP project data and contains 47 patients, 94 sequencing data(B + T) files.
```bash
0   P001_B    P001_B_WES_hg38
1   P002_B    P002_B_WES_hg38
2   P003_B    P003_B_WES_hg38
3   P004_B    P004_B_WES_hg38
4   P005_B    P005_B_WES_hg38
5   P001_T    P001_T_WES_hg38
6   P002_T    P002_T_WES_hg38
7   P003_T    P003_T_WES_hg38
8   P004_T    P004_T_WES_hg38
9   P005_T    P005_T_WES_hg38
...
```

**Index File 2nd**

The second file is used for Somatic Mutation Calling and CNV Analysis modules, and it contains only two columns, representing row number(0-index), and patient prefix resepectively.
The second index file only needs to contain one index row for each patient, which is different to the first one. 

Below is the example format and the demo index file could be found [here](scripts/list_SP_savi_CNV.txt).

```bash
0   P001
1   P002
2   P003
4   P004
5   P005
...
```
# Part1 Sequencing Data Quality Control

# Part2 Genome Alignment

# Part3 Somatic Mutation Calling 

# Part4 Copy Number Variation Analysis


# Author
Jihong Tang &lt;jtangbd@connect.ust.hk&gt;instructed by @[Prof. Jiguang Wang](https://github.com/JiguangWang), @[Dr. Quanhua Mu](https://github.com/qhmu) and @[Dong Song](https://github.com/ForceField17)
# Reference 
- Part0 Introduction and Preparation
  - [Conda User Guide](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html)
  - [Brew User Guide](https://docs.brew.sh/)
  - [Tom Battaglia's RNAseq-workflow](https://github.com/twbattaglia/RNAseq-workflow)
- Part1 Sequencing Data Quality Control
- Part2 Genome Alignment
- Part3 Somatic Mutation Calling
- Part4 Copy Number Variation Analysis

