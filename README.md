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
  - [1.1 Setting Up Tools](#11-setting-up-tools)
  - [1.2 FastQ Files Preparation](#12-fastq-files-preparation)
  - [1.3 Automatic QC with fastp](#13-automatic-qc-with-fastp)
    - [1.3.1 Working Scripts](#131-working-scripts)
    - [1.3.2 Scripts Usage Examples](#132-scripts-usage-examples)
- [Part2 Genome Alignment](#part2-genome-alignment)
  - [2.1 Setting Up Tools](#21-setting-up-tools)
  - [2.2 Genome Alignment with bwa](#22-genome-alignment-with-bwa)
    - [2.2.1 Working Scripts](#221-working-scripts)
    - [2.2.2 Scripts Usage Example](#222-scripts-usage-example)
  - [2.3 Mark Duplicates with Picard](#23-mark-duplicates-with-picard)
    - [2.3.1 Working Scripts](#231-working-scripts)
    - [2.3.2 Scripts Usage Example](#232-scripts-usage-example)
- [Part3 Somatic Mutation Calling](#part3-somatic-mutation-calling)
  - [3.1 Setting Up Tools](#31-setting-up-tools)
  - [3.2 Somatic Mutations Calling with SAVI](#32-somatic-mutations-calling-with-savi)
- [Part4 Copy Number Variation Analysis](#part4-copy-number-variation-analysis)
  - [4.1 Preparation of input files](#41-preparation-of-input-files)
    - [4.1.1 Conda environment](#411-conda-environment)
    - [4.1.2 Access file for reference genome](#412-access-file-for-reference-genome)
    - [4.1.3 Baits file](#413-baits-file)
    - [4.1.4 Target and Antitarget bed files](#414-target-and-antitarget-bed-files)
  - [4.2 CNVkit Commands](#42-cnvkit-commands)
    - [4.2.1 Coverage](#421-coverage)
    - [4.2.2 Reference;Fix;Segment;Export](#422-referencefixsegmentexport)
  - [4.3 CNV Profile Generation](#43-cnv-profile-generation)
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

**`Index File 1st`**

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

**`Index File 2nd`**

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
## 1.1 Setting Up Tools
The main bioinformatics tool used in this module is fastp. The fastp is a tool designed to provide fast all-in-one preprocessing for FastQ files.
It can automatically filter reads, detect reads adapters, make adapter trimming and output QC results based on the demand. It covers the functions of well known QC tools like FastQC and Trim Galero, and makes the whole preprocessing work automatically in one step, which is very convenient.
More details about fastp could be found in their [user manual page](https://github.com/OpenGene/fastp). 

fastp could be easily downloaded and installed with bioconda:
```bash
# note: the fastp version in bioconda may be not the latest
conda install -c bioconda fastp
```

or download the prebuilt binary versions for usage:
```bash
# download the latest build
wget http://opengene.org/fastp/fastp
chmod a+x ./fastp

# or download specified version, i.e. fastp v0.23.1
wget http://opengene.org/fastp/fastp.0.23.1
mv fastp.0.23.1 fastp
chmod a+x ./fastp
```
## 1.2 FastQ Files Preparation
FastQ file is often used to store the high throughput sequencing reads, and contains useful information including quality scores for each base. 
Deatiled description about FastQ files could be found [here](https://compgenomr.github.io/book/fasta-and-fastq-formats.html). In order to start the genomics analysis,
we need to prepare the fastq files for the project. There are mainly two methods to get the fastq files for your project. The first one is producing your own data by sequencing(eg. sequencing service from companies). 
The second one is downloading data from public database, and it is very useful when you would like to take advantage of the previous work from the community.

However, some datasets in the public database will not contain fastq files, instead, bam files are often contained. In order to control the analysis parameters and pipeline, we could extract the fastq files from the downloaded bam files, and then deal with them with your own pipeline. 
In order to do so, we could use several different tools, and *SamToFastq* in Picard is taken as example here. 

```bash
#!/usr/bin/bash

picard=/scratch/PI/jgwang/dsongad/software/Picard/picard.jar


NAME=$1
NEW=$2

# The example script is setted for the HPC3 User in WangLab HKUST
/scratch/PI/jgwang/dsongad/software/jdk-15.0.2/bin/java -Djava.io.tmpdir=tmp_$NAME \
        -Xmx36g -XX:+UseParallelGC -XX:ParallelGCThreads=2 \
        -jar ${picard} SamToFastq \
        INPUT=../${NAME}.bam \
        F=${NEW}.R1.fastq \
        F2=${NEW}.R2.fastq
```

## 1.3 Automatic QC with fastp
In order to facilitate the parallel processing of multiple samples and make the pipeline reproducible, we could summarize the commands in scripts. All the scripts we created in this whole pipeline are template structured, which will read parameters from the command line if we use them individually. 
In this way, we only need to change several simple parameters to control the whole pipeline to deal with different samples. 

The arrangement logic of the scripts for part 1-3 is the same and can be summarized as follows. For each step, one bash script and one perl script will be used. The bash script contains the detailed command to use different bioinformatics tools, eg. fastp in this step. The required parameters of the bash script often determines which sample will be processed each time. 

The perl script contains the template structure to read in the index files we created in the above [section](#024-creating-cohort-index-files), the detailed command to run the bash script and how to control the parallel processing procedure. The perl script often requires the path of the index file, the working directory, the start index and end index to the working samples as parameters. 

### 1.3.1 Working Scripts
**[`Script1: do_fastp.sh(Bash script)`](/scripts/part1-Sequencing-Data-QC/do_fastp.sh)**

The script will read in two parameters, and use them as the raw name prefix of fastq files and the changed name prefix. The paired fastq files are expected to end with suffix *.R1.fq.gz* and *.R2.fq.gz* to use this script, and the script can also be modified based on your own file naming style. 

The QC results will be stored in the subfoler *qc_res*, and the raw fastq files will be removed after the automatic preprocessing. 
```bash
#!/usr/bin/bash

rawName=$1
cleanName=$2_fp

#mkdir -p qc_res

fastp -w 2 -i ./$rawName\.R1.fq.gz -o ./$cleanName\.R1.fq.gz \
 -I ./$rawName\.R2.fq.gz -O ./$cleanName\.R2.fq.gz \
 -q 15 -u 40 --length_required 45 \
 --detect_adapter_for_pe \
 -h ./qc_res/fastp_$cleanName\.html -j ./qc_res/fastp_$cleanName\.json -R ./qc_res/report_$cleanName\.txt

rm ./$rawName\.R1.fq.gz ./$rawName\.R2.fq.gz
```

**[`Script2: run_fastp.pl(Perl script)`](/scripts/part1-Sequencing-Data-QC/runfastp.pl)**

The script will read in three parameters, and use them as index file path, start index and end index respectively. The previous created index file will be used here to control the parallel processing procedure. The logic is that we could use start index and end index to choose the specific samples to process each time.

For example, when we want to run the first two samples to test the scripts, we could set the start index as 0, the end index as 1. By this way, the perl script will loop from row index 0 to 1 in the index file to read in the index information and pass them to the bash script. 
Then the bash script will process the two samples one by one based on the parameters received. 
```perl
#!/usr/bin/perl

#please run this job in a new clear fold
die "usage: perl runfastp.pl <Data_Index> <start> <end> \n" if @ARGV!= 3;

$index_table = $ARGV[0];



open INDEX, $index_table;
$i = 0;
$id = 0;
while($line = <INDEX>){
	chomp($line);
	my @temp=split('\t',$line);
	$sample[$i][0] = $temp[1];
	$sample[$i][1] = $temp[2];
	$i++;
}
close INDEX;


$start = $ARGV[1];
$end = $ARGV[2];

for($j=$start;$j<=$end;$j++){         
        `sh ./do_fastp.sh $sample[$j][0] $sample[$j][0]`; # Important command, call the bash script based on parameters
	$completed_No = $j;
        print "Job completed!!! No.$completed_No : $sample[$j][1]\t$sample[$j][0]\n\n";
}
```

### 1.3.2 Scripts Usage Examples
For simple usage or script test, we could directly run the perl script with index paramters. It is highly recommended to record all the log information by redirecting them to another log file. 

The example command used to call perl script for simple test usage could be shown as follows. The command will call the scripts to do fastp work on the first sample only and store the log information into the subfoler *log*.

```bash
perl run_fastp.pl list_SP_fq_bwa.txt 0 0 2>./log/log.fp_0_0
```

**[`Script3: slurm.fp_all(Slurm script)`](/scripts/part1-Sequencing-Data-QC/slurm.fp_all)**

For the slurm cluster users to run multiple samples, the example sbatch script could be arranged as follows. 
```bash
#!/bin/bash
##SBATCH -J bwa24  #Slurm job name
#SBATCH --mail-user=njutangjihong@gmail.com 
#SBATCH --mail-type=end
##SBATCH -p x-gpu-share
##SBATCH -o tra_MD.out
##SBATCH -e tra_MD.err
##SBATCH -N 1 -n 24
#SBATCH --cpus-per-task=20
##SBATCH --exclusive

perl runfastp.pl list_SP_fq_bwa.txt $1 $2 2>./log/log.fp_$1_$2
```

Taken the demo project as example, there are totally 94(0-93) samples in the index files, and we could process them parallelly on four computing nodes as follows:
```bash
sbatch -p x-gpu-share -J fp0 slurm.fp_all 0 23 
sbatch -p x-gpu-share -J fp24 slurm.fp_all 24 47 
sbatch -p x-gpu-share -J fp48 slurm.fp_all 48 71
sbatch -p x-gpu-share -J fp72 slurm.fp_all 72 93 
```
# Part2 Genome Alignment

## 2.1 Setting Up Tools
Genome alignment and following duplicated marking are essential steps in the genomics analysis pipeline to get the bam files. The Bam files are the binary version of Sam files, and the Sam files are tab-delimited text file that contains sequence alignment data. 
Detailed information about SAM/BAM files could be found [here](http://samtools.github.io/hts-specs/).

The bioinformatics tools used in this module are bwa and Picard. Bwa is one of the most famous alignment tools, and is considered as the most useful one. Detailed information about the bwa tools can be found [here](https://github.com/lh3/bwa).

Bwa can also be installed easily through conda:
```bash
conda install -c bioconda bwa
```
In addition, bwa can be installed as follows:
```bash
git clone https://github.com/lh3/bwa.git
cd bwa; make
./bwa index ref.fa
./bwa mem ref.fa read-se.fq.gz | gzip -3 > aln-se.sam.gz
./bwa mem ref.fa read1.fq read2.fq | gzip -3 > aln-pe.sam.gz
```

Picard is a set of Java command line tools for manipulating high-throughput sequencing (HTS) data and formats. Well noted manual and notes have been provided by the developers, and they can be found [here](https://broadinstitute.github.io/picard/).
It is highly recommended to follow their instructions to download and install the Picard and set up java in order to use specific commands. 

For WangLab HKUST HPC3 users, the Picard tool can be used by linking to the following path:
```bash
/scratch/PI/jgwang/dsongad/software/Picard/picard.jar
```
## 2.2 Genome Alignment with bwa
The scripts arrangement logic is the same as previously stated in the section [1.3 Automatic QC with fastp](#13-automatic-qc-with-fastp). 
### 2.2.1 Working Scripts
**[`Script1: do_bwa_dna.sh(Bash script)`](/scripts/part2-Genome-Alignment/do_bwa_dna.sh)**
The script is expected to receive four parameters from user input, representing sample ID, R1 fastq file name, R2 fastq file name and user defined result bam prefix name.
The alignment result bam will be sorted and indexed after the bwa commands, and only the sorted bam will be remained, accompanied by its own index .bai file. 
```bash
#!/usr/bin/bash
A='WangLab'
Sample_ID=$1
FAQ1=$2
FAQ2=$3
result=$4


echo [Directory] `pwd`
echo [Machine] `uname -n`
echo [Start] `date`
echo [args] $*
time1=$( date "+%s" )


#line_q1=`zcat $FAQ1 | wc -l`
#line_q2=`zcat $FAQ2 | wc -l`

bwa mem -t 20 -T 0 -R "@RG\tCN:$A\tID:$Sample_ID\tSM:$Sample_ID" /scratch/PI/jgwang/jtangbd/reference/TCGA_GRCh38_bwa_ref/GRCh38.d1.vd1.fa $FAQ1 $FAQ2 2> ./tmp/tmp.$result | samtools view -Sbh1 -@ 8 - -o $result\.bam

#### mapping duplicates
###counting the num of bam reads

samtools sort -@ 8 -O bam -o $result.sorted.bam -T the_temp_$result $result\.bam
samtools index $result\.sorted.bam 
rm $result\.bam
#rm $FAQ1
#rm $FAQ2
#samtools idxstats result_test.sorted.bam |head 
#samtools view -H /home/qmu/projects/AsianpGBM/TCGAGBM/C484.TCGA-06-6699-10A-01D-1845-08.2_gdc_realn.bam |grep "^@PG"
#samtools flagstat result_test.sorted.bam
#samtools view -f 4 result_test.sorted.bam chr7 | wc -l
#samtools mpileup ../raw_data/results_B_DNA/IGCT_S002_B.sorted.bam -r chr10:63207680 -f /home/dsongad/software/my_lib/STAR_index_hg38_RefSeq_hg38/hg38.fa |head -n20

time2=$( date "+%s" )
echo [deltat] $(( $time2 - $time1 ))
echo [End]
```

**[`Script2: run_bwa_mapping.pl(Perl script)`](/scripts/part2-Genome-Alignment/run_bwa_mapping.pl)**
The script will read in four parameters, and use them as index file path, fastq file path, start index and end index respectively. 
The previous created index file will be used here to control the parallel processing procedure. The logic is that we could use start index and end index to choose the specific samples to process each time. 

The parameter `<fastq-dir>` is essential to be provided to locate the fastq files, and user can choose to delete the raw fastq files after alignment for storage room saving. 
```perl
#!/usr/bin/perl

#please run this job in a new clear fold
die "usage: perl run_bwa_mapping.pl <Data_Index> <fastq-dir> <START> <END> \n" if @ARGV!= 4;

$index_table = $ARGV[0];
$root = $ARGV[1];


open INDEX, $index_table;
$i = 0;
$id = 0;
while($line = <INDEX>){
	chomp($line);
	my @temp=split('\t',$line);
	$sample[$i][0] = $temp[1];
	$sample[$i][1] = $temp[2];
	$i++;
}
close INDEX;

#`mkdir Dir_tmp_1`;
#chdir "./Dir_tmp_1";
#
$start = $ARGV[2];
$end = $ARGV[3];

for($j=$start;$j<=$end;$j++){
        `sh ./do_bwa_dna.sh $sample[$j][0] $root/$sample[$j][0]_fp.R1.fq.gz $root/$sample[$j][0]_fp.R2.fq.gz $sample[$j][1]`;
	$completed_No = $j;
        print "bwa-MEM completed!!! No.$completed_No : $sample[$j][1]\t$sample[$j][0]\n";
}
#chdir "../";
```

### 2.2.2 Scripts Usage Example
For simple usage or script test, we could directly run the perl script with index paramters. It is highly recommended to record all the log information by redirecting them to another log file. 

The example command used to call perl script for simple test usage could be shown as follows. The command will call the scripts to do bwa mapping on the first sample only and store the log information into the subfoler *log*.

```bash
perl run_bwa_mapping.pl list_SP_fq_bwa.txt 0 0 2>./log/log.bwa_0_0
```

**[`Script3: slurm.bwa_all(Slurm script)`](/scripts/part1-Sequencing-Data-QC/slurm.bwa_all)**
For the slurm cluster users to run multiple samples, the example sbatch script could be arranged as follows.
```bash
#!/bin/bash
##SBATCH -J bwa24  #Slurm job name
#SBATCH --mail-user=njutangjihong@gmail.com 
#SBATCH --mail-type=end
##SBATCH -p x-gpu-share
##SBATCH -o tra_MD.out
##SBATCH -e tra_MD.err
##SBATCH -N 1 -n 24
#SBATCH --cpus-per-task=20
#SBATCH --exclusive

perl run_bwa_mapping.pl list_SP_fq_bwa.txt ../all_fq $1 $2 2>./log/log.bwa_$1_$2
```
Taken the demo project as example, there are totally 94(0-93) samples in the index files, and we could process them parallelly on four computing nodes as follows:
```bash
sbatch -p x-gpu-share -J bwa0 slurm.bwa_all 0 23 
sbatch -p x-gpu-share -J bwa24 slurm.bwa_all 24 47 
sbatch -p x-gpu-share -J bwa48 slurm.bwa_all 48 71
sbatch -p x-gpu-share -J bwa72 slurm.bwa_all 72 93 
```
## 2.3 Mark Duplicates with Picard
The scripts arrangement logic is the same as previously stated in the section [1.3 Automatic QC with fastp](#13-automatic-qc-with-fastp). 
### 2.3.1 Working Scripts
**[`Script1: do_mark_duplicates.sh(Bash script)`](/scripts/part2-Genome-Alignment/do_mark_duplicates.sh)**
The script is expected to receive only one parameter from user input, the prefix name of the bam file. After the mark duplicates work, the bam files will be ended with .MD. surfix, and the .MD.bam files will be sorted and indexed again. 
The raw bam files will be removed and only the .MD.bam and .MD.bam.bai files will be remained. 

One important note: The mark duplicates step is implemented in the same directory after the bwa mapping step, and must be used in the same folder if no change is made to the script.
```bash
#!/usr/bin/bash

#set -euxo pipefail
#export LD_LIBRARY_PATH=/home/qmu/tools/lib:/opt/intel/advisor_xe_2016.1.40.463413/lib64:/home/qmu/tools/xz-5.2.3/tmp/lib 
#picard=/home/share/jgwang/softwares/picard/build/libs/picard.jar
picard=/scratch/PI/jgwang/dsongad/software/Picard/picard.jar

NAME=$1

mkdir md_tmp/tmp_$NAME

java -Djava.io.tmpdir=md_tmp/tmp_$NAME \
	-Xmx120g -XX:+UseParallelGC -XX:ParallelGCThreads=4 \
	 -jar ${picard} MarkDuplicates \
	INPUT=${NAME}.sorted.bam \
	OUTPUT=${NAME}.sorted.MD.bam \
	METRICS_FILE=./md_metrix/${NAME}.metrics.txt \
	ASSUME_SORTED=true \
	REMOVE_DUPLICATES=true \
	VALIDATION_STRINGENCY=LENIENT

rm ${NAME}.sorted.bam
rm ${NAME}.sorted.bam.bai
samtools index ${NAME}.sorted.MD.bam
#rm ${NAME}.sorted.MD.bam.txt
```
**[`Script2: run_mark_duplicates.pl(Perl script)`](/scripts/part2-Genome-Alignment/run_mark_duplicates.pl)**
The script will read in four parameters, and use them as index file path, fastq file path, start index and end index respectively. 
The previous created index file will be used here to control the parallel processing procedure. The logic is that we could use start index and end index to choose the specific samples to process each time. 

The parameter setting is the same to the script `run_bwa_mapping.pl` in order to simplify the usage of the index files. 

```perl
#!/usr/bin/perl

#please run this job in a new clear fold
die "usage: perl run_mark_duplicates.pl <Data_Index> <fastq-dir> <START> <END> \n" if @ARGV!= 4;

$index_table = $ARGV[0];
$root = $ARGV[1];


open INDEX, $index_table;
$i = 0;
$id = 0;
while($line = <INDEX>){
	chomp($line);
	my @temp=split('\t',$line);
	$sample[$i][0] = $temp[1];
	$sample[$i][1] = $temp[2];
	$i++;
}
close INDEX;

#`mkdir Dir_tmp_1`;
#chdir "./Dir_tmp_1";
#
$start = $ARGV[2];
$end = $ARGV[3];

for($j=$start;$j<=$end;$j++){          #run the first 100 samples as a try
	`sh ./do_mark_duplicates.sh $sample[$j][1]`;
	$completed_No = $j;
        print "MarkingDup completed!!! No.$completed_No : $sample[$j][0]\t$sample[$j][1]\n\n";
}
#chdir "../";
```
### 2.3.2 Scripts Usage Example
For simple usage or script test, we could directly run the perl script with index paramters. It is highly recommended to record all the log information by redirecting them to another log file. 

The example command used to call perl script for simple test usage could be shown as follows. The command will call the scripts to do picard duplicates marking on the first sample only and store the log information into the subfoler *log*.

```bash
perl run_mark_duplicates.pl list_SP_fq_bwa.txt 0 0 2>./log/log.md_0_0
```

**[`Script3: slurm.md_all(Slurm script)`](/scripts/part1-Sequencing-Data-QC/slurm.md_all)**
For the slurm cluster users to run multiple samples, the example sbatch script could be arranged as follows.


```bash
#!/bin/bash
##SBATCH -J md24  #Slurm job name
#SBATCH --mail-user=njutangjihong@gmail.com 
#SBATCH --mail-type=end
##SBATCH -p x-gpu-share
##SBATCH -o tra_MD.out
##SBATCH -e tra_MD.err
##SBATCH -N 1 -n 24
#SBATCH --cpus-per-task=20
#SBATCH --exclusive

perl run_mark_duplicates.pl list_SP_fq_bwa.txt ./ $1 $2 2>./log/log.md_$1_$2
```
Taken the demo project as example, there are totally 94(0-93) samples in the index files, and we could process them parallelly on four computing nodes as follows:
```bash
sbatch -p x-gpu-share -J md0 slurm.md_all 0 23 
sbatch -p x-gpu-share -J md24 slurm.md_all 24 47 
sbatch -p x-gpu-share -J md48 slurm.md_all 48 71
sbatch -p x-gpu-share -J md72 slurm.md_all 72 93 
```

# Part3 Somatic Mutation Calling 

## 3.1 Setting Up Tools

## 3.2 Somatic Mutations Calling with SAVI

```bash
#!/bin/bash

export LIBRARY_PATH=/scratch/PI/jgwang/dsongad/software/ncurses/lib/:$LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scratch/PI/jgwang/dsongad/software/ncurses/lib/
export PATH=/scratch/PI/jgwang/dsongad/software/samtools-1.2:$PATH
export PATH=$PATH:/scratch/PI/jgwang/dsongad/software/snpEff
export LIBRARY_PATH=/scratch/PI/jgwang/dsongad/software/ncurses/lib/:$LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scratch/PI/jgwang/dsongad/software/ncurses/lib/
#conda activate savi_hg38_py2.7

OUTPUT=$1
BD=$2
T1=$3


V=/scratch/PI/jgwang/dsongad/software/my_lib/tmp_vcf_lib
S=/scratch/PI/jgwang/dsongad/IGCT/run_savi/SAVI
R=/scratch/PI/jgwang/dsongad/software/my_lib/star_genome_d1_vd1_gtfv22/GRCh38.d1.vd1.fa

for i in {1..22} {X,Y,M}; do
	mkdir -p ${OUTPUT}/${i}
        ${S}/savi_SONG_allregion.py --bams ${BD},${T1} \
                --names Nor,T1 \
                --ann GRCh38.p7.RefSeq --memory 40 \
                --ref ${R} -c 2:1 \
                --region chr${i} -v --superverbose \
                --noncoding --mapqual 20 --mindepth 4 \
                --outputdir ${OUTPUT}/${i} \
                --annvcf ${V}/Cosmic_2020_hg38/Song_hg38_Cosmic.sorted.vcf.gz,${V}/Cosmic_2020_hg38/Song_hg38_ClinVar.sorted.vcf.gz,${V}/Cosmic_2020_hg38/Song_hg38_PancancerGermline.sorted.vcf.gz,${V}/Cosmic_2020_hg38/Song_hg38_cbio.vcf.gz,${V}/normal_hg38/Song_hg38_All_SNP.sorted.vcf.gz,${V}/normal_hg38/Song_hg38_TOPMED.sorted.vcf.gz,${V}/normal_hg38/IGCTfamB.sorted.vcf.gz,${V}/normal_hg38/IGCTpatB.sorted.vcf.gz,${V}/normal_hg38/Song_hg38_GATKnormal.vcf.gz,${V}/normal_hg38/Song_hg38_HGDP_1KG.sorted.vcf.gz,${V}/normal_hg38/Song_hg38_MuTectPON.vcf.gz,${V}/normal_hg38/Song_hg38_dbSNP_ALFA.sorted.vcf.gz,${V}/normal_hg38/Song_hg38_gnomAD.sorted.vcf.gz,${V}/normal_hg38/Song_hg38_219normals.vcf.gz,${V}/normal_hg38/Song_China_map.sorted.vcf.gz,${V}/normal_hg38/Song_hg38_ExAC.sorted.vcf.gz,${V}/normal_hg38/Song_hg38_NGDC_SNP.vcf.gz,${V}/normal_hg38/Song_hg38_meganormal.vcf.gz 1> ${OUTPUT}/${i}/err.out 2>${OUTPUT}/${i}/err.log &
done
wait

echo SAVI done!

echo Start merging...

O1=${OUTPUT}/merge.report.coding.PDfilter.txt
O2=${OUTPUT}/merge.report.coding.txt


head -n 1 ${OUTPUT}/1/report/report.coding.PDfilter.txt > ${O1}
head -n 1 ${OUTPUT}/1/report/report.coding.txt > ${O2}

for i in {1..22} {X,Y,M}; do
        F1=${OUTPUT}/${i}/report/report.coding.PDfilter.txt
	F2=${OUTPUT}/${i}/report/report.coding.txt
	
        if [ -e ${F1} ]; then
                tail -n +2 ${F1} >> ${O1}
        else
                echo ${F1} does not exist
        fi

        if [ -e ${F2} ]; then
                tail -n +2 ${F2} >> ${O2}
        else
                echo ${F2} does not exist
        fi

done
```

```perl
#!/usr/bin/perl

#please run this job in a new clear fold
die "usage: perl run_savi.pl <Data_Index> <START> <END> \n" if @ARGV!= 3;

$index_table = $ARGV[0];


open INDEX, $index_table;
$i = 0;
$id = 0;
while($line = <INDEX>){
	chomp($line);
	my @temp=split('\t',$line);
	$sample[$i][0] = $temp[1];
	$i++;
}
close INDEX;

#`mkdir Dir_tmp_1`;
#chdir "./Dir_tmp_1";
#
$start = $ARGV[1];
$end = $ARGV[2];

for($j=$start;$j<=$end;$j++){
	`sh do_savi_WES_2sample_GRCh38.p7.RefSeq.sh  $sample[$j][0] ../all_mapping_hg38/$sample[$j][0]_B_WES_hg38.sorted.MD.bam ../all_mapping_hg38/$sample[$j][0]_T_WES_hg38.sorted.MD.bam`;
	$completed_No = $j;
        print "bwa-MEM completed!!! No.$completed_No : $sample[$j][0]\n";
}
#chdir "../";
```

```bash
#!/bin/bash
##SBATCH -J savi24  #Slurm job name
#SBATCH --mail-user=njutangjihong@gmail.com 
#SBATCH --mail-type=end
##SBATCH -p x-gpu-share
##SBATCH -o tra_MD.out
##SBATCH -e tra_MD.err
##SBATCH -N 1 -n 24
#SBATCH --cpus-per-task=20
#SBATCH --exclusive

perl run_savi.pl list_SP_savi_CNV.txt $1 $2 2> log/log.savi_$1_$2
```

```bash
sbatch -p x-gpu-share -J savi0 slurm.bwa_all 0 11
sbatch -p x-gpu-share -J savi24 slurm.savi_all 12 23 
sbatch -p x-gpu-share -J savi48 slurm.savi_all 24 35
sbatch -p x-gpu-share -J savi72 slurm.savi_all 36 46 
```
# Part4 Copy Number Variation Analysis

[5._109.07.08-2_湯硯安博士-51.pdf](CNVkit%20pipeline%208760af4ce2ce4cfaa5410789947fa423/5._109.07.08-2_%25E6%25B9%25AF%25E7%25A1%25AF%25E5%25AE%2589%25E5%258D%259A%25E5%25A3%25AB-51.pdf)

[https://cnvkit.readthedocs.io/en/stable/pipeline.html](https://cnvkit.readthedocs.io/en/stable/pipeline.html)

## 4.1 Preparation of input files

### 4.1.1 Conda environment

### 4.1.2 Access file for reference genome

Baits files and Access files,

For access file, hg38 and hg19 should use different calculated bed files. The calculate command is shown as follows;

```bash
cnvkit.py access /scratch/PI/jgwang/jtangbd/reference/hg19bwa/hg19.fa \
-o access.hg19.bed
```

### 4.1.3 Baits file

### 4.1.4 Target and Antitarget bed files

```bash
cnvkit.py autobin ../all_mapping_hg19/*B*.bam \
-t S07604514_Regions_hg19_220825.bed \
-g access.hg19.bed \
--target-output-bed SP_20220910_binauto.target.bed \
--antitarget-output-bed SP_20220910_binauto.antitarget.bed \
2>log/log_autobin_auto
```

## 4.2 CNVkit Commands

![CNVkit Pipeline](/cnvkit_workflow.webp)

### 4.2.1 Coverage

For each bam file, we need to calculate the coverage for target and antitarget regions, therefore, there will be two .cnn file as results for each bam file. 

```bash
#!/bin/bash

PathBam=/scratch/PI/jgwang/jtangbd/projects/SP_WES/all_mapping_hg19
PathCNV=/scratch/PI/jgwang/jtangbd/projects/SP_WES/cnvkit_hg19

for i in `cut -f2 list_SP_savi_CNV.txt`;do
        #i=`echo $n | cut -d "\t" -f 2`
        cnvkit.py coverage ${PathBam}/$i\_B_WES_hg19.sorted.MD.bam SP_20220910_binauto.target.bed -p 20 -o ${PathCNV}/coverage/$i\_B.target.cnn
done

for i in `cut -f2 list_SP_savi_CNV.txt`;do
        #i=`echo $n | cut -d "\t" -f 2`
        cnvkit.py coverage ${PathBam}/$i\_B_WES_hg19.sorted.MD.bam SP_20220910_binauto.antitarget.bed -p 20 -o ${PathCNV}/coverage/$i\_B.antitarget.cnn
done

for i in `cut -f2 list_SP_savi_CNV.txt`;do
        #i=`echo $n | cut -d "\t" -f 2`
        cnvkit.py coverage ${PathBam}/$i\_T_WES_hg19.sorted.MD.bam SP_20220910_binauto.target.bed -p 20 -o ${PathCNV}/coverage/$i\_T.target.cnn
done

for i in `cut -f2 list_SP_savi_CNV.txt`;do
        #i=`echo $n | cut -d "\t" -f 2`
        cnvkit.py coverage ${PathBam}/$i\_T_WES_hg19.sorted.MD.bam SP_20220910_binauto.antitarget.bed -p 20 -o ${PathCNV}/coverage/$i\_T.antitarget.cnn
done
```

### 4.2.2 Reference;Fix;Segment;Export

The total pipeline for step 2 to 5 could be summarized in the following bash scripts

```bash
#!/bin/bash

PathCoverage=/scratch/PI/jgwang/jtangbd/projects/SP_WES/cnvkit_hg19/coverage
FileReference=/scratch/PI/jgwang/jtangbd/reference/hg19bwa/hg19.fa
# step0 - make target and antitarget files
# step1 - cnvkit coverage for all samples

for i in `cut -f2 list_SP_savi_CNV.txt`;do
        # step2 - cnvkit reference for paired sample
        cnvkit.py reference ${PathCoverage}/$i\_B*.cnn -f ${FileReference} -o ./reference/$i\_reference.cnn

        # step3 - cnvkit fix
        cnvkit.py fix  ${PathCoverage}/$i\_T.target.cnn  ${PathCoverage}/$i\_T.antitarget.cnn ./reference/$i\_reference.cnn -i $i -o ./subtraction/$i\_fix.cnr

        # step4 - cnvkit segment
        cnvkit.py segment subtraction/$i\_fix.cnr -p 20 --drop-low-coverage -m cbs -t 0.000001 -o ./segmentation/$i\_segment.cbs.cns

        # step5 - cnvkit export and transfer results to seg files
        cnvkit.py call ./segmentation/$i\_segment.cbs.cns -o ./results/$i\_call.cbs.cns
        cnvkit.py export seg ./results/$i\_call.cbs.cns -o ./results/$i\_cbs.seg

        #tail -n+2 ./results/$i\_cbs.seg >> ./results/CNV_CGPA_103samples_pair_ref.seg
done
```

## 4.3 CNV Profile Generation

```bash
# Generate hg19 gene lists from gencode gtf files 
cat gencode.v41lift37.annotation.gtf | sed '/^#/d' | awk '$3=="gene"' | \
awk -F "\t" '{split($9, a, ";"); OFS="\t"; print $1, $4, $5, $7, a[2], a[3]}' | \
awk -F "\t" '{split($5, b, " "); split($6, c, " "); OFS="\t"; print $1,$2,$3,$4,b[2],c[2]}' | \
sed 's/"//g' > GRCh37_gencode_v41_genes.txt
```

```bash
# Generate Protein coding genes without sex chromosome genes
cat GRCh37_gencode_v41_genes.txt | sed '/chr[X|Y]/d' | awk '$5=="protein_coding"' \
> GRCh37_gencode_v41_genes_Protein_noSex.txt
```

```bash
# sort CNV result file for a single sample 
cat GRCh37_gencode_v41_genes_Protein_noSex.txt | sort -k1,1 -k2,2n > GRCh37_gencode_v41_genes_Protein_noSex_sort.txt

cat ../results/P370376_cbs.seg | awk '{OFS="\t"; print $2,$3,$4,$5,$6,$1}' | \
sort -k1,1 -k2,2n > P370376_sort.bed
```

```bash
# Map for single sample
bedtools map -a GRCh37_gencode_v41_genes_Protein_noSex_sort.txt \
-b P370376_sort.bed -c 5 | awk '{print $6"\t"$7}' | sort -k1 | \
sed "1s/^/gene\t$var1\n/" > P370376_CNVprofile.txt
```

The total bash script:

```bash
#!/bin/bash

PathSeg=/scratch/PI/jgwang/jtangbd/projects/SP_WES/cnvkit_hg19/results

for i in `cut -f2 list_SP_savi_CNV.txt`;do
        #i=`echo $n | cut -d "\t" -f 2`
        cat ${PathSeg}/${i}\_cbs.seg | awk '{OFS="\t"; print $2,$3,$4,$5,$6,$1}' | \
        sort -k1,1 -k2,2n > ./sort/${i}\_sort.bed
        # Map for single sample
        var1=${i}
        bedtools map -a GRCh37_gencode_v41_genes_Protein_noSex_sort.txt -b ./sort/${i}\_sort.bed -c 5 | \
        awk '{print $6"\t"$7}' | sort -k1 | \sed "1s/^/gene\t$var1\n/" > ./profile/${i}\_CNVprofile.txt
done

mkdir -p tmp
cut -f1 ./profile/P370376_CNVprofile.txt > gn.txt 
for i in `cut -f2 list_SP_savi_CNV.txt`;do
        cut -f2 ./profile/${i}\_CNVprofile.txt > ./tmp/tmp.${i}
done
paste -d "\t" ./tmp/tmp* > tmp.txt 
paste -d "\t" gn.txt tmp.txt > SP_46samples_CNVPatients.txt
rm gn.txt tmp.txt
```


# Author
Jihong Tang &lt;jtangbd@connect.ust.hk&gt;instructed by @[Prof. Jiguang Wang](https://github.com/JiguangWang), @[Dr. Quanhua Mu](https://github.com/qhmu) and @[Dong Song](https://github.com/ForceField17)
# Reference 
- Part0 Introduction and Preparation
  - [Conda User Guide](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html)
  - [Brew User Guide](https://docs.brew.sh/)
  - [Tom Battaglia's RNAseq-workflow](https://github.com/twbattaglia/RNAseq-workflow)
- Part1 Sequencing Data Quality Control
  - Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu; fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560
- Part2 Genome Alignment
  - Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168]
  - “Picard Toolkit.” 2019. Broad Institute, GitHub Repository. https://broadinstitute.github.io/picard/; Broad Institute
- Part3 Somatic Mutation Calling
- Part4 Copy Number Variation Analysis

