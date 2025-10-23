
# _longTBprofiler_

LongTBprofiler is workflow/framework utilizing long-read Oxford Nanopore data to classify Mycobacterium tuberculosis isolates into lineages and predict drug resistance.

# Software requirements

The analysis is performed in an Unix enviroment (Mac/Linux). The main tools utilized in this pipeline are:
- Nanoplot (https://github.com/wdecoster/NanoPlot)
- Nanocomp (https://github.com/wdecoster/nanocomp)
- Minimap2 (https://github.com/lh3/minimap2)
- freebayes (https://github.com/freebayes/freebayes) 
- TB-Profiler (https://github.com/jodyphelan/TBProfiler)
- TbLG (https://github.com/dbespiatykh/tblg)

Detailed descriptions of commands offered by these tools are available at their corresponding github pages, their websites, as well as available after typing the following in the command line: "tool name [function] -h". For example:
  ```
  tb-profiler profile -h 
  ```

# Installation guide
First we need to install mamba. Mamba is a package manager, that will help us easily install all of the bioinformatic tools required for the longTBprofiler.

## Mamba

I. Go to: https://conda-forge.org/download/

II. Choose the correct version for your system (x86_64/arm64)*, then right click and "copy link adress"

*u can check your system architecture by:
  ```
  arch
  ```

III. Go to home directory in terminal (with command "cd")

  ```
  cd
  ```

IV. Download the miniforge 

  For x86_64 architecture
  ```
  wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
  ```

IV. After a moment u have your miniforge downloaded. Check the name of file by viewing the content of current directory with "ls"

  ```
  ls
  ```


V. Copy the name of the downloaded file. It should be "Miniforge3-Linux-x86_64.sh " or something similar.

VI. Run the downloaded script to install mamba.

  ```
  bash Miniforge3-Linux-x86_64.sh
  ```

VII. Accept the terms of service when prompted. 

VIII. Do not allow for initiation of mamba for each shell session.

XIX. Check if mamba has been successfully installed.

  ```
  mamba
  ```

It should print the version of mamba, show the location of installation and list options. If you see these, your mamba has been installed correctly and u can proceed. If not, please inform one of the instructors. 
  ```
  Version: 2.3.0



/home/walter/miniconda3/bin/mamba [OPTIONS] [SUBCOMMAND]


OPTIONS:
  -h,     --help              Print this help message and exit
          --version

Configuration options:
          --rc-file FILE1 FILE2...
                              Paths to the configuration files to use
          --no-rc             Disable the use of configuration files
          --no-env            Disable the use of environment variables

  ```


## Packages 

1. Download the txt file ("secret-ingredients") from the provided google drive link and move it to your home directory in ubuntu.
  
2. Create an enviroment containing required packages.

  ```
  mamba create -n longmtb --file secret-ingredients.txt
  ```
  ```
  mamba activate longmtb
  ```
  ```
  mamba install pip
  ```
  ```
  conda install bioconda::nanocomp
  ```
  ```
  conda install -c bioconda nanoplot
  ```
  ```
  pip install tblg
  ```

3. Accept package installations/updates if prompted. 

# Analysis

Now, we have environment ready with all the neccessary packages. Finally we can start analysis. 

Transfer your files to correct directories.

1. To check home directory, type:
```
pwd
```
2. Create subdirectories inside your home directory for your fastq files
```
mkdir -p ./data/WGS
mkdir -p ./data/target
```

3. Move your WGS and targeted sequencing (amplicon) fastq files to correct directories.

## Quality check - generating report for long read data

*change the sample name if its different than barcode83, 51..

1. Activate mamba environment with all the packages
  ```
  mamba activate longmtbc
  ```

2.  Generate quality report using nanocomp and nanoplot

  #### for single sample 
  ```
  NanoPlot --fastq ./data/WGS/barcode83.fastq.gz -o ./QC/WGS/nanoplot_quality_report
  NanoPlot --fastq ./data/target/barcode51.fastq.gz -o ./QC/target/nanoplot_quality_report
  ```
  #### for multiple samples
  ```
  NanoComp --fastq ./data/WGS/barcode83.fastq.gz ./data/target/barcode51.fastq.gz -o ./QC/WGS/nanocomp_quality_report
  ```

## Read Alignment/mapping

   I. Download Reference genome from NCBI

  ```
  mkdir -p ./data/reference
  ```
  ```
  cd ./data/reference
  ```
  ```
  wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz
  ```
  ```
  mv GCF_000195955.2_ASM19595v2_genomic.fna.gz reference.fna.gz
  ```
  ```
  gunzip -d reference.fna.gz
  ```
  ```
  head -10 reference.fna
  ```
  ```
  cd ..
  ```
  
  II. Map samples to the reference

  ```
  mkdir ./mapped_varcall
  mkdir -p ./mapped_varcall/WGS
  mkdir -p ./mapped_varcall/target
  ```
  

  ### Minimap2

  #### for WGS data
  
  ```
  minimap2 -a -x lr:hq ./data/reference/reference.fna ./data/WGS/barcode83.fastq.gz | samtools sort -O bam -o ./mapped_varcall/WGS/barcode83_sorted.bam 
  ```
  #### for targeted sequencing data
  
  ```
  minimap2 -a -x lr:hq ./data/reference/reference.fna ./data/target/barcode51.fastq.gz | samtools sort -O bam -o ./mapped_varcall/target/barcode51_sorted.bam
  ```
  To check if mapping was sucessfull (u should see sorted.bam files)
  ```
  ls 
  ./mapped_varcall/WGS
  ./mapped_varcall/target
  ```

*if the terminal show you some problem with permissions, type:

  ```
  sudo chmod -R a+rwx /path/to/your/directory
  ```  
 
  III. Indexing the mapped reads

  ### Samtools

  #### for WGS data
  ``` 
  samtools index ./mapped_varcall/WGS/barcode83_sorted.bam
  ```
  #### for targeted sequencing data
  ```
  samtools index ./mapped_varcall/target/barcode51_sorted.bam
  ```

   to check if indexing was succesfull (u should see .bam.bai file)
  ```
  ls 
  ./mapped_varcall/WGS
  ./mapped_varcall/target
  ```

## Variant calling

  ### freebayes

  #### for WGS data
  ```
  freebayes --ploidy 1 -f ./data/reference/reference.fna ./mapped_varcall/WGS/barcode83_sorted.bam > ./mapped_varcall/wgs/barcoded83.vcf 
  ```
  #### for targeted sequencing data
  ```
  freebayes --ploidy 1 -f./data/reference/reference.fna ./mapped_varcall/target/barcode51_sorted.bam > ./mapped_varcall/target/barcoded51.vcf
  ```


## Lineage drug resistance typing

We will use TB-Profiler for drug resistance prediction and lineage identification. Additionally, we will use TbLG for more precise lineage identification. 

  Prepare directories for typing results
  ```
  mkdir ./TB_Profiler
  mkdir -p ./TB_Profiler/WGS  
  mkdir -p ./TB_Profiler/target
  mkdir ./tblg
  mkdir -p ./tblg/WGS
  mkdir -p ./tblg/target
  ```

  ### TB-Profiler

  #### for WGS data
  ```
  tb-profiler profile --platform nanopore --call_whole_genome -1 ./data/WGS/barcode83.fastq.gz --dir ./TB_Profiler/WGS/ --csv
  ```

  #### for targeted sequencing data
  ```
  tb-profiler profile --platform nanopore --call_whole_genome -1 ./data/target/barcode51.fastq.gz --dir ./TB_Profiler/target/ --csv

  ```

  View your drug resistance and lineage typing results in:
  ```
  ls ./TB_Profiler/WGS/
  ls ./TB_Profiler/target/
  ```
  
  

  ### TbLG 

  #### for WGS data
  ```
  cd ./tblg/WGS
  tblg ../../mapped_varcall/WGS/barcoded83.vcf -o csv
  ```
  #### for targeted sequencing data
  ```
  cd ./tblg/target
  tblg ../../mapped_varcall/target/barcoded51.vcf -o csv
  ```



