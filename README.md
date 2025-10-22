
# _longTBprofiler_

LongTBprofiler is workflow/framework utilizing long-read Oxford Nanopore data to classify Mycobacterium tuberculosis isolates into lineages and predict drug resistance.
# Software requirments

The analysis is performed in an Unix enviroment (Mac/Linux). The main tools utilized in this pipeline are:
- Nanoplot (https://github.com/wdecoster/NanoPlot)
- Nanocomp
- Minimap2
- freebayes
- TB-Profiler
- TbLR


# Installation guide
First we need to install mamba. Mamba is a package manager that will help us easily install all bioinformatic tools for LongTBprofiler pipeline.

## Mamba

I. Go to: https://conda-forge.org/download/

II. Choose the correct version for your system (x86_64/arm64), then right click and "copy link adress"

III. Go to home directory in terminal (with command "cd")

  ```
  cd
  ```

IV. Download the miniforge

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
  bash "Miniforge3-Linux-x86_64.sh"
  ```

VII. Accept the terms of service when prompted. 

VIII. Allow for initiation of mamba for each shell session.

XIX. Check if mamba has been successfully installed.

  ```
  mamba
  ```

It should print the version of mamba, show the location of installation and list options. If you see these, your mamba has been installed correctly and u can proceed. If not, please inform one of the instructors. 


## Packages 

1. Download the txt file ("secret-ingredients") from the provided google drive link and move it to your home directory in ubuntu.
  
2. Create an enviroment containing required packages.

  ```
  mamba create -n longmtb --file secret-ingredients.txt
  mamba activate longmtb
  mamba install pip nanocomp nanoplot
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

  # for single sample 
  ```
  NanoPlot --fastq ./data/WGS/barcode83.fastq.gz -o ./QC/WGS/nanoplot_quality_report
  NanoPlot --fastq ./data/target/barcode51.fastq.gz -o ./QC/target/nanoplot_quality_report
  ```
  # for multiple samples
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
  gunzip -d reference.fna
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
  ```
  Mapping WGS data
  ```
  ```
  minimap2 -a -x lr:hq ./data/reference/reference.fna ./data/WGS/barcode83.fastq.gz | samtools sort -O bam -o ./mapped_varcall/WGS/barcode83_sorted.bam 
  ```
  ```
  Mapping targeted sequencing data
  ```
  ```
  minimap2 -a -x lr:hq ./data/reference/reference.fna ./data/target/barcode51.fastq.gz | samtools sort -O bam -o ./mapped_varcall/target/barcode51_sorted.bam 
  ```

*if the terminal show you some problem with permissions, type:
  ```
  sudo chmod -R a+rwx /path/to/your/directory
  ```  
 
  III. Indexing the mapped reads

  WGS
  ``` 
  samtools index ./mapped_varcall/WGS/barcode83_sorted.bam
  ```
  targeted sequencing
  ```
  samtools index ./mapped_varcall/target/barcode51_sorted.bam
  ```
  To check if mapping was sucessfull (u should see sorted.bam files)
  
  ```
  ls 
  ./mapped_varcall/WGS
  ./mapped_varcall/target
  ```
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

  #### for WGS
  ```
  freebayes --ploidy 1 -f ./data/reference/reference.fna ./mapped_varcall/WGS/barcode83_sorted.bam > ./mapped_varcall/wgs/barcoded83.vcf 
  ```
  #### for targeted sequencing
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

  for WGS
  ```
  tb-profiler profile --platform nanopore --call_whole_genome -1 ./data/WGS/barcode83.fastq.gz --dir ./TB_Profiler/WGS/
  ```

  for target
  ```
  tb-profiler profile --platform nanopore --call_whole_genome -1 ./data/target/barcode51.fastq.gz --dir ./TB_Profiler/target/

  ```

  View your drug resistance and lineage typing results in:
  ```
  ls ./TB_Profiler/WGS/
  ls ./TB_Profiler/target/
  ```
  
  

  ### TbLG 

  for WGS
  ```
  cd ./tblg/WGS
  tblg ../../mapped_varcall/WGS/barcoded83.vcf -o csv
  ```
  for targeted
  ```
  cd ./tblg/target
  tblg ../../mapped_varcall/target/barcoded51.vcf -o csv
  ```



