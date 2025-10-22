# longTBprofiler

LongTBprofiler is workflow/framework utilizing long read Oxford Nanopore data to classify Mycobacterium tuberculosis isolates into lineages and predict antibiotic resistance. 

## Software requirments 

The analysis is performed in an Unix enviroment (Mac/Linux). The required tools for this pipeline are:
- Nanoplot
- Nanocomp
- minimap2
- samtools
- ivar
- mtbtyper (with dependencies) 
- TB-Profiler


## Installation guide
First we need to install mamba. Mamba is a package manager that will help us easily install all bioinformatic tools for LongTBprofiler pipeline.

### Mamba.

I. Go to: https://conda-forge.org/download/

II. Choose the correct version for your system (x86_64/arm64), then right click and "copy link adress"

III. Go to home directory in terminal (with command "cd")

  ```
  cd

  ```

IV. Download the most recent miniforge version, suitable for your software

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

It stould print the version of mamba, show the location of installation and list options. If you see these, your mamba has been installed correctly and u can proceed. If not, please inform one of the instructors. 


#### Packages 

1. Copy the content of "packages.yaml" file from this github repository.
2. In the terminal, create a packages.yaml file and paste the copied content.

  ```
  nano packages.yaml
  ```
3. Create an enviroment containing packages listed in the yaml file

  ```
  mamba create -f packages.yaml
  ```
4. Accept package installations/updates if prompted. 

Note that it might take several minutes and u might not see any logs in the meantime.


## Analysis

Create a directory for an analysis.

  ```
  mkdir tutorial_tb_longreads
  cd
  ```
Move your fastq files to this directory.
In case you don't have your own data and what to try this pipeline, you can download sample from NCBI SRA database. However, you need to have sra-tool kit installed. 

  ``` # downloading example sample
  prefetch SRR35794931
  find . -mindepth 2 -type f -name "*.sra" -exec mv {} . \;
  fasterq-dump SRR35794931.sra
  gzip SRR35794931.fastq 
  ```

### Preprocessing - quality check, reads removal

*change the sample name if its different than barcode83, 51..

1. Activate mamba environment with all the packages
  ```
  mamba activate longtbprofiler_env
  ```

2.  Generate quality report using nanocomp and nanoplot

  ```

  NanoComp --fastq ./data/WGS/barcode83.fastq.gz -o ./QC/WGS/nanocomp_quality_report
  NanoComp --fastq ./data/target/barcode51.fastq.gz -o ./QC/target/nanocomp_quality_report

  NanoPlot --fastq ./data/WGS/barcode83.fastq.gz -o ./QC/WGS/nanoplot_quality_report
  NanoPlot --fastq ./data/target/barcode51.fastq.gz -o ./QC/target/nanoplot_quality_report

  ```

3. Perform Read Mapping

   I. Download Reference genome from NCBI

  ```
  mkdir reference
  cd reference
  wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz
  mv GCF_000195955.2_ASM19595v2_genomic.fna.gz reference.fna.gz
  gunzip -d reference.fna
  cat reference.fna
  cd ..
  
  ```

  II. Map samples to the reference

  ```
  mkdir mapped
  minimap2 -a -x lr:hq ./data/reference/reference.fna ./data/WGS/barcode83.fastq.gz | samtools sort -O bam -o ./mapped_varcall/WGS/barcode83_sorted.bam -
  minimap2 -a -x lr:hq ./data/reference/reference.fna ./data/target/barcode51.fastq.gz | samtools sort -O bam -o ./mapped_varcall/target/barcode51_sorted.bam -

  ```
 
  III. Samtools

  ``` 
  samtools index ./mapped_varcall/WGS/barcode83_sorted.bam
  samtools index ./mapped_varcall/target/barcode51_sorted.bam

  ```

  IV. lofreq 

  ```
  lofreq call -f ./data/reference/reference.fna -o ./mapped_varcall/WGS/barcode83.vcf ./mapped_varcall/WGS/barcode83_sorted.bam
  lofreq call -f ./data/reference/reference.fna -o ./mapped_varcall/target/barcode51.vcf ./mapped_varcall/target/barcode51_sorted.bam


  lofreq call -f ref.fa -o vars.vcf aln.bam


  V. Mtbtyper

  VI. Tblg

  VII. TB-profiler 
