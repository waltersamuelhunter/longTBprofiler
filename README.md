
# longTBprofiler

LongTBprofiler is workflow/framework utilizing long read Oxford Nanopore data to classify Mycobacterium tuberculosis isolates into lineages and predict antibiotic resistance. 

## Software requirments 

The analysis is performed in an Unix enviroment (Mac/Linux). The main tools utilized in this pipeline are:
- Nanoplot
- Nanocomp
- tb-profiler
- TbLR


## Installation guide
First we need to install mamba. Mamba is a package manager that will help us easily install all bioinformatic tools for LongTBprofiler pipeline.

### Mamba.

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


#### Packages 

1. Download the 
   
2. Download the txt file ("secret-ingredients") from the provided google drive link and move it to your home directory in ubuntu.
  
3. Create an enviroment containing required packages.

  ```
  mamba create -name longmtbc --file secret-ingredients.txt
  mamba activate longmtbc
  mamba install pip
  pip install tblg
  ```

4. Accept package installations/updates if prompted. 

## Analysis

Now, we have environment ready with all the neccessary packages. Finally we can start analysis. 

Move your fastq files to the home directory.
To check home directory, type:
```
pwd
```



### Preprocessing - quality check, reads removal

*change the sample name if its different than barcode83, 51..

1. Activate mamba environment with all the packages
  ```
  mamba activate longtbprofiler_env
  ```

2.  Generate quality report using nanocomp and nanoplot

  ```

  # for single sample 
  NanoPlot --fastq ./data/WGS/barcode83.fastq.gz -o ./QC/WGS/nanoplot_quality_report
  NanoPlot --fastq ./data/target/barcode51.fastq.gz -o ./QC/target/nanoplot_quality_report

  # for multiple samples
  NanoComp --fastq ./data/WGS/barcode83.fastq.gz ./data/target/barcode51.fastq.gz -o ./QC/WGS/nanocomp_quality_report

  ```

### Perform Read Mapping


   I. Download Reference genome from NCBI

  ```
  mkdir reference
  cd reference
  wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz
  mv GCF_000195955.2_ASM19595v2_genomic.fna.gz reference.fna.gz
  gunzip -d reference.fna
  head reference.fna
  cd ..
  
  ```

  II. Map samples to the reference

  ```
  mkdir mapped_varcall
  ```
  ```
  minimap2 -a -x lr:hq ./data/reference/reference.fna ./data/WGS/barcode83.fastq.gz | samtools sort -O bam -o ./mapped_varcall/WGS/barcode83_sorted.bam 
  ```
  ```
  minimap2 -a -x lr:hq ./data/reference/reference.fna ./data/target/barcode51.fastq.gz | samtools sort -O bam -o ./mapped_varcall/target/barcode51_sorted.bam 
  ```

*if the terminal show you problem with permissions, type:
  ```
  sudo chmod -R a+rwx /path/to/your/directory
  ```
  
 
  III. Samtools

  ``` 
  samtools index ./mapped_varcall/WGS/barcode83_sorted.bam
  ```
  samtools index ./mapped_varcall/target/barcode51_sorted.bam
  ```

### Variant calling

  I. lofreq 

freebayes 

freebayes --ploidy 1-f./data/reference/reference.fna ./mapped_varcall/wgs/ barcoded83_sorted.bam > ./mapped_varcall/wgs/barcoded83.vcf 

freebayes --ploidy 1-f./data/reference/reference.fna ./mapped_varcall/target/barcoded51_sorted.bam > ./mapped_varcall/target/barcoded51.vcf


  ```
  lofreq call -f ./data/reference/reference.fna -o ./mapped_varcall/WGS/barcode83.vcf ./mapped_varcall/WGS/barcode83_sorted.bam
  ```
  ```
  lofreq call -f ./data/reference/reference.fna -o ./mapped_varcall/target/barcode51.vcf ./mapped_varcall/target/barcode51_sorted.bam
  ```



### Lineage drug resistance typing 

  Prepare directories for typing results
  ```
  mkdir TB_Profiler
  cd TB_Profiler 
  mkdir 



  II. Tblg

  Create a new environemnt for Tuberculosis Lineage Genotyping (TbLG)

  ```
  mamba create -n tblg_env python pip
  pip install tblg
  mamba activate tblg_env
  pip install tblg
  ```


  Run the lineage typing analysis.

  ```
  tblg --o 
  ```

  III. TB-profiler

  Create a new environment for TB-profiler

  ```
  mamba create -n tbprofiler_env tb-profiler
  ```
  ```
  tb-profiler
  ```
