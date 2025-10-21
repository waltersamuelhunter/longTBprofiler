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


### Installation guide
First we need to install mamba. Mamba is a package manager that will help us easily install all bioinformatic tools for LongTBprofiler pipeline.

#### Mamba.

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






