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

#### Mamba.

I. Go to: https://conda-forge.org/download/
II. choose the correct version for your system (x86_64/arm64), then right click and "copy link adress"
III. open terminal, go to home directory (with command "cd")
IV. type "wget", paste the copied link and then enter
IV. After a moment u have your miniforge downloaded. Check the name of file by viewing the content of current directory with "ls"
V. Copy the name of the downloaded file. It should be "Miniforge3-Linux-x86_64.sh " or something similar.
VI. Write

  ```
  bash "Miniforge3-Linux-x86_64.sh"

  ```





