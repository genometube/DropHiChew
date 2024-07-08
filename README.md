# DropHiChew
Droplet based (10x) high efficient chromatin conformation capture

This repository contains all the code that was used to preprocess and analyze DropHiChew data from raw fastq to the tabular data. At the start of each file, dependencies like essential libraries or fixed paths to software and data files are declared, which are required for the scripts to function properly.

For loop velocity analysis in R. Loop velocity can refer to https://github.com/Dekayzc/Loopvelocity


## Speed and resource of DropHiChew data preprocessing
Sample: One lane of MGI T7 PE100 sequencing raw fastq.gz data, 5000 million reads, 7000 cells of HEK293T.

Steps and Performance

### Demultiplexing: Approximately 1.5 days

Resources: 10 threads in parallel (Intel Xeon Platinum 8475B 2.7 GHz), approximately 40GB of memory, and runs on a 2T Alibaba Cloud PL0 ESSD which delivers up to 10,000 IOPS.

### Hi-C Pro Processing: Approximately 3 days

Resources: 20 parallel tasks (Intel Xeon Platinum 8475B 2.7 GHz), each task using 2 threads and 10GB of memory, and runs on Alibaba Cloud NAS.


