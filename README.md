# PNS-glia

## System requirement
The code has been implemented and tested on both MacOS 11.6.1 (Big Sur)  and Ubuntu/Linux 4.15.0-46-generic kernel.

For R codes, the following library version were tested:\
R 4.0.2\
Seurat 3.2.3 & 4.0.4\
ggplot2 3.3.2\
cowplot 1.1.0\
Matrix 1.3-4\
monocle3 0.2.3.0\
edgeR 3.20.9\
limma 3.34.9\

For Python, the following version has been tested:\
Python 3.7 & 3.8

Libraries required by the Python code:\
pysam\
pandas\
seaborn\
matplotlib\
numpy\
scipy\
tqdm\
joblib

For the Circos plot, circos-0.69-6 was used.

## Installation guide

User will have to install R and Python3 with the specific version number stated above. RStudio is recommended for running the R codes. The Python libraries can be installed through the Python Package Installer pip3. No compilation is required for both R and Python codes and can be run directly.

For all python codes, the arguments will be shown when executed, for example : 
Bash - python3 module-based_analysis.py

usage: module-based_analysis.py [-h] [-s SCRNA] [-c CLUSTER] [-g GTF]

optional arguments:
  -h, --help            show this help message and exit
  -s SCRNA, --scRNA SCRNA
                        ESSENTIAL Single cell RNA-seq matrix
  -c CLUSTER, --cluster CLUSTER
                        ESSENTIAL Cluster identity from Seurat
  -g GTF, --gtf GTF     ESSENTIAL GTF file

