# snp-annotation
Contact: m.lundberg@uq.net.au


## Prerequisites:
- Script was tested with the following versions installed:
 - Anaconda2
 - Pandas: 0.22.0
 - Numpy: 1.14.3
 - Pysam: 0.14.1
 
## Getting Started
In order to download the SNP annotation tool, you should clone this repository via the commands
```  
git clone git@github.com:MischaLundberg/snp-annotation.git
cd snp-annotation
```

In order to install the Python dependencies, you will need the [Anaconda](https://store.continuum.io/cshop/anaconda/) Python distribution and package manager. After installing Anaconda, run the following commands to create an environment with TEpoly's dependencies:

```
conda env create --file environment.yml
source activate snpannotation

## to deactivate the environment, type
#conda deactivate
```

In case you are updating your current version of the CpG Methylation pipeline, it would be best practice to also update your environment to the updated prerequisetes.
If you are using a Anaconda environment, you can do so by typing
```
conda env update --name snpannotation --file environment.yml
```


If you receive any errors while running the SNP annotation tool, please ensure your versioning for the prerequisites is according to the tested versions.