# iQTL pipline for ATAC-seq
## Allelic ATAC-seq mapping for CCGG.

Find the Hi-C pipeline [here](https://github.com/Xieeeee/AlleliC/).
Find the RNA-seq pipeline [here](https://github.com/lindsayhrlee/iQTL_RNA).

### Requirements

- trim_galore v0.6.10
- bowtie2 v2.5.4
- Picard v2..27.5
- Samtools v1.9
- Pysam v0.18.0
- Python v3.6.8
    * Numpy v1.19.5
    * Pandas v1.1.5


### Running the pipeline
We recommend to create a Conda environment that contains all the required software above. This pipeline can be run simply by:

```
./atac_phasing.sh [mouse] [CC1] [CC2] [Name] [prefix_fastq]
```

 The required input variables are:
1. mouse = Mouse Name
2.	CC1 = Maternal Genome Name 
3.	CC2 = Paternal Genome Name 
4.	Name = CC1xCC2_mouse
5.	prefix_fastq = the prefix of the fastq file

Be sure to add/edit the directories in the shellscript (especially the location of picard and on line 86).
### Contact Us

For any questions regarding this software, contact Ming Hu (hum@ccf.org) or Lindsay Lee (leeh7@ccf.org).



