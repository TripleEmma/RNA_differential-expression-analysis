# RNA_differential-expression-analysis
We want to detect genes differentially expressed in two color morphs of species *Gomphocerus sibiricus*.

It is an opportunity for me to learn the principle behind this kind of analysis and get first hand practices.


1. Attain assembly from raw fq data (fq files have been check with FastQC) via Trinity

    1.1 potential problem is that there too many sequences

    https://github.com/trinityrnaseq/trinityrnaseq/wiki
   

2. Quantification via salmon or RSEM

    2.1 results of RSEM and salmon is highly correlated (r = 0.88, p = 0)
   
    https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-323
    

3. Differential expression analysis via DEseq2 (R package)

   https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

   https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8
