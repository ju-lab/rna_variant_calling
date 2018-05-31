# rna_variant_calling
Calls SNVs from an RNA-seq aligned with STAR using GATK Best Practice
http://www.broadinstitute.org/gatk/guide/article?id=3891

Note that this pipeline does not do indel-realignment nor does BQSR, but can be added if necessary. 


```
$ python rna_variant_calling.py -h
usage: rna_variant_calling.py [-h] -i INPUT_STAR_BAM

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_STAR_BAM, --input_star_bam INPUT_STAR_BAM
```
will generate `std.STAR.sorted.rg.md.split.bam` which can be used as input for GATK HaplotypeCaller.

Final output vcf `std.STAR.sorted.rg.md.split.bam.rnavar.vcf` can be filtered for enough depth by using the following command. 

```bash
bcftools view -i 'DP>10' std.STAR.sorted.rg.md.split.bam.rnavar.vcf
```


## Requirements
* STAR aligned BAM, with 2-pass mode
* PICARD
* GATK v3.5
* Java v1.8
* samtools
