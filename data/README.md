# bwa-mem data
* all files are included in docker image, 'schmitzr1984/tcss702-bwa-mem-reference:latest'

### Samples from GATK tutorials
https://console.cloud.google.com/storage/browser/gatk-tutorials/workshop_1903/3-somatic/bams

* normal.bam,

https://storage.cloud.google.com/gatk-tutorials/workshop_1903/3-somatic/bams/normal.bam
* tumor.bam,

https://storage.cloud.google.com/gatk-tutorials/workshop_1903/3-somatic/bams/tumor.bam

* normal and tumor bam files were processed using biobambam2, converting them from bam to fastq
fastq files were then compressed using:
```
tar -cJf normal.tar.xz normal_1.fq normal_2.fq
tar -cJf tumor.tar.xz tumor_1.fq tumor_2.fq
```

### GDC.h38.d1.vd1 BWA Index Files
> **_NOTE:_** the filename is incorrect when downloaded from source, should be GRCh38.d1.vd1_BWA.tar.gz

https://api.gdc.cancer.gov/data/25217ec9-af07-4a17-8db9-101271ee7225


# nlp data
* news_data.tar.xz is derived from,

https://www.kaggle.com/datasets/therohk/million-headlines


# sebs data
* Sample data from serverless benchmark data,

https://github.com/spcl/serverless-benchmarks-data/tree/master/500.scientific/504.dna-visualisation