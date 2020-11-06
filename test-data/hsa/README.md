
```bash
# MAQCA_1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR367/007/SRR3670977/SRR3670977.fastq.gz
seqtk sample -s100 SRR3670977.fastq.gz 100000 > SRR3670977_sub.fastq
rm SRR3670977.fastq.gz
# MAQCA_1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR367/008/SRR3670978/SRR3670978.fastq.gz
seqtk sample -s100 SRR3670978.fastq.gz 100000 > SRR3670978_sub.fastq
rm SRR3670978.fastq.gz
# MAQCB_1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR367/005/SRR3670985/SRR3670985.fastq.gz
seqtk sample -s100 SRR3670985.fastq.gz 100000 > SRR3670985_sub.fastq
rm SRR3670985.fastq.gz
# MAQCB_1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR367/006/SRR3670986/SRR3670986.fastq.gz
seqtk sample -s100 SRR3670986.fastq.gz 100000 > SRR3670986_sub.fastq
rm SRR3670986.fastq.gz
```
