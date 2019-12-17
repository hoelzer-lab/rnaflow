https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5982834/
RNA-Seq data are deposited with NCBI GEO under the accession number GSE74809

I just rndm selected some samples. not necessarily replicates, ...

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR293/007/SRR2932637/SRR2932637.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR293/008/SRR2932638/SRR2932638.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR293/009/SRR2932639/SRR2932639.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR293/004/SRR2932644/SRR2932644.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR293/005/SRR2932645/SRR2932645.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR293/006/SRR2932646/SRR2932646.fastq.gz

adding paired-end data

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR504/000/SRR5043290/SRR5043290_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR504/000/SRR5043290/SRR5043290_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR504/001/SRR5043291/SRR5043291_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR504/001/SRR5043291/SRR5043291_2.fastq.gz
seqtk sample -s100 SRR5043290_1.fastq.gz 0.04 > SRR5043290_sub_1.fastq.gz
seqtk sample -s100 SRR5043290_2.fastq.gz 0.04 > SRR5043290_sub_2.fastq.gz
seqtk sample -s100 SRR5043291_1.fastq.gz 0.04 > SRR5043291_sub_1.fastq.gz
seqtk sample -s100 SRR5043291_2.fastq.gz 0.04 > SRR5043291_sub_2.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR268/005/ERR2686025/ERR2686025_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR268/005/ERR2686025/ERR2686025_2.fastq.gz
seqtk sample -s100 ERR2686025_1.fastq.gz 0.03 > ERR2686025_sub_1.fastq.gz
seqtk sample -s100 ERR2686025_2.fastq.gz 0.03 > ERR2686025_sub_2.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR268/004/ERR2686034/ERR2686034_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR268/004/ERR2686034/ERR2686034_2.fastq.gz
seqtk sample -s100 ERR2686034_1.fastq.gz 0.03 > ERR2686034_sub_1.fastq.gz
seqtk sample -s100 ERR2686034_2.fastq.gz 0.03 > ERR2686034_sub_2.fastq.gz