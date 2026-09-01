

Download

  rsync -avP --include='*.vcf.gz' --include='*.vcf.gz.tbi' rsync://ftp.ensembl.org/ensembl/pub/current_variation/vcf/homo_sapiens/  ./vcf/

  rsync -avr --progress rsync://ftp.ensembl.org/ensembl/pub/release-115/variation/indexed_vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz  .vep/
  tar -xzf .vep/homo_sapiens_vep_115_GRCh38.tar.gz -C .vep/


  curl  https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
  mv Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz .vep/

Decompress and recompress with bgzip (gzip won't work with VEP)

  gunzip .vep/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
  bgzip .vep/Homo_sapiens.GRCh38.dna.primary_assembly.fa

Index it

  samtools faidx .vep/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

Should say release/115 :

  vep --help 2>&1 | grep "ensembl-vep"


Run vep to get annotated vcf files

  sbatch vep_sbatch.sh

Turn into parquet

  python vcf2parquet.py 
