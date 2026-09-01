

Download

  rsync -avP --include='*.vcf.gz' --include='*.vcf.gz.tbi' rsync://ftp.ensembl.org/ensembl/pub/current_variation/vcf/homo_sapiens/  ./vep_data/vcf/

  rsync -avr --progress rsync://ftp.ensembl.org/ensembl/pub/release-115/variation/indexed_vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz  .vep/
  tar -xzf vep_data/.vep/homo_sapiens_vep_115_GRCh38.tar.gz -C vep_data/.vep/


  curl  https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
  mv Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz vep_data/.vep/

Decompress and recompress with bgzip (gzip won't work with VEP)

  gunzip vep_data/.vep/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
  bgzip vep_data/.vep/Homo_sapiens.GRCh38.dna.primary_assembly.fa

Index it

  samtools faidx vep_data/.vep/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

Should say release/115 :

  vep --help 2>&1 | grep "ensembl-vep"


Run vep to get annotated vcf files

  bash scripts/vep_sbatch.sh

Turn into parquet

  python scripts/vcf2parquet.py 