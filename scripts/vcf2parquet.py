from pathlib import Path
import re, cyvcf2, pyarrow as pa, pyarrow.parquet as pq
from tqdm.auto import tqdm

PARQUET_DIR = Path('parquet')

def get_csq_fields(vcf: cyvcf2.VCF) -> list[str]:
    for line in str(vcf.raw_header).split('\n'):
        if 'ID=CSQ' in line:
            match = re.search(r'Format: ([^"]+)', line)
            if match:
                return match.group(1).strip().split('|')
    return []

for chrom in tqdm([f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']):
    vcf = cyvcf2.VCF(f'vcf_annotated/homo_sapiens-{chrom}.vep.vcf.gz')
    csq_fields = get_csq_fields(vcf)

    rows = []
    for v in vcf:
        csq_raw = v.INFO.get('CSQ', '')
        if not csq_raw:
            continue
        for transcript in csq_raw.split(','):
            record = dict(zip(csq_fields, transcript.split('|')))
            # Skip non-genic consequences
            if not record.get('SYMBOL'):
                continue
            rows.append({
                'chrom'            : v.CHROM,
                'pos'              : v.POS,
                'ref'              : v.REF,
                'alt'              : ','.join(v.ALT),
                'variant_id'       : v.ID or '.',
                'gene_symbol'      : record.get('SYMBOL', ''),
                'gene_id'          : record.get('Gene', ''),
                'transcript_id'    : record.get('Feature', ''),
                'consequence_terms': record.get('Consequence', ''),
                'impact'           : record.get('IMPACT', ''),
                'biotype'          : record.get('BIOTYPE', ''),
                'canonical'        : 1 if record.get('CANONICAL') == 'YES' else 0,
                'hgvsc'            : record.get('HGVSc', ''),
                'hgvsp'            : record.get('HGVSp', ''),
                'sift'             : record.get('SIFT', ''),
                'polyphen'         : record.get('PolyPhen', ''),
                'af_gnomad'        : record.get('gnomADe_AF', '') or record.get('AF', ''),
                'clin_sig'         : record.get('CLIN_SIG', ''),
            })

    if not rows:
        continue

    pq.write_to_dataset(
        pa.Table.from_pylist(rows),
        root_path=PARQUET_DIR,
        partition_cols=['chrom'],
        compression='zstd',
        existing_data_behavior='overwrite_or_ignore',
    )