import pandas as pd
import os
import sys
import subprocess
from tqdm.notebook import tqdm
from Bio  import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import shutil
import time
from requests.exceptions import HTTPError
import glob

import geneinfo as gi
gi.email("kaspermunch@birc.au.dk")

from tqdm import tqdm

import Bio.Data.CodonTable
Bio.Data.CodonTable.standard_dna_table

from Bio.Data.CodonTable import unambiguous_dna_by_id


def altcodons(codon, table):
    """
    List codons that code for the same aminonacid / are also stop.
    """
    tab = Bio.Data.CodonTable.standard_dna_table
#    tab = unambiguous_dna_by_id[table]

    if codon in tab.stop_codons:
        return tab.stop_codons

    try:
        aa = tab.forward_table[codon]
    except:
        return []

    return [
        k
        for (k, v) in tab.forward_table.items()
        if v == aa and k[0] == codon[0] and k[1] == codon[1]
    ]


def degeneration(codon, table):
    """
    Determine how many codons code for the same amino acid / are also stop
    """
    return len(altcodons(codon, table))


def is_x_degenerated(x, codon, table):
    """
    Determine if codon is x-fold degenerated.
    """
    alt = altcodons(codon, table)
    third = {c[-1] for c in alt}
    return x <= len(alt) and x <= len(third)


def degenerated_subseq(seq, x, table):
    """
    Get a subsequence consisting of the x-fold degenerated codons only.
    """
    data = ""
    for i in range(0, len(seq), 3):
        codon = seq[i : i + 3]#.tostring()
        if is_x_degenerated(x, codon, table):
            data += codon
    return data

def fourfold_degenerate_GC_content(seq):
    """
    Get GC content third position of of four-fold degenerate codons
    (with alternative codons allowing all three bases at the third position)
    """
    four_fold_deg_codons = degenerated_subseq(seq, 4, 'Standard')
    thirds = four_fold_deg_codons[2:len(four_fold_deg_codons):3]
    return (thirds.count('G') + thirds.count('C')) / len(thirds), len(thirds)


def get_orf(geneid):

    # clean up previous call
    if os.path.exists('README.md'):
        os.remove('README.md')
    if os.path.exists('md5sum.txt'):
        os.remove('md5sum.txt')        
    if os.path.exists('ncbi_dataset.zip'):
        os.remove('ncbi_dataset.zip')
    if os.path.exists('ncbi_dataset'):        
        shutil.rmtree('ncbi_dataset')
    time.sleep(1)

    # download rna and protein
    subprocess.call(f'datasets download gene gene-id {geneid}'.split(),
                     stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
                     ) 
    time.sleep(1)
    # unzip
    subprocess.call(f'unzip ncbi_dataset.zip'.split(), 
                    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
                    )
    time.sleep(1)
    
    # read protein seq
    with open("ncbi_dataset/data/protein.faa") as f:
        protein_record = next(SeqIO.parse(f, "fasta"))
    orf = None
    # read rna transcripts
    with open("ncbi_dataset/data/rna.fna") as f:
        # loop transcripts to find the one corresponding to the protein sequence
        for rna_record in (SeqIO.parse(f, "fasta")):
            # loop reading frames
            for frame in range(3):
                rna = str(rna_record[frame:].seq)
                # pad with N to make length multiplum of 3
                rna = rna + 'N' * (3 - len(rna) % 3)
                # see if translated rna is part of protein sequence
                trans = str(SeqRecord(Seq(rna), id='rna').translate().seq)
                idx = trans.find(str(protein_record.seq))
                if idx >= 0:
                    # start of orf in rna
                    offset = frame + idx * 3
                    # cut out orf
                    orf = str(rna[offset:len(protein_record.seq) * 3 + offset])
                    break
            if orf is not None:
                break

    return orf


if __name__ == "__main__":
        
    orthlogs = pd.read_csv(sys.argv[1])
    output_file = sys.argv[2]

    if not os.path.exists('tmp'):
        os.makedirs('tmp')

    for row in tqdm(orthlogs.itertuples(), total=len(orthlogs)):

        out_file_name = f'tmp/human_ncbi_{row.HsEntrez}.tsv'
        if os.path.exists(out_file_name):
            continue

        tries = 3
        # try op to three times in case there is some hickup
        while tries:
            try:
                human_ncbi_id = row.HsEntrez
                human_orf = get_orf(human_ncbi_id)
                if human_orf is None:
                    print('skipping human', human_ncbi_id, file=sys.stderr)
                    break
                yeast_ncbi_id = gi.ensembl2ncbi(row.ScENSP)
                yeast_orf = get_orf(yeast_ncbi_id)    
                if yeast_orf is None:
                    print('skipping yeast', yeast_ncbi_id, file=sys.stderr)
                    break
            except HTTPError as e:
                print(e, file=sys.stderr)
                tries -= 1
                continue
            except FileNotFoundError as e:
                print(e, file=sys.stderr)
                break
            except Exception as e:
                print(e, file=sys.stderr)
                tries -= 1
                continue

            # write line to a file for the gene pair
            with open(out_file_name, 'w') as f:
                print(human_ncbi_id, *fourfold_degenerate_GC_content(human_orf), 
                        yeast_ncbi_id, *fourfold_degenerate_GC_content(yeast_orf), sep='\t', file=f)

            break

    with open(output_file, 'w') as out:
        for path in glob.glob('tmp/*.tsv'):
            with open(path) as f:
                for line in f:
                    print(line.strip(), file=out)
