from gwf import Workflow, AnonymousTarget
import os
import pandas as pd
import glob

gwf = Workflow()

starting_dir = os.path.expanduser("~/primatediversity/data/gVCFs_recalling_10_12_2024/{}/filteredVCF/bcf_step1/")
metadata_path = os.path.expanduser("~/primatediversity/data/gVCFs_recalling_10_12_2024_metadata/")
zarr_dir = "zarr_data/"

def generate_zarr(focus_chr, all_chr, out_path):
    inputs = [all_chr]
    o = out_path+focus_chr+"/call_genotype"
    outputs = [o]
    options = {
        "cores": 4,
        "memory": "30g",
        "walltime": "120:00:00",
        "account": "baboondiversity"
    }
    spec = """
    bcftools view -Ou -r {focus_chr} {all_chr} > {temp_bcf}
    bcftools index {temp_bcf}
    vcf2zarr convert {temp_bcf} {out_path}
    rm {temp_bcf}
    rm {temp_bcf}.csi
    """.format(focus_chr=focus_chr, all_chr=all_chr,
               temp_bcf=out_path+focus_chr+"_temp.bcf", out_path=out_path+focus_chr)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def get_ID_vcf(idx, target):
    filename = target.spec.split("/")[-2].split(".")[0]+"_"+target.spec.split("/")[-1].split(".")[0]
    return 'vcf_zarr_{}'.format(filename)

# This section goes through all metadata folders and identifies the corresponding VCFs to use.
metadata_folders = glob.glob(metadata_path+"*_individuals.txt")

for folder in metadata_folders:
    metadata_df = pd.read_csv(folder, sep="\t")
    short_form = folder.split("/")[-1].split("_")[0]
    regions_df = pd.read_csv(metadata_path+"{}_regions_and_batches.txt".format(short_form), sep="\t")
    for GVCF_FOLDER in metadata_df.GVCF_FOLDER.unique():
        #print(GVCF_FOLDER)
        reference = metadata_df.loc[metadata_df.GVCF_FOLDER == GVCF_FOLDER].REFERENCE_FOLDER.unique()
        all_chr_file = starting_dir.format(GVCF_FOLDER)+"{}_all_chr.sorted.bcf".format(GVCF_FOLDER)
        if not os.path.exists(all_chr_file):
            print("Problem with ", GVCF_FOLDER)
            continue
        zarr_input = []
        # Chromosome X
        zarr_input.extend(regions_df.loc[(regions_df.END >= 1000000) &
                                      (regions_df.FEMALE_PLOIDY == 2) &
                                      (regions_df.MALE_PLOIDY == 1) &
              (regions_df.REFERENCE_FOLDER == reference[0])].CONTIG_ID.unique())
        # Autosomes
        zarr_input.extend(regions_df.loc[(regions_df.END >= 1000000) &
                                      (regions_df.FEMALE_PLOIDY == 2) &
                                      (regions_df.MALE_PLOIDY == 2) &
              (regions_df.REFERENCE_FOLDER == reference[0])].CONTIG_ID.unique()[:1])
        out_path = zarr_dir+GVCF_FOLDER+"/"
        os.makedirs(out_path, exist_ok=True)
        gwf.map(generate_zarr, zarr_input, name=get_ID_vcf,
                extra={"all_chr": all_chr_file, "out_path": out_path})
