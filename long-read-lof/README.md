# Mao et al. 2024 (8 NHPs) — Data S1–S5 incl. the disrupted-gene table (Data S2)
wget -O mao2024_suppl.zip \
  "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10947866/supplementaryFiles"

# Yoo et al. 2025 (6 apes) — Supplementary Tables incl. VIII.34 (gain/dup/loss)
wget -O yoo2025_suppl.zip \
  "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12058530/supplementaryFiles"

unzip -o mao2024_suppl.zip -d mao2024_suppl
unzip -o yoo2025_suppl.zip -d yoo2025_suppl

BASE="https://genome.senckenberg.de/download/TOGA/human_hg38_reference/Primates"
for sp in Pan_troglodytes Pan_paniscus Gorilla_gorilla Pongo_abelii \
          Pongo_pygmaeus Symphalangus_syndactylus Macaca_mulatta; do
  wget -r -np -nH --cut-dirs=4 -A "loss_summ_data.tsv.gz,orthologyClassification.tsv.gz" \
       "$BASE/$sp*/"
done
