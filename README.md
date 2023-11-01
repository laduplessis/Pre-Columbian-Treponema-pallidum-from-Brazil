#  Redefining the treponemal history through pre-Columbian genomes from Brazil

**Kerttu Majander** , **Marta Pla-Diaz**, Louis du Plessis, Natasha Arora, Jose Filippini, Luis Pezo Lanfranco, Sabine Eggers, **Fernando González-Candelas**, **Verena J. Schuenemann**

---

This repository contains the data files, configuration files and scripts necesssary to reproduce the molecular clock dating analyses presented in " Redefining the treponemal history through pre-Columbian genomes from Brazil" (Majandar, Pla-Diaz et al). Some scripts may need some adjustments depending on the local setup. 

## Table of contents

1. Abstract
2. Molecular clock dating workflow
  1. Dependencies
  2. Data
  3. Create alignments
  4. BEAST2 analyses
  5. BEAST2 date-randomization test
  6. Post-processing and figures



## Abstract
_The origins of treponemal diseases have long remained enigmatic, especially considering the sudden onset of the first syphilis epidemic in late 15th century Europe, and its hypothesised arrival from the Americas with Columbus’ expeditions. Recently, ancient DNA evidence has revealed various treponemal infections circulating in early modern Europe and colonial-era Mexico3–6. However, no genetic evidence of any treponematosis has been recovered from either the Americas or the Old World that can be reliably dated prior to the first transatlantic contact. Here, we present the first-ever treponemal genomes from nearly 2,000-year-old human remains from Brazil. We reconstruct four ancient genomes of a prehistoric treponemal pathogen, most closely related to the bejel-causing agent, Treponema pallidum ssp. endemicum. Contradicting bejel’s modern-day geographical niche in the arid regions of the world, the results call to question the prior paleopathological characterisation of treponeme subspecies and showcase their adaptive potential. A high-coverage genome is used to improve molecular clock dating estimations, placing the divergence of all modern T. pallidum subspecies firmly in pre-Columbian times. Overall, our study demonstrates the archaeogenetic opportunities to uncover key events in pathogen evolution and emergence, paving the way to unprecedented hypotheses on the spread of treponematoses across time._

---



## Molecular clock dating workflow

### Dependencies

- Python packages:
  - BioPython
  - numpy
  - scipy
  - yaml
- R-packages:
  - tidyverse
  - ggridges
  - cowplot
  - treedataverse
  - coda
  - [beastio](https://github.com/laduplessis/beastio) commit #ff276c2
- snp-sites
- BEAST v2.6.7

### Data
Compressed input alignments and metadata for the molecular clock dating analyses are in [`data/`]().

To create the alignments without ZH1540 for downstream analyses, uncompress fasta files and run: 

```

sed -e '/ZH1540/,+1d' all_01_02_23_alignment_norec_07_02_23.fasta > noZH1540_01_02_23_alignment_norec_07_02_23.fasta
sed -e '/ZH1540/,+1d' all_01_02_23_alignment_norec_no16S23S.fasta > noZH1540_01_02_23_alignment_norec_no16S23S.fasta

```

### Create alignments

Add metadata to sequence headers in alignment files:

```bash

# Alignment without recombining sites
python scripts/processdata.py -i data/Genome_metadata_information.csv -a data/all_01_02_23_alignment_norec_07_02_23.fasta -s "Taxon,Accession number,Subspecies,Year" -o results/alignments/ -p all_01_02_23_alignment_norec

# Alignment without recombining sites, 16S, 23S
python scripts/processdata.py -i data/Genome_metadata_information.csv -a data/all_01_02_23_alignment_norec_no16S23S.fasta -s "Taxon,Accession number,Subspecies,Year" -o results/alignments/ -p all_01_02_23_alignment_norec_no16S23S

# Alignment without recombining sites and ZH1540
python scripts/processdata.py -i data/Genome_metadata_information.csv -a data/noZH1540_01_02_23_alignment_norec_07_02_23.fasta -s "Taxon,Accession number,Subspecies,Year" -o results/alignments/ -p noZH1540_01_02_23_alignment_norec

```

Use `snp-sites` to get the SNP-alignment and constant site breakdown for each dataset:

```bash

# Alignment without recombining sites
snp-sites -o results/alignments/all_01_02_23_snpalignment_norec.fas  results/alignments/all_01_02_23_alignment_norec.fas
snp-sites -C -o results/alignments/all_01_02_23_snpalignment_norec.csv  results/alignments/all_01_02_23_alignment_norec.fas

# Alignment without recombining sites, 16S, 23S
snp-sites -o results/alignments/all_01_02_23_snpalignment_norec_no16S23S.fas  results/alignments/all_01_02_23_alignment_norec_no16S23S.fas
snp-sites -C -o results/alignments/all_01_02_23_snpalignment_norec_no16S23S.csv  results/alignments/all_01_02_23_alignment_norec_no16S23S.fas

# Alignment without recombining sites and ZH1540
snp-sites -o results/alignments/noZH1540_01_02_23_snpalignment_norec.fas  results/alignments/noZH1540_01_02_23_alignment_norec.fas
snp-sites -C -o results/alignments/noZH1540_01_02_23_snpalignment_norec.csv  results/alignments/noZH1540_01_02_23_alignment_norec.fas

```
Only use SNP alignments for subsequent analyses.

Create NEXUS files with SNP sequences, sampling dates and tip date priors.

```bash

# Alignment without recombining sites
python scripts/fasta2nexus.py -i results/alignments/all_01_02_23_snpalignment_norec.fas -o results/alignments/all_01_02_23_snpalignment_norec.nexus -d 3 -c 2

# Alignment without recombining sites, 16S, 23S
python scripts/fasta2nexus.py -i results/alignments/all_01_02_23_snpalignment_norec_no16S23S.fas -o results/alignments/all_01_02_23_snpalignment_norec_no16S23S.nexus -d 3 -c 2

# Alignment without recombining sites and ZH1540
python scripts/fasta2nexus.py -i results/alignments/noZH1540_01_02_23_snpalignment_norec.fas -o results/alignments/noZH1540_01_02_23_snpalignment_norec.nexus -d 3 -c 2

```

### BEAST2 analyses

- Modify NEXUS files by hand to add in extra taxon sets.
- The NEXUS files above can be dragged into BEAUti v2.6 to immediately load the alignment and set tip date priors and taxon sets. 
- Create the XML files in BEAUti v2.6 and adjust by hand:
   - Set to date-forward
	- GTR+G4+I+F
	- UCLD / UCED
	- BSP with 10 groups / Exponential growth / Constant size
	- **Remember to do ascertainment bias correction (manually in xml)!**
	- **Remember to untick monophyly constraints on taxon sets!**
	- **Remember to set maximum population size!**
	- Mean clock rate prior: lognormal with M=1E-7 s/s/y (in real space) and use the same starting value (narrow S=0.25, wide S=1)
- Run [xml]() files in BEAST2.


### BEAST2 date-randomization test

- Create config files in `results/beast2_dateshuffling/config/` by hand. 
- Run Python scripts to produce shuffled XML files:

```bash

python scripts/ShuffleBEASTXML.py -c results/beast2_dateshuffling/config/all_01_02_23_alignment_norec.gtrgi.bsp10.uced.narrow.cfg

```

Run [xml]() files in BEAST2. 


### Post-processing and figures

- Run RMarkdown reports in `reports/`
- Change output to `pdf_document` and `device` to `pdf` or `cairo_pdf` to produce pdf figures.





