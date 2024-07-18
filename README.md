# Summer 2024 UG research project
# Ctenid Tx toxins evolution and protein sturcture analyses

### Michael S. Brewer
### Jordan Baccus
### Amir El

## Dependencies

- Ubuntu 22.04.4 LTS (GNU/Linux 5.15.153.1-microsoft-standard-WSL2 x86_64)
- HYPHY 2.5.62(MP) for Linux on x86_64 x86 SSE4 SIMD zlib (v1.2.13)
- MAFFT v7.525 (2024/Mar/13)
- pal2nal.pl v14
- IQ-TREE multicore version 2.3.4 COVID-edition for Linux x86 64-bit built Jun 18 2024
- biopython v1.78

### Install from conda

```
conda install -c bioconda -c conda-forge mafft pal2nal iqtree biopython hyphy wkhtmltopdf pandoc
```

### Install from GitHub

```
wget https://raw.githubusercontent.com/jvollme/fasta2phylip/master/fasta2phylip.py
```

## Sequence collection and curration

### Download ickFinderCDS.correctedCDS.csv from tijeco/killerknots GitHub

```
mkdir ctenid_sequences

cd ctenid_sequences
```

### Extract non-male ICK entries for only ctenid taxa

```
grep -v "_male_" ../ickFinderCDS.correctedCDS.csv | grep "^Anahita\|^Ctenus\|^Isoctenus\|^Leptoctenus\|^Phoneutria" > ickFinderCDS.correctedCDS_Ctenids.csv
```

### Create fasta files for CDS, full peptides, and mature peptides

- CDS

```
cut -d"," -f1,2 ickFinderCDS.correctedCDS_Ctenids.csv | sed 's/^/>/g' | sed 's/,/\n/g' > ickFinderCDS.correctedCDS_Ctenids_CDS.fasta
```

- Mature peptides

```
cut -d"," -f1,5 ickFinderCDS.correctedCDS_Ctenids.csv | sed 's/^/>/g' | sed 's/,/\n/g' > ickFinderCDS.correctedCDS_Ctenids_MaturePeptides.fasta
```

- Full peptides

```
cut -d"," -f1,6 ickFinderCDS.correctedCDS_Ctenids.csv | sed 's/^/>/g' | sed 's/,/\n/g' > ickFinderCDS.correctedCDS_Ctenids_FullPeptides.fasta
```

#### Download peptide sequences for PnTx1, PnTx2, PnTx3, PnTx4 from Genbank

| PnTx family | PnTx ID | Accession Number |
| ----------- | ------- | ---------------- |
| PnTx1       | PnTx1   | P17727.2         |
| PnTx2       | PnTx2-6 | P29425.2         |
| PnTx3       | PnTx3-6 | P81792.2         |
| PnTx4       | PnTx4-5 | P59367.2         |

> \>sp|P17727.2|TXL1_PHONI RecName: Full=Mu-ctenitoxin-Pn1a; Short=Mu-CNTX-Pn1a; AltName: Full=Toxin > Tx1; Short=PNTx1; Short=PhTx1; Flags: Precursor
> MKLLGIFLVASFAFVLSFGEEMIEGENPLEDQRAELTSCFPVGHECDGDASNCNCCGDDVYCGCGWGRWN
> CKCKVADQSYAYGICKDKVNCPNRHLWPAKVCKKPCRRNCGG
>  
> \>sp|P29425.2|TX36A_PHONI RecName: Full=Delta-ctenitoxin-Pn2a; Short=Delta-CNTX-Pn2a; AltName: > Full=Delta-CNTX-Pn1c; AltName: Full=Neurotoxin Tx2-6; Short=PnTx2-6; Flags: Precursor
> MKVAILFLSILVLAVASESIEESRDDFAVEELGRATCAGQDQPCKETCDCCGERGECVCGGPCICRQGYF
> WIAWYKLANCKK
>  
> \>sp|P81792.2|TX90B_PHONI RecName: Full=Omega-ctenitoxin-Pn4a; Short=Omega-CNTX-Pn4a; AltName: > Full=CTK 01512-2; AltName: Full=Neurotoxin Tx3-6; Short=PnTx3-6; AltName: Full=Ph-alpha-1-beta > toxin; Flags: Precursor
> MKCAVLFLSVIALVHIFVVEAEEEPDSDALVPQERACIPRGEICTDDCECCGCDNQCYCPPGSSLGIFKC
> SCAHANKYFCNRKKEKCKKA
>  
> \>sp|P59367.2|TX35C_PHONI RecName: Full=GAMMA-ctenitoxin-Pn1a; Short=GAMMA-CNTX-Pn1a; AltName: > Full=Insecticidal neurotoxin Tx4(5-5); Short=PnTx4(5-5); Short=PnTx4-5-5; AltName: Full=Toxin Pn4A; > Flags: Precursor
> MKVAIVFLSLLVLAFASESIEENREEFPVEESARCADINGACKSDCDCCGDSVTCDCYWSDSCKCRESNF
> KIGMAIRKKFC

- Add these to new full peptide file and reduce names to "Full=" names only

```
cp ickFinderCDS.correctedCDS_Ctenids_FullPeptides.fasta ickFinderCDS.correctedCDS_Ctenids_FullPeptides_withPnTxFromGenbank.fasta
```
- Add GenBank sequences

## Create phylogeny of Full peptides including GenBank sequences

### Multiple sequence alignment

```
mafft ickFinderCDS.correctedCDS_Ctenids_FullPeptides_withPnTxFromGenbank.fasta > ickFinderCDS.correctedCDS_Ctenids_FullPeptides_withPnTxFromGenbank_MAFFT.aln.fasta
```

- Convert FASTA to PHYLIP

```
../fasta2phylip.py -i ickFinderCDS.correctedCDS_Ctenids_FullPeptides_withPnTxFromGenbank_MAFFT.aln.fasta -o ickFinderCDS.correctedCDS_Ctenids_FullPeptides_withPnTxFromGenbank_MAFFT.aln.phy -r
```

### Phylogeny inference using IQTREE2

#### Unrooted tree

```
iqtree -s ickFinderCDS.correctedCDS_Ctenids_FullPeptides_withPnTxFromGenbank_MAFFT.aln.phy -T AUTO -B 1000 -m MFP
```

#### Tree rooted using non-reversible model

```
iqtree -s ickFinderCDS.correctedCDS_Ctenids_FullPeptides_withPnTxFromGenbank_MAFFT.aln.phy -B 1000  --model-joint NONREV -T AUTO --prefix nonrev_aa
```

#### Rootstrapping to explore support for root position

```
iqtree -s ickFinderCDS.correctedCDS_Ctenids_FullPeptides_withPnTxFromGenbank_MAFFT.aln.phy --model-joint NONREV --root-test -zb 1000 -au -te nonrev_aa.treefile --prefix nonrev_aa_test
```

## Create phylogeny of CDS sequences

### Multiple sequence alignment

- Full peptides without GenBank sequences

```
 mafft ickFinderCDS.correctedCDS_Ctenids_FullPeptides.fasta > ickFinderCDS.correctedCDS_Ctenids_FullPeptides_MAFFT.aln.fasta
```

- Backtranslate CDS sequences onto Full peptide alignment

```
pal2nal.pl ickFinderCDS.correctedCDS_Ctenids_FullPeptides_MAFFT.aln.fasta ickFinderCDS.correctedCDS_Ctenids_CDS.fasta -output fasta > ickFinderCDS.correctedCDS_Ctenids_CDS_PAL2NAL.aln.fasta
```

- Convert FASTA to PHYLIP

```
../fasta2phylip.py -i ickFinderCDS.correctedCDS_Ctenids_CDS_PAL2NAL.aln.fasta -o ickFinderCDS.correctedCDS_Ctenids_CDS_PAL2NAL.aln.phy -r
```

### Phylogeny inference using IQTREE2

#### Unrooted tree

```
iqtree -s ickFinderCDS.correctedCDS_Ctenids_CDS_PAL2NAL.aln.phy -T AUTO -B 1000 -st CODON1 --modelomatic
```

#### Tree rooted using non-reversible model

```
iqtree -s ickFinderCDS.correctedCDS_Ctenids_CDS_PAL2NAL.aln.phy -B 1000 --model-joint UNREST -T AUTO --prefix nonrev_dna
```

#### Rootstrapping to explore support for root position

```
iqtree -s ickFinderCDS.correctedCDS_Ctenids_CDS_PAL2NAL.aln.phy --model-joint UNREST --root-test -zb 1000 -au -te nonrev_dna.treefile --prefix nonrev_dna_test
```

## PnTx3 analyses

- Explore Full peptide phylogeny and extract clade associated with PnTx3
    - We chose 17 sequences
    - Place terminal names in text file (PnTx3.fasta.txt)

### Extract PnTx3 associated sequences from CDS, Full Peptide, and Mature Peptides

```
python ../ToxinSorter.py PnTx3.fasta.txt ickFinderCDS.correctedCDS_Ctenids_CDS.fasta ickFinderCDS.correctedCDS_Ctenids_FullPeptides.fasta ickFinderCDS.correctedCDS_Ctenids_MaturePeptides.fasta
```

### PnTx3 phylogenetic analyses

#### CDS sequences

##### Multiple sequence alignment

- Full peptides

```
 mafft ickFinderCDS.correctedCDS_Ctenids_FullPeptides.fastaToxinSorterout.fasta > ickFinderCDS.correctedCDS_Ctenids_FullPeptides_ToxinSorterout.MAFFT.aln.fasta
```

- Backtranslate CDS sequences onto Full peptide alignment

```
pal2nal.pl ickFinderCDS.correctedCDS_Ctenids_FullPeptides_ToxinSorterout.MAFFT.aln.fasta ickFinderCDS.correctedCDS_Ctenids_CDS.fastaToxinSorterout.fasta -output fasta > ickFinderCDS.correctedCDS_Ctenids_FullPeptides_ToxinSorterout.PAL2NAL.aln.fasta
```

- Convert FASTA to PHYLIP

```
../fasta2phylip.py -i ickFinderCDS.correctedCDS_Ctenids_FullPeptides_ToxinSorterout.PAL2NAL.aln.fasta -o ickFinderCDS.correctedCDS_Ctenids_FullPeptides_ToxinSorterout.PAL2NAL.aln.phy -r
```

##### Phylogeny inference using IQTREE2

###### Unrooted tree

```
iqtree -s ickFinderCDS.correctedCDS_Ctenids_FullPeptides_ToxinSorterout.PAL2NAL.aln.phy -T AUTO -B 1000 -st CODON1 --modelomatic
```

#### CDS sequences without _Ctenus corniger_

- Full peptides

```
 mafft ickFinderCDS.correctedCDS_Ctenids_FullPeptides_ToxinSorterout_noCtenuscorniger.fasta > ickFinderCDS.correctedCDS_Ctenids_FullPeptides_ToxinSorterout_noCtenuscorniger.MAFFT.aln.fasta
```

- Backtranslate CDS sequences onto Full peptide alignment

```
pal2nal.pl ickFinderCDS.correctedCDS_Ctenids_FullPeptides_ToxinSorterout_noCtenuscorniger.MAFFT.aln.fasta ickFinderCDS.correctedCDS_Ctenids_CDS.fastaToxinSorterout_noCtenuscorniger.fasta -output fasta > ickFinderCDS.correctedCDS_Ctenids_FullPeptides_ToxinSorterout_noCtenuscorniger.PAL2NAL.aln.fasta
```

- Convert FASTA to PHYLIP

```
../fasta2phylip.py -i ickFinderCDS.correctedCDS_Ctenids_FullPeptides_ToxinSorterout_noCtenuscorniger.PAL2NAL.aln.fasta -o ickFinderCDS.correctedCDS_Ctenids_FullPeptides_ToxinSorterout_noCtenuscorniger.PAL2NAL.aln.phy -r
```

##### Phylogeny inference using IQTREE2

###### Unrooted tree

```
iqtree -s ickFinderCDS.correctedCDS_Ctenids_FullPeptides_ToxinSorterout_noCtenuscorniger.PAL2NAL.aln.phy -T AUTO -B 1000 -st CODON1 --modelomatic
```

###### Tree rooted using non-reversible model

```
iqtree -s ickFinderCDS.correctedCDS_Ctenids_FullPeptides_ToxinSorterout_noCtenuscorniger.PAL2NAL.aln.phy -B 1000 --model-joint UNREST -T AUTO --prefix  nonrev_dna_ToxinSorterout_noCtenusconiger
```

###### Rootstrapping to explore support for root position

```
iqtree -s ickFinderCDS.correctedCDS_Ctenids_FullPeptides_ToxinSorterout_noCtenuscorniger.PAL2NAL.aln.phy --model-joint UNREST --root-test -zb 1000 -au -te nonrev_dna.treefile --prefix nonrev_dna_test_ToxinSorterout_noCtenusconiger
```

##### Tests for Selection

###### HyPhy - FEL

```
hyphy
```
- Options
    - 1 (selection)
    - 2 (FEL)
    - ickFinderCDS.correctedCDS_Ctenids_FullPeptides_ToxinSorterout_noCtenuscorniger.PAL2NAL.aln.phy (alignment)
    - ickFinderCDS.correctedCDS_Ctenids_FullPeptides_ToxinSorterout_noCtenuscorniger.PAL2NAL.aln.phy.treefile (treefile)

|     Codon      |   Partition    |     alpha      |      beta      |      LRT       |Selection detected?|
|:--------------:|:--------------:|:--------------:|:--------------:|:--------------:|:-----------------:|
|       4        |       1        |        9.115   |        0.548   |        5.367   |  Neg. p = 0.0205  |
|       6        |       1        |        8.129   |        0.530   |        3.686   |  Neg. p = 0.0549  |
|       15       |       1        |       12.016   |        0.000   |        3.904   |  Neg. p = 0.0482  |
|       17       |       1        |       25.560   |        0.323   |        8.762   |  Neg. p = 0.0031  |
|       35       |       1        |        3.667   |        0.000   |        6.485   |  Neg. p = 0.0109  |
|       40       |       1        |        5.715   |        0.748   |        2.763   |  Neg. p = 0.0964  |
|       41       |       1        |        0.000   |        3.229   |        2.745   |  Pos. p = 0.0975  |
|       44       |       1        |       12.096   |        1.330   |        2.773   |  Neg. p = 0.0959  |
|       49       |       1        |       12.508   |        0.000   |       11.044   |  Neg. p = 0.0009  |
|       50       |       1        |        1.320   |        0.000   |        2.940   |  Neg. p = 0.0864  |
|       54       |       1        |       25.749   |        0.442   |       13.288   |  Neg. p = 0.0003  |
|       61       |       1        |        1.275   |        0.000   |        2.895   |  Neg. p = 0.0888  |
|       62       |       1        |        1.866   |        0.000   |        6.271   |  Neg. p = 0.0123  |
|       63       |       1        |        6.710   |        0.336   |        8.589   |  Neg. p = 0.0034  |
|       66       |       1        |        2.761   |        0.000   |        6.388   |  Neg. p = 0.0115  |
|       72       |       1        |        3.216   |        0.470   |        3.568   |  Neg. p = 0.0589  |
|       77       |       1        |        5.853   |        0.000   |        7.096   |  Neg. p = 0.0077  |
|       85       |       1        |        1.386   |        0.000   |        2.993   |  Neg. p = 0.0836  |
|       89       |       1        |        3.642   |        0.322   |        3.360   |  Neg. p = 0.0668  |
|       92       |       1        |        4.357   |        0.000   |        6.354   |  Neg. p = 0.0117  |
|       94       |       1        |       13.547   |        0.573   |        2.764   |  Neg. p = 0.0964  |
|       95       |       1        |        0.000   |        1.437   |        2.766   |  Pos. p = 0.0963  |

**Found _2_ sites under pervasive positive diversifying and _20_ sites under negative selection at p <= 0.1**

###### HyPhy - MEME

```
hyphy
```
- Options
    - 1 (selection)
    - 1 (MEME)
    - ickFinderCDS.correctedCDS_Ctenids_FullPeptides_ToxinSorterout_noCtenuscorniger.PAL2NAL.aln.phy (alignment)
    - ickFinderCDS.correctedCDS_Ctenids_FullPeptides_ToxinSorterout_noCtenuscorniger.PAL2NAL.aln.phy.treefile (treefile)

|   Codon    | Partition  |   alpha    |non-syn rate (beta) distribution, rates : weights|    LRT     |Episodic selection detected?| # branches |         List of most common codon substitutions at this site          |
|:----------:|:----------:|:----------:|:-----------------------------------------------:|:----------:|:--------------------------:|:----------:|:---------------------------------------------------------------------:|
|     41     |     1      |    0.000   |             0.00/24.59 : 0.00/1.00              |    3.274   |      Yes, p =  0.0924      |     1      |                  [1]aaT>aaA,AaT>GaG,gAt>gCt,GCt>AAt                   |

**Found _1_ sites under episodic diversifying positive selection at p <= 0.1**

### PnTx3-6 and 17 homologs protein modelling

#### AlphaFold 3D predictions

- Used AlphaFold 3.0 on Google Deepmind 
    - https://alphafoldserver.com
    - All sequences modelled on 03 July 2024
        - Used default parameters
        - Results placed in "4alphafold" directory

#### Protein-protein structural alignment

- Used US-align web server
    - https://zhanggroup.org/US-align/
    - All sequences aligned on XX July 2024
        - Used default parameters
        - Model results placed in "4usalign" directory

| Sample Name                                                                                       | Aligned Length | RMSD | TM-score |
| ------------------------------------------------------------------------------------------------- | -------------- | ---- | -------- |
| ctenus_captiosis_female_tjc00311_clean_trinity_trinity_trinity_dn2670_c0_g2_i1_model_0            | 50             | 1.84 | 0.66941  |
| phoneutria_nigriventer_female_h_ibsp_221739_clean_trinity_trinity_trinity_dn417_c0_g1_i2_model_0  | 52             | 2.5  | 0.63385  |
| phoneutria_nigriventer_sra_clean_trinity_trinity_trinity_dn26_c0_g1_i1_model_0                    | 50             | 2.36 | 0.61999  |
| ctenus_corniger_sra_clean_trinity_trinity_trinity_dn94115_c0_g1_i1_model_0                        | 51             | 2.46 | 0.58333  |
| ctenus_exlineae_female_tjc00244_clean_trinity_trinity_trinity_dn518_c0_g2_i2_model_0              | 51             | 2.47 | 0.60423  |
| leptoctenus_byrrhus_female_tjc00136_clean_trinity_trinity_trinity_dn2803_c0_g1_i1_model_0         | 52             | 2.77 | 0.58482  |
| ctenus_hibernalis_female_4926_bcb_0001_clean_trinity_trinity_trinity_dn178_c1_g2_i1_model_0       | 50             | 2.41 | 0.61414  |
| ctenus_hibernalis_female_tjc00092_clean_trinity_trinity_trinity_dn186_c0_g2_i1_model_0            | 54             | 2.6  | 0.6438   |
| ctenus_exlineae_female_tjc00245_clean_trinity_trinity_trinity_dn5225_c0_g1_i1_model_0             | 53             | 2.69 | 0.61664  |
| isoctenus_sp_female_f_ibsp_221743_clean_trinity_trinity_trinity_dn2656_c0_g1_i1_model_0           | 52             | 2.62 | 0.60889  |
| ctenus_corniger_sra_clean_trinity_trinity_trinity_dn26760_c0_g1_i1_model_0                        | 46             | 2.01 | 0.60282  |
| ctenus_corniger_sra_clean_trinity_trinity_trinity_dn27263_c0_g1_i1_model_0                        | 47             | 2.48 | 0.58569  |
| ctenus_corniger_sra_clean_trinity_trinity_trinity_dn27263_c0_g2_i1_model_0                        | 46             | 2.27 | 0.6092   |
| ctenus_corniger_sra_clean_trinity_trinity_trinity_dn8463_c0_g1_i1_model_0                         | 49             | 2.65 | 0.60351  |
| isoctenus_aff_ordinario_female_g_ibsp_221742_clean_trinity_trinity_trinity_dn167_c0_g1_i5_model_0 | 52             | 2.78 | 0.58605  |
| phoneutria_nigriventer_female_i_ibsp_221740_clean_trinity_trinity_trinity_dn49_c0_g1_i4_model_0   | 54             | 2.51 | 0.64526  |
| phoneutria_nigriventer_female_j_ibsp_221741_clean_trinity_trinity_trinity_dn218_c0_g1_i10_model_0 | 51             | 2.9  | 0.5718   |

##### Test for phylogenetic signal in RMSD scores

- Used R v4.3.1
    - readr_2.1.5
    - phylosignal_1.3.1
    - phylobase_0.8.12 
    - adephylo_1.1-16
    - ade4_1.7-22
    - ape_5.8          
    - stringr_1.5.1
    - tibble_3.2.1
    - dplyr_1.1.4
- Methods and results in PnTx3-6_PhylogeneticSignalRTests.Rmd
    - Rendered in PnTx3-6_PhylogeneticSignalRTests.nb.html

```
pandoc PnTx3-6_PhylogeneticSignalRTests.nb.html --pdf-engine wkhtmltopdf -o PnTx3-6_PhylogeneticSignalRTests.nb.pdf
```
