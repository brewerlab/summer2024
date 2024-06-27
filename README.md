# Summer 2024 UG research project
# Ctenid Tx toxins evolution and protein sturcture analyses

### Michael S. Brewer
### Jordan Baccus
### Amir El

## Sequence collection and curration

### Download ickFinderCDS.correctedCDS.csv from tijeco/killerknots GitHub

```
mkdir ctenid_sequences
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

\>sp|P17727.2|TXL1_PHONI RecName: Full=Mu-ctenitoxin-Pn1a; Short=Mu-CNTX-Pn1a; AltName: Full=Toxin Tx1; Short=PNTx1; Short=PhTx1; Flags: Precursor
MKLLGIFLVASFAFVLSFGEEMIEGENPLEDQRAELTSCFPVGHECDGDASNCNCCGDDVYCGCGWGRWN
CKCKVADQSYAYGICKDKVNCPNRHLWPAKVCKKPCRRNCGG
 
\>sp|P29425.2|TX36A_PHONI RecName: Full=Delta-ctenitoxin-Pn2a; Short=Delta-CNTX-Pn2a; AltName: Full=Delta-CNTX-Pn1c; AltName: Full=Neurotoxin Tx2-6; Short=PnTx2-6; Flags: Precursor
MKVAILFLSILVLAVASESIEESRDDFAVEELGRATCAGQDQPCKETCDCCGERGECVCGGPCICRQGYF
WIAWYKLANCKK
 
\>sp|P81792.2|TX90B_PHONI RecName: Full=Omega-ctenitoxin-Pn4a; Short=Omega-CNTX-Pn4a; AltName: Full=CTK 01512-2; AltName: Full=Neurotoxin Tx3-6; Short=PnTx3-6; AltName: Full=Ph-alpha-1-beta toxin; Flags: Precursor
MKCAVLFLSVIALVHIFVVEAEEEPDSDALVPQERACIPRGEICTDDCECCGCDNQCYCPPGSSLGIFKC
SCAHANKYFCNRKKEKCKKA
 
\>sp|P59367.2|TX35C_PHONI RecName: Full=GAMMA-ctenitoxin-Pn1a; Short=GAMMA-CNTX-Pn1a; AltName: Full=Insecticidal neurotoxin Tx4(5-5); Short=PnTx4(5-5); Short=PnTx4-5-5; AltName: Full=Toxin Pn4A; Flags: Precursor
MKVAIVFLSLLVLAFASESIEENREEFPVEESARCADINGACKSDCDCCGDSVTCDCYWSDSCKCRESNF
KIGMAIRKKFC

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
fasta2phylip.py -i ickFinderCDS.correctedCDS_Ctenids_FullPeptides_withPnTxFromGenbank_MAFFT.aln.fasta -o ickFinderCDS.correctedCDS_Ctenids_FullPeptides_withPnTxFromGenbank_MAFFT.aln.phy -r
```

### Phylogeny inference using IQTREE2

```
iqtree -s ickFinderCDS.correctedCDS_Ctenids_FullPeptides_withPnTxFromGenbank_MAFFT.aln.phy -T AUTO -B 1000 -m MFP
```


