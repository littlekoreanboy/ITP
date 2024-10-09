### Ion Transport Peptide

## Important Reading
```bash
Doi : https://doi.org/10.7554/eLife.97043.1
```

## Dr. Zandawala meeting
```bash
* Neuropeptide and Receptors are known as "Lock" and "Key" mechanism where a specific Key can only unlock a specific Lock. Therefore in terms of Neuropeptide and Receptors, a specific Neuropeptide can only modulate the activity of a Receptor that is specific for that Neuropeptide
* Protostomes and Deuterostomes can be rounghly described as Invertebrates vs Vertebrates. Take for example : Flies are invertebrates, because they lack a spinal column, and Humans are vertebrates, because we have a spinal column
* Protostomes and Deuterostomes are the 2 major clades (groups) of bilaterians: bilaterally symmetrical animals
* In co-evolution, neuropeptide ITP and its receptor ITPR SHOULD have evolved together due to their interdependence, where changes in one likely drive adaptations in the other
* It is known that Invertebrates express ITP. So to identify the receptor for ITP, we will find receptors that are only found in invertebrates. ITPR gene has not been identified. Therefore finding ITPR sequence is important
* The conserved amino-acid cysteines are found in all ITP. Specifically, there are 6 cysteine residues
```
## What is ITP
```bash
* ITP stands for Ion Transport Peptide
* ITP has the following similar structure : (Signlaing Peptide) - X - KR(Propeptide 1)KR - RR(Propeptide 2)RR - X, where X represents any amino acid sequence 
* KR, RR, RK, RRR are direct splice sites that results in one or more Propeptide isoform. Each precursor can give rise to one or more mature neuropeptides
* Additonally, different species can have more than one mature neuropeptides
* Interestingly, A long Propeptide (After it have been spliced) is associated as a Hormone
``` 

## Project Task - (Timeline : 10/3/2024 - 10/7/2024)
```bash
# Our Goal is to BLAST search 20 ITP protein sequences for each arthropod order. We will use Multiple Sequence Alignment to access their evolutionary relationship, common patterns, and conservation

$DATA
* Orders of Arthropods : https://en.wikipedia.org/wiki/List_of_arthropod_orders
* NCBI Protein sequences
* Data links:
1) https://flybase.org/reports/FBgn0035023.htm
2) https://www.uniprot.org/uniprotkb/E1JGW3/entry
3) https://www.uniprot.org/uniprotkb/Q0E8W5/entry
4) https://www.uniprot.org/uniprotkb/Q9W151/entry

$INPUT
* Orders of Arthropods

$OUTPUT
* Table summarizing presence or absence of ITP
* Aligned Sequence + Accession ID
* Weblogo 
```

## Additional Program Readings
```bash
https://biopython.org/docs/latest/api/Bio.Blast.NCBIXML.html
https://biopython.org/docs/latest/api/Bio.Blast.NCBIWWW.html#Bio.Blast.NCBIWWW.qblast
https://weblogo.threeplusone.com/manual.html
https://pypi.org/project/pyMSAviz/
```
## Before running
* Please run this following conda if you have conda installed:
```bash
conda create -y -n ITP -c bioconda -c conda-forge --file installation.txt && conda activate ITP
```
* This conda command will install all the necessary libraries. 

## Basic Utalization
* Please Run the following command:
```bash
$python run.py
```
* `run.py` is the main source. To change vairables, such as chaning the BLAST program or protein/nucleotide sequences, please edit the `config.yaml` file. For example, here is what it looks like:
```bash
itp :
  program : 'blastp'
  database : 'nr'
  sequence : 'MCSRNIKISVVLFLVLIPIFAALPHNHNLSKRSNFFDLECKGIFNKTMFFRLDRICEDCYQLFRETSIHRLCKQECFGSPFFNACIEALQLHEEMDKYNEWRDTLGRK'
  order : 
    Diptera : 7147
    Orthoptera : 6993
    Hemiptera : 7524
    Lepidoptera : 7088
    Blattodea : 85823
  hit_list : 20
  raw_folder: 'raw'
  outfile : 'raw/raw_itp_BLAST_results.fa'
  muscleOutput : 'raw/raw_muscled_itp_BLAST_results.fa'
  itp_conserved_regions:
    cys_1 : 57
    cys_2 : 73
    cys_3 : 76
    cys_4 : 89
    cys_5 : 93
    cys_6 : 102
  output_folder: 'output'
  sorted_muscleFile : 'output/files/sorted_muscled_itp_BLAST_results.fa'
  tableFile : 'output/files/table.csv'
  weblogo : 'output/figures/weblogo.pdf'
  pymsaviz : 'output/figures/pymsaviz.pdf'
```

## To do
```bash
### Try and Find tachykinin
### Change e-values 0.001 or less have to test this out
### Pubmeb, blast search, research paper to varify the sequence is a NO and nothing else
### Get authentc no's to test hypothesis 
```
