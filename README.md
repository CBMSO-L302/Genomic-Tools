# Genomic tools

The genomic tools we present here were created with the purpose of extracting information from eukaryotic DNA files in FASTA format.

## coord_prot_to_gff.py

Usage:

In terminal, we launch this python script using the following arguments:

`python coord_prot_to_gff.py [1] [2] [3]`

`[1]`: **FASTA file path**

`[2]`: has 2 options:

- **WANTED file path** in .txt 
with the id names of the selected fasta sequences from the FASTA file in `[1]`

    NOTE: don't add '>' to the names

- **'all'**
If you want to extract all the proteins from all the sequences in the FASTA file


`[3]`: has 2 options:

- **'joined'** if you want all the files to be extracted in the same gff file (the file will be automatically named 'joined.gff.txt')
        
- **'separated'** if you want to extract all the files separatedly. 
The files will be automatically named with the name of the fasta sequence and the .gff.txt extension

_Example of usage 1:_

`python coord_prot_to_gff.py LdHU3.fasta 'all' 'joined'`

This will read ALL the sequences from the LdHU3.fasta file and will have the output joined.gff.txt

_Example of usage 2:_

`python coord_prot_to_gff.py LdHU3.fasta wanted.txt 'separated'`

This will read the target sequences specified in wanted.txt and output them in separated files with the format 'sequencename'.gff.txt

Example of WANTED file (wanted.txt): 
In case I wanted the chromosomes 01, 07 and 20 from the LdHU3 genome, the wanted.txt file would be:

```
LdHU3.01
LdHU3.07
LdHU3.20
```

The colums of the gff file can be modified in the script.
