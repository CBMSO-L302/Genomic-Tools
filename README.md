# Genomic tools

The genomic tools we present here were created with the purpose of extracting information from eukaryotic DNA files in FASTA format.

## coord_prot_to_gff.py

Recommendation:

Before the usage of these scripts, save the path in a variable, for example, `$SCRIPT`:

`SCRIPT=/Path_to_script`

Where "Path_to_script" should be substituted by the folder path where you have the script.

_Usage:_

In terminal, we launch this python script using the following arguments:

`python $SCRIPT/coord_prot_to_gff.py [1] [2] [3]`

`[1]`: **FASTA file path** with the **DNA** sequence(s)

`[2]`: has 2 options:

- **file path** in .txt with the id names of the selected fasta sequences from the FASTA file in `[1]`

    NOTE: don't add '>' to the names

- **'all'**
If you want to extract all the proteins from all the sequences in the FASTA file


`[3]`: has 2 options:

- **'joined'** if you want all the files to be extracted in the same gff file (the file will be automatically named 'joined.gff.txt')
        
- **'separated'** if you want to extract all the files separatedly. 
The files will be automatically named with the name of the fasta sequence and the .gff.txt extension

_Example of usage 1:_

`python $SCRIPT/coord_prot_to_gff.py LdHU3.fasta 'all' 'joined'`

This will read ALL the sequences from the LdHU3.fasta file and will have the output joined.gff.txt

_Example of usage 2:_

`python $SCRIPT/coord_prot_to_gff.py LdHU3.fasta wanted.txt 'separated'`

This will read the target sequences specified in wanted.txt and output them in separated files with the format 'sequencename'.gff.txt

Example of WANTED file (_wanted.txt_): 
If I needed the chromosomes 01 07 and 20 from the LdHU3 genome, the _wanted.txt_ file would be:

```
LdHU3.01
LdHU3.07
LdHU3.20
```
_Example of output:_

_LdHU3.36.gff.txt_
```
LdHU3.36	Predicted Protein	CBMSO	119	370	.	+	.	LdHU3.36:119..370
LdHU3.36	Predicted Protein	CBMSO	153	356	.	-	.	LdHU3.36:153..356..r
LdHU3.36	Predicted Protein	CBMSO	166	426	.	-	.	LdHU3.36:166..426..r
LdHU3.36	Predicted Protein	CBMSO	288	485	.	+	.	LdHU3.36:288..485
LdHU3.36	Predicted Protein	CBMSO	355	501	.	+	.	LdHU3.36:355..501
LdHU3.36	Predicted Protein	CBMSO	407	499	.	-	.	LdHU3.36:407..499..r
```

The colums of the gff file can be modified within the script.

The output file will be created inside the current working directory.

## DNA_to_proteins_fasta.py

The usage is exactly as **coord_prot_to_gff.py**, with the form:

`python $SCRIPT/DNA_to_proteins_fasta.py [1] [2] [3]`

In this case, the extension output will be **.prot.fasta**

_Example of usage 1:_

`python $SCRIPT/DNA_to_proteins_fasta.py LdHU3.fasta wanted.txt 'joined'`

It will look for the 'wanted' sequences in the LdHU3.fasta file and will have the output joined.prot.fasta

_Example of usage 2:_

`python $SCRIPT/DNA_to_proteins_fasta.py LdHU3.fasta 'all' 'separated'`

It will extract all the proteins from all the reading frames from each sequence in the fasta file and extract them in separated files, one per sequence given in the LdHU3.fasta file.

_Example of output:_

_LdHU3.02.prot.fasta_
```
>LdHU3.02:577..666
MGLESARVVVVHVLRVCVCTLVGRVQDAH
>LdHU3.02:784..942
METPVKKKKVGHCGYLSQLVLLATRKRVDSHSRQACVREAGEKRLSKVLRQR
>LdHU3.02:1756..2268
MAIDNLVEFVPPHSIQAFLCVVDAQRRWDNPFVRNGRRQVARFHPPTAHMAQVAQVLWHLQVRVVVALHDVRHVWEGSVRVVQPYHDLLSHAYFRRETPRYRIAHFFFVVGLRSVKKSSHSELRWCSASCHNAADTWHTYCAAAVKAGVTGSDNWAERAARGEVPSTRLA
>LdHU3.02:2479..3132
MRCVLRITLRSSEHVELQQRASKVVLGARNLSVRQPRRPLRVPQRVALRLVRGRYAHHDQVPRSLTCAPRLRRVDEAALLRVRLDHPLPPRQFSCYVHVVRHARHGKQHGPHEAQRRRCCREHRKPFQPLVSRPVSGSPRIRSNRRSVRHIVEAHRAVLVDHAQRPRRQLITVQSHQRVIFRGSKTALCGRQPLSSSLVVHENRKRRWKLRFTGRLL
>LdHU3.02:3225..3341:r
MHGRMGYHRSFPLGYADILDDPRHNACAYKRNAVFWLR
>LdHU3.02:2781..2900:r
MRGDPLTGRETRGWKGFLCSRQHLRRWASWGPCCLPWRA
>LdHU3.02:1563..1634:r
MRACSCSTRTFSWDGNCKVFLLG
>LdHU3.02:1068..1142:r
MCHRVWKVVQHVSFLSPSTLGRQA
```

The output file will also be created inside the current working directory.




