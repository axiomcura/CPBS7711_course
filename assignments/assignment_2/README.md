# Table of Contents
[Net_sim.py ](#netsimpy)\
[Dependencies ](#dependencies)\
 [Using pip](#using-pip)\
 [Using conda](#using-conda)\
[Work Flow ](#work-flow)\
 [Arguments ](#arguments)\
 [Required Inputs ](#required-inputs)\
  [Formats of the files:](#formats-of-the-files)\
  [gmt format ](#gmt-format)\
  [String Database ](#string-database)\
 [Commands ](#command)\
 [Results ](#results)

---
# Net_sim.py 

net_sim.py is a light weight script that generates random subnetworks
from a given `.gmt` file that contains a tab delimited file format of
gene sets. It samples randomly generated subnetwork, where it contains
one gene per locus, and compares it to the larger gene network contained
in the `STRING.txt` database. This allows understanding how significant
these subnetworks are when comparing two different states.

# Dependencies 

This program requires two external package

``` {.}
pandas v1.2.5 
numpy  v1.20
```

Below is how to install these packages using either `pip` or `condo` as
your python package manager .

## Using pip

``` {.}
pip install -Iv pandas==1.2.5
pip install -Iv numpy==1.20
```

## Using conda

``` {.}
conda install pandas=1.2.5
conda install numpy=1.20
```

If your installations are more specialized, here are the documentations
for both package managers:

`conda install`:
<https://docs.conda.io/projects/conda/en/latest/commands/install.html>

`pip install` : <https://pip.pypa.io/en/stable/cli/pip_install/>

# Work Flow 

## Arguments 

To have open the help documentation is accessed by typing
`python net_sim.py -h`

``` {.markdown}
usage: net_sim.py [-h] [-i INPUT] [-db FILE] [-s PARAM] [-x PARAM] [--n_test PARAM]

Simple program that creates random subnetworks and is compared to STRING database network

optional arguments:
  -h, --help            show this help message and exit
  -db FILE, --database FILE
                        path to database. Default path is `./Data/STRING.txt
  -s PARAM, --shuffles PARAM
                        Number of shuffles done in permutations test. Default is 500
  -x PARAM, --iterations PARAM
                        Number of sampling iterations. Default is 500
  --n_test PARAM        Number of permutation test conducted

Required Arguments:
  -i INPUT, --input INPUT
                        input file
```

## Required Inputs 

`net_sim.py` requires the user to have a database and an input file In
`.gmt` format. These two files will be responsible for generating the
sub networks.

### Formats of the files:

Bellow are the formats of each file that is being used in the program:

### gmt format 

``` {.}
Fanconi anemia locus 0  Locus for PALB2 NUPR1 CTB-134H23.2
SLC5A11 KIAA0556  CD19  SH2B1
Fanconi anemia locus 1  Locus for FANCF CSTF3 FBXO3 SLC17A6
CCDC73  CAPRIN1 RCN1  BDNF  METT
Fanconi anemia locus 2  Locus for RAD51C, BRIP1 CD79B CACNG1
TANC2 SMG8  RP11-15E18.4  TEX2
Fanconi anemia locus 3  Locus for FANCC TMOD1 MSANTD3-TMEFF1
TEX10 HIATL2  LPPR1 MRPL50
Fanconi anemia locus 4  Locus for FANCA DBNDD1  AC133919.6
RP11-356C4.2  MC1R  SPIRE2  C16orf
Fanconi anemia locus 5  Locus for UBE2T PPFIA4  SNRPE LGR6
ZC3H11A TMEM183A  PPP1R12B  SOX13
Fanconi anemia locus 6  Locus for FANCD2  VGLL4 TAMM41
ZFYVE20 EAF1  IQSEC1  SLC6A6  VHL FANC
Fanconi anemia locus 7  Locus for FANCE KCNK5 MTCH1 PI16
BRPF3 TMEM217 RPL10A  GLP1R ARMC1
Fanconi anemia locus 8  Locus for BRCA2 PROSER1 MRPS31  UFM1
AKAP11  N4BP2L2 TNFSF11 FOXO1
Fanconi anemia locus 9  Locus for ERCC4 GPR139  C16orf88
RPS15A  BFAR  CTD-2349B8.1  NDE1
Fanconi anemia locus 10 Locus for FANCI POLG  AC016251.1
FANCI GABARAPL3 FES C15orf38
Fanconi anemia locus 11 Locus for SLX4  FAM86A  USP7  ROGDI
ALG1  CDIP1 UBN1  RP11-297M9.1
```

Each row is a locus and the associated genes downstream from the locus
name

### String Database 

``` {.}
ARF5  DVL2  0.166000
ARF5  DYRK4 0.166000
ARF5  PPP5C 0.254968
ARF5  MAP4K5  0.157276
ARF5  RALBP1  0.156000
ARF5  PKP2  0.160210
ARF5  ACAP1 0.328000
ARF5  MAP2K5  0.242000
ARF5  MYO15A  0.272395
ARF5  MAPK13  0.190000
ARF5  STX1B 0.263160
ARF5  MAPK12  0.19000
...
```

Each row represents that interactions of two different genes and their
score.

## Commands 

To execute this `net_sim.py`, enter:

``` {.}
python net_sim.py -i ./Data/Input.gmt.txt
```

If you have your own dataset, you can replace the value in `-i` to the
locations where you `.gmt` input file is.

``` {.}
python net_sim.py -i /path/to/your/input/file.gmt
```

One can add the optional arguments as well. To get more information
about these optional arguments, use the help flag to display their
function and default values

This example shows that I want to increase the number of test iterations
from 10 (default) to 15 and want to increase the number of shuffling
from 500 (default) to 1000 in the permutation test.

``` {.}
python net_sim.py -i ./Data/Input.gmt.txt -s 1000 --n_test 15
```

Keep in mind that these optional arguments are not required to be
explicitly when executing the program. Internally, these arguments have
a default value that the algorithm uses.

The user only explicitly type the optional arguments if they want to
change their default values, hence the example above.

## Results 

A `.csv` should be generated into the current directory that contains
the number of permutation test conducted and its p_value.

Opening `.csv` file in terminal by using `cat results.csv` in
terminal/console. These file can also be opened in excel as well.

``` {.}
test,p_val
0,0.972
1,0.888
2,0.956
3,0.866
4,0.932
5,0.94
6,0.912
7,0.932
8,0.998
9,0.994
```
