# Synthedia

A wide range of software is available to identify and quantify peptides in bottom-up proteomics data acquired by LC-MS/MS and, frequently, analysis of the same input data file with two different packages yields different results. These differences likely originate from differences in algorithms for preprocessing data, matching MSn spectra to peptides, quantifying peak areas, matching assignments between samples as well as differences in the myriad parameters that must be set by the user to initiate an analysis. Thorough analysis of these variables is, however, complicated since the ‘true’ sample composition of real data files is almost never completely known.

Synthedia creates synthetic LC-MS/MS data that mimics real data but with a composition that is exactly known. Currently, synthedia support the creation of Data-Independent Acquisition (DIA) style data wherein fixed, large m/z windows are sequentially isolated for fragmentation. We have focused on creating DIA data to date since the complexity of the analysis preformed by processing tools is substantial and the impact of different acquisition methodologies on the eventual outcome is somewhat more difficult to predict.

Synthedia can be configured to produce synthetic DIA data that mimics a wide range of acquisition strategies. For example, data can be simulated that models the same set of peptides eluted over progressively shorter gradients which would allow an assessment of how processing software copes with increasing complexity of data. Data can be modelled with different mass spectral resolutions, chromatographic peak tailing, DIA window strategies among others. Synthedia can also help to assess quantitation by creating multiple sets of data files which simulate replicates in a multi-group comparison. These modules include abilities to set within and between-group variability as well as simulate missing data in both a random and group-wise manner. 

## Getting started

Synthedia can be used by downloading and installing the package or through our web server which can be accessed [here](http://45.113.234.202:8502/). Please note that producing synthetic data is a time-consuming task and capacity on our server is limited. In most cases, it will be faster to install synthedia on your own computer as this can take advantage of multiple processing cores.

## Installation

Clone the repo:
```
git clone https://github.com/mgleeming/synthedia.git
cd synthedia
```
Create a virtual environment (optional):
```
virtualenv venv
source venv/bin/activate
```

Note that the [pyOpenMS](https://pyopenms.readthedocs.io/en/latest/index.html) dependency of synthedia requires python 3.7, 3.8 or 3.9. If you have multiple python versions installed, you can construct virtual environments with different versions by:
```
virtualenv venv --python=/path/to/python3.X
```
On linux systems, the path is usually ```/usr/bin/pythonX.X```

Install synthedia:
```
pip install .
```

## Input data types

To create synthetic data, synthedia requires information about peptide fragmentation patterns and some metric as to relative retention time. These can be supplied by either a:
  - [Prosit](https://www.proteomicsdb.org/prosit/)-predicted spectral library
  - MaxQuant 'txt' directory

#### Prosit libraries


#### MaxQuant 'txt' directories
## Usage

## Viewing mzML files

