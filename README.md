# SJARACNe
[![Build Status](https://travis-ci.com/jyyulab/SJARACNe.svg?branch=master)](https://travis-ci.com/jyyulab/SJARACNe)

SJARACNe is a scalable solution of ARACNe that dramatically improves the computational 
performance, especially on the memory usage to allow even researchers with modest 
computational power to generate networks from thousands of samples. The algorithm uses adaptive 
partitioning mutual information to calculate the correlation between all pairs of genes to 
reconstruct the regulatory network.


## Download
```git clone https://github.com/jyyulab/SJARACNe  # Clone the repo```


## Prerequisites
* [Python 3.6.1](https://www.python.org/downloads/)
	* [numpy==1.14.2](https://www.scipy.org/scipylib/download.html)
	* [scipy==1.0.1](https://www.scipy.org/install.html)
	* [pandas==0.22.0](https://pandas.pydata.org/)
	* [cwltool==1.0.20190618201008](https://github.com/common-workflow-language/cwltool/releases)
* [cwlexec==0.2.2](https://github.com/IBMSpectrumComputing/cwlexec/releases) (required for running on IBM LSF)


## Installation
### Using conda to create a virtual environment (recommended)
The recommended method of setting up the required Python environment and dependencies is to use the
[conda](https://conda.io/en/latest/) dependency manager:

```bash
$ conda create -n py36 python=3.6.1
$ source activate py36
$ conda install --file requirements.txt
```


### Using pip
First install [Python 3.6.1](https://www.python.org/downloads/) and then use the following command to install package requirements.

```$ pip install -f requirements.txt```


### Install from source
```bash
$ git clone https://github.com/jyyulab/SJARACNe
$ cd SJARACNe
$ python setup.py install
```


## Usage
```$ sjaracne 
usage: sjaracne [-h] {local,lsf} ...

SJARACNe is a scalable tool for gene network reverse engineering.

optional arguments:
  -h, --help   show this help message and exit

Subcommands:
  {local,lsf}  platforms
    local      run cwltool in a local workstation
    lsf        run cwlexec as in a IBM LSf cluster
```
```sjaracne``` workflow is implemented with [CWL](https://www.commonwl.org/). It supports multiple
 computing platforms. We have tested it locally using [cwltool](https://github.com/common-workflow-language/cwltool) 
 and on an IBM LSF cluster using [cwlexec](https://github.com/IBMSpectrumComputing/cwlexec). 
 For the convenience, a python wrapper is developed for you to choose computing platform using ```subcommand```.

The local mode (```sjaracne local```) runs in parallel by default using cwltool's ```--parallel``` option. 
To run it in serial, use ```--serial``` option.

To use LSF mode, editing the LSF-specific configuration file ```SJARACNe/config/config_cwlexec.json``` to change the default queue and
adjust memory reservation for each step is necessary. Consider increasing memory reservation for bootstrap step and consensus step if
the dimension of your expression matrix file is large.


### Inputs
The main input for SJARACNe is a tab-separated genes/protein by cells/samples expression matrix
with the first two columns being ID and symbol. The second required input file is the list of
significant genes/proteins IDs to be considered as hubs in the reconstructed network. An output directory is required
for storing output files. Additional parameters (e.g., LSF queue) for running on different platforms are required. 
Those are available in the helping information of the corresponding subcommands, e.g., ```sjaracne lsf -h```.


### Outputs
The main output of SJARACNe is a network file, which is a tab delimited text file with the following columns: source,
target, mutual information, Pearson and Spearman correlations coefficients, regression line slope and p-value. SJARACNe
also outputs two meta information files: parameter_info_.txt and bootstrap_info_.txt, which stores SJARACNe 
input parameters and bootstrap parameters respectively.


## Examples to create a transcription factor network
### Running on an IBM LSF cluster
```sjaracne local -e ./test_data/inputs/BRCA100.exp -g ./test_data/inputs/tf.txt -n 2 -o ./test_data/outputs/cwl/cwlexec```


### Running on a single machine (Linux/OSX) 
```sjaracne lsf -j ./SJARACNe/config/config_cwlexec.json -e ./test_data/inputs/BRCA100.exp -g ./test_data/inputs/tf.txt -n 2 -o ./test_data/outputs/cwl/cwltool```


## Reference
Alireza Khatamian, Evan O. Paull, Andrea Califano* & Jiyang Yu*. SJARACNe: a scalable 
software tool for gene network reverse engineering from big data. Bioinformatics (2018). *Corresponding authors.

