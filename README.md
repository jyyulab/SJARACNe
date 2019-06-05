# SJARACNe
SJARACNe is a scalable solution of ARACNe that dramatically improves the computational 
performance, especially on the memory usage to allow even researchers with modest 
computational power to generate networks from thousands of samples.


## Download
```git clone https://github.com/jyyulab/SJARACNe  # Clone the repo```


## Prerequisites
* [Python 3.6.1](https://www.python.org/downloads/)
	* [numpy==1.14.2](https://www.scipy.org/scipylib/download.html)
	* [python-igraph==0.7.1.post6](https://igraph.org/python/)
	* [scipy==1.0.1](https://www.scipy.org/install.html)
	* [XlsxWriter==1.0.2](https://xlsxwriter.readthedocs.io/)
	* [pandas==0.22.0](https://pandas.pydata.org/)


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
```sjaracne``` workflow is implemented with [CWL](https://www.commonwl.org/). It should run in multiple
 computing platforms. We have tested it locally using [cwltool](https://github.com/common-workflow-language/cwltool) 
 and on an IBM LSF cluster using [cwlexec](https://github.com/IBMSpectrumComputing/cwlexec). 
 For the convenience, a python wrapper is developed for you to choose computing platform using ```subcommand```.
 Please refer to link? for the more information, e.g., inputs and outputs.


## Examples to create a transcription factor network
### Running on an IBM LSF cluster
```sjaracne local -e ./test_data/inputs/BRCA100.exp -g ./test_data/inputs/tf.txt -n 2 -o ./test_data/outputs/cwl/cwltool/SJARACNE_out.final```

### Running on a single machine (Linux/OSX) 
```sjaracne lsf -e ./test_data/inputs/BRCA100.exp -g ./test_data/inputs/tf.txt -n 2 -o ./test_data/outputs/cwl/cwltool/SJARACNE_out.final -q [queue_name]```

## Reference
Alireza Khatamian, Evan O. Paull, Andrea Califano* & Jiyang Yu*. SJARACNe: a scalable 
software tool for gene network reverse engineering from big data. Bioinformatics (2018). *Corresponding authors.
