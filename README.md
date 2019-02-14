# SJARACNe
SJARACNe is a scalable solution of ARACNe that dramatically improves the computational 
performance, especially on the memory usage to allow even researchers with modest 
computational power to generate networks from thousands of samples.


## Download
```git clone https://github.com/jyyulab/SJARACNe  # Clone the repo```


## Prerequisites
* [Python 3.6.1](https://www.python.org/downloads/)
	* [numpy==1.14.2](https://www.scipy.org/scipylib/download.html)
	* [argparse==1.1](https://docs.python.org/3/library/argparse.html)
	* [python-igraph==0.7.1](https://igraph.org/python/)
	* [scipy==1.0.1](https://www.scipy.org/install.html)
	* [XlsxWriter==1.0.2](https://xlsxwriter.readthedocs.io/)
	* [pandas==0.22.0](https://pandas.pydata.org/)


## Installation
### Using conda to create a virtual environment (recommended)
The recommended method of setting up the required Python environment and dependencies is to use the
[conda](https://conda.io/en/latest/) dependency manager:

```bash
$ conda create -n py36 python=3.6.1
$ conda activate py36
$ conda install --file requirements.txt
```

### Using pip
First install [Python 3.6.1](https://www.python.org/downloads/) and then use the following command to install package requirements.

```$ pip install -f requirements.txt```


### Install from source
The linux and OSX pre-built distribution are provided in `SJARACNe/bin` and the program will use the corresponding 
distribution with respect the operating system. Alternatively, you may use the Makefile to compile the code and build your own distribution.
You can install SJARACNe directly from the source using `setup.py`:

```bash
$ git clone https://github.com/jyyulab/SJARACNe
$ cd SJARACNe
$ python setup.py install
```


## Usage
```$ sjaracne [project_name] [expression_matrix] [hub_genes] [output_directory]```

### Options
* ```--bootstrap, default=100, Number of bootstrap networks.```
* ```--c_threshold, default=1e-5, P-value threshold in building consensus network.```
* ```--p_threshold, default=1e-7, P-value threshold in building bootstrap netwroks.```
* ```--depth, default=40, help=Maximum partitioning depth.```
* ```--run, default=False, help=Whether run the pipeline or just generate and stop.```
* ```--host, default=LOCAL, help=Whether to run on clusters or localhost. [LOCAL | LSF].```

### Notes:
1. Setting the host option to LSF will change the run option to False.
2. Absolute / relative filepaths _without_ any environmental variables (e.g. `$HOME`) must be used.


## Examples
### EASY RUN
```$ sjaracne [project_name] [expression_matrix] [hub_genes] [output_directory]```

The above command will create 4 directories under the provided out_directory parameter as follows:

* ```[out_directory]/SJARACNE_[project_name]/SJARACNE_log```
* ```[out_directory]/SJARACNE_[project_name]/SJARACNE_out.final```
* ```[out_directory]/SJARACNE_[project_name]/SJARACNE_scripts```

There will be shell script files corresponding to the provided input files in the scripts 
directory in the following order:

* ```00_cleanup_[project_name].sh</code>```
* ```00_pipeline_[project_name].sh</code>```
* ```01_prepare_[project_name].sh```
* ```02_bootstrap_[project_name].sh```
* ```03_getconsensusnetwork_[project_name].sh```
* ```04_getenhancedconsensusnetwork_[project_name].sh```

The command will run scripts 02-04 automatically and generate the final results.


### Example of Running the Scripts on IBM LSF Cluster
```$ sjaracne TF ./test_data/inputs/BRCA100.exp ./test_data/inputs/tf.txt ./test_data/outputs/ --host LSF```

### Example of Running Signaling Network on a Single Machine (Linux/OSX)
```$ sjaracne SIG ./test_data/inputs/BRCA100.exp ./test_data/inputs/sig.txt ./test_data/outputs/```

### Example of Running Transcription Factor Network on a Single Machine (Linux/OSX)
```$ sjaracne TF ./test_data/inputs/BRCA100.exp ./test_data/inputs/tf.txt ./test_data/outputs/```

### Expected Output
Expected output for the example data with 100 bootstraps is available under 
```test_data/outputs/SJARACNE_TF/SJARACNE_out.final``` directory.


## Reference
Alireza Khatamian, Evan O. Paull, Andrea Califano* & Jiyang Yu*. SJARACNe: a scalable software tool for gene network reverse engineering from big data. Bioinformatics (2018). * Corresponding authors.
