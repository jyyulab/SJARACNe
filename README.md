# SJARACNe

## Reference

Alireza Khatamian, Evan O. Paull, Andrea Califano* & Jiyang Yu*. SJARACNe: a scalable software tool for gene network reverse engineering from big data. Bioinformatics (in press) (2018). * Corresponding authors.

## Download

<code>git clone https://github.com/jyyulab/SJARACNe.git</code>

## Requirements

* [Python 3.6.1](https://www.python.org/downloads/)
	* numpy==1.14.2
	* argparse==1.1
	* igraph==0.1.11
	* scipy==1.0.1
	* XlsxWriter==1.0.2
	* pandas==0.22.0

## Install Requirements
First install [Python 3.6.1](https://www.python.org/downloads/) and then use the following command to install package requirements.
</br></br>
<code>pip install -f dependencies.txt</code>

## Install SJARACNe
The linux and OSX pre-built distribution are provided and the program will use the corresponding distribution with respect the operating system.</br>

Use the Makefile to compile the code and build your own distribution.</br>

## Set Environemnt

<code>$export PYTHON_PATH=$(which python3)</code></br>
<code>$export SJARACNE_PATH=$(pwd)</code></br>

## Main Command

<code>$python3 generate_pipeline.py [project_name] [expression_matrix] [hub_genes] [out_directory]</code></br>

### Options

* <code>--bootstrap, default=100, Number of bootstrap networks.</code></br>
* <code>--c_threshold, default=1e-5, P-value threshold in building consensus network.</code></br>
* <code>--p_threshold, default=1e-7, P-value threshold in building bootstrap netwroks.</code></br>
* <code>--depth, default=40, help=Maximum partitioning depth.</code></br>
* <code>--run, default=False, help=Whether run the pipeline or just generate and stop.</code></br>
* <code>--host, default=LOCAL, help=Whether to run on clusters or localhost. [LOCAL | CLUSTER].</code></br>

NOTE: Setting the host option to CLUSTER will change the run option to False.

## EASY RUN

<code>$python3 generate_pipeline.py [project_name] [expression_matrix] [hub_genes] [out_directory] --run True</code></br>

<p>The above command will create 4 directories under the provided out_directory parameter as follows:</p></br>

* <code>[out_directory]/sjaracne_[project_name]\_out\_</code></br>
* <code>[out_directory]/sjaracne_[project_name]\_log\_</code></br>
* <code>[out_directory]/sjaracne_[project_name]\_out\_.final</code></br>
* <code>[out_directory]/sjaracne_[project_name]\_scripts\_</code></br></br>

<p>There will be shell script files corresponding to the provided input files in the scripts directory in the following order:</p></br>

* <code>00_cleanup_[project_name].sh</code></br>
* <code>00_pipeline_[project_name].sh</code></br>
* <code>01_prepare_[project_name].sh</code></br>
* <code>02_bootstrap_[project_name].sh</code></br>
* <code>03_getconsensusnetwork_[project_name].sh</code></br>
* <code>04_getenhancedconsensusnetwork_[project_name].sh</code></br>

<p>The command will run scripts 02-04 automatically and generate the final results.</p></br>

<p>To generate only scripts and stop, run the command as follows:</p></br>
<code>$python3 generate_pipeline.py [project_name] [expression_matrix] [hub_genes] [out_directory]</code></br>

## Running on Single Machine (Linux/OSX)
<code>$python3 generate_pipeline.py [project_name] [expression_matrix] [hub_genes] [out_directory] --run True</code></br>

<p>OR</p>

<code>$python3 generate_pipeline.py [project_name] [expression_matrix] [hub_genes] [out_directory] --run False</code></br>
<code>$sh [out_directory]/sjaracne_[project_name]_scripts_/00_pipeline_[project_name].sh</code>

## Running on Clusters

<code>$python3 generate_pipeline.py [project_name] [expression_matrix] [hub_genes] [out_directory] --run False --host CLUSTER</code></br>
<p>To run the pipeline on a cluster, use the script files under the scripts directory and submit the scripts 02 to 04 to the clusters.</p>

### Example of Running the Scripts on IBM LSF Cluster
<p>Run the following command for each script file, wait until the jobs are done, and do the same thing for the next script file:<p>
<code>cat [script_file] | while read line; do bsub $line; done</code>

## Example of Running Signaling Network on a Single Machine

<code>$export PYTHON_PATH=$(which python3)</code></br>
<code>$export SJARACNE_PATH=$(pwd)</code></br>
<code>$python3 generate_pipeline.py SIG data/BRCA100.exp data/sig.txt data/output/ --run True</code></br>

## Example of Running Transcription Factor Network on a Single Machine

<code>$export PYTHON_PATH=$(which python3)</code></br>
<code>$export SJARACNE_PATH=$(pwd)</code></br>
<code>$python3 generate_pipeline.py TF data/BRCA100.exp data/tf.txt data/output/ --run True</code></br>

## Expected Output
<p>Expected output for the example data with 3 bootstraps is available under data/output/ directory.</p>

