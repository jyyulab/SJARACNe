# SJARACNE

## Download

<code>git clone https://github.com/jyyulab/SJARACNE.git</code>

## Requirements

<p> Python 3.6.1 </p></br>


## Running pipeline

### Set Environemnt

<code>$export PYTHON_PATH=[path_to_python3]</code></br>
<code>$export SJARACNE_PATH=[path_to_package]</code></br>

### EASY RUN

<code>$python3 generate_pipeline.py [project_name] [expression_matrix] [hub_genes] [out_directory]</code></br>

<p>The above command will create 4 directories under the provided out_directory parameter as follows:</p></br>
<code>[out_directory]/sjaracne_[project_name]_out_</code></br>
<code>[out_directory]/sjaracne_[project_name]_log_</code></br>
<code>[out_directory]/sjaracne_[project_name]_out_.final</code></br>
<code>[out_directory]/sjaracne_[project_name]_scripts_</code></br>
<p>There will be shell script files corresponding to the provided input files in the scripts directory in the following order:</p>
<code>00_cleanup_[project_name].sh</code></br>
<code>00_pipeline_[project_name].sh</code></br>
<code>01_prepare_[project_name].sh</code></br>
<code>02_bootstrap_[project_name].sh</code></br>
<code>03_getconsensusnetwork_[project_name].sh</code></br>
<code>04_getenhancedconsensusnetwork_[project_name].sh</code></br>

#### Running on Single Machine
<p>To run the pipeline on a single machine, only run the 00_pipeline_[project_name].sh from the scripts directory</p>
<code>sh 00_pipeline_[project_space].sh</code></br>

#### Running on Clusters

<p>To run the pipeline on a cluster, use the script files under the scripts directory and submit the scripts 01 to 04 to the clusters.</p>

### Options

<code>--bootstrap, default=100, Number of bootstrap networks.</code></br>
<code>--c_threshold, default=1e-5, P-value threshold in building consensus network.</code></br>
<code>--p_threshold, default=1e-7, P-value threshold in building bootstrap netwroks.</code></br>
<code>--depth, default=40, help=Maximum partitioning depth.</code></br>

## Example of Running Signaling Network on a Single Machine

<code>$export PYTHON_PATH=[path_to_python3]</code></br>
<code>$export SJARACNE_PATH=[path_to_package]</code></br>
<code>$python generate_pipeline.py SIG data/BRCA100.exp data/sig.txt data/</code></br>
<code>$cd data/sjaracne_SIG_scripts_/</code></br>
<code>$sh 00_pipeline_SIG.sh</code></br>

## Example of Running Transcription Factor Network on a Single Machine

<code>$export PYTHON_PATH=[path_to_python3]</code></br>
<code>$export SJARACNE_PATH=[path_to_package]</code></br>
<code>$python generate_pipeline.py TF data/BRCA100.exp data/tf.txt data/</code></br>
<code>$cd data/sjaracne_TF_scripts_/</code></br>
<code>$sh 00_pipeline_TF.sh</code></br>
