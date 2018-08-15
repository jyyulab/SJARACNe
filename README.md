# SJARACNE

## Download

<code>git clone https://github.com/jyyulab/SJARACNE.git</code>

## Requirements

<p> Python 3.6.1 </p></br>

<p>For python dependencies look at the dependencies.txt file</p>


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
<p>The command will run scripts 02-04 automatically and generate the final results.</p>
<p>To generate only scripts run the command as follows:</p>
<code>$python3 generate_pipeline.py [project_name] [expression_matrix] [hub_genes] [out_directory] --run False</code></br>

#### Running on Single Machine
<code>$python3 generate_pipeline.py [project_name] [expression_matrix] [hub_genes] [out_directory]</code></br>

OR
<code>$python3 generate_pipeline.py [project_name] [expression_matrix] [hub_genes] [out_directory] --run False</code></br>
<code>$sh [out_directory]/sjaracne_[project_name]_scripts_/00_pipeline_[project_name].sh</code>


#### Running on Clusters

<code>$python3 generate_pipeline.py [project_name] [expression_matrix] [hub_genes] [out_directory] --run False</code></br>
<p>To run the pipeline on a cluster, use the script files under the scripts directory and submit the scripts 01 to 04 to the clusters.</p>


### Options

<code>--bootstrap, default=100, Number of bootstrap networks.</code></br>
<code>--c_threshold, default=1e-5, P-value threshold in building consensus network.</code></br>
<code>--p_threshold, default=1e-7, P-value threshold in building bootstrap netwroks.</code></br>
<code>--depth, default=40, help=Maximum partitioning depth.</code></br>
<code>--run, default=True, help=Whether run the pipeline or just generate and stop.</code></br>

## Example of Running Signaling Network on a Single Machine

<code>$export PYTHON_PATH=[path_to_python3]</code></br>
<code>$export SJARACNE_PATH=[path_to_package]</code></br>
<code>$python3 generate_pipeline.py SIG data/BRCA100.exp data/sig.txt data/</code></br>

## Example of Running Transcription Factor Network on a Single Machine

<code>$export PYTHON_PATH=[path_to_python3]</code></br>
<code>$export SJARACNE_PATH=[path_to_package]</code></br>
<code>$python3 generate_pipeline.py TF data/BRCA100.exp data/tf.txt data/</code></br>

