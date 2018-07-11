# SJARACNE

## Download

<code>git clone https://github.com/jyyulab/SJARACNE.git</code>

## Requirements

<p> Python 3.6.1 </p></br>


## Running pipeline

### EASY RUN

<code>$python3 generate_pipeline.py [project_name] [expression_matrix] [hub_genes] [out_directory]</code></br>

### Options

<code>--bootstrap, default=100, Number of bootstrap networks.</code>
<code>--c_threshold, default=1e-5, P-value threshold in building consensus network.</code>
<code>--p_threshold, default=1e-7, P-value threshold in building bootstrap netwroks.</code>
<code>--depth, default=40, help=Maximum partitioning depth.</code>
