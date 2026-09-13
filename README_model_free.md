# Model-Free Measures

This workflow computes model-free measures from empirical BOLD time series:

- Global metastability from the Kuramoto order parameter
- Whole-brain turbulence across lambda scales
- Network/node-level turbulence
- Information cascade across adjacent lambda scales

## Inputs

Place Schaefer-1000 time-series files under:

```text
time_series/sub-*/
```

The expected files are selected by `create_timeseries_data_schaefer1000.py` using:

```python
ATLAS_SUFFIX = "schaefer1000-7networks_ts.txt"
CONDITIONS = ["memory", "counting"] + [f"j{i}" for i in range(1, 9)]
```

Also keep this coordinate file in the repository root:

```text
schaefer1000_7N_coords_decolab.csv
```

## Step 1: Build MATLAB Input

Run from the repository root:

```bash
python3 create_timeseries_data_schaefer1000.py
```

This creates:

```text
timeseries_data_all.mat
```

Optional alignment check:

```bash
python3 check_schaefer1000_alignment.py
```

## Step 2: Preprocess BOLD To Hilbert Phase

In MATLAB:

```matlab
BOLD2hilbertupdated
```

This creates:

```text
preprocessed_data.mat
phase_for_turbulence.mat
```

Main parameters in `BOLD2hilbertupdated.m`:

```matlab
min_T = 10;
TR = 2.9;
low_cut = 0.008;
high_cut = 0.08;
excTp = 3;
```

## Step 3: Global Metastability

```matlab
metastability_kuramoto
```

This appends `OP_all`, `Sync`, and `Meta` to `preprocessed_data.mat`, then saves:

```text
metastability_statistics.mat
```

Main file to edit:

```text
metastability_kuramoto.m
```

The condition order is assumed to be:

```text
Memory, Counting, J1, J2, J3, J4, J5, J6, J7, J8
```

## Step 4: Whole-Brain Turbulence

```matlab
whole_brain_turbulence
```

This creates:

```text
turbulence_results.mat
```

Main lambda range is set in `whole_brain_turbulence_helper.m`:

```matlab
LAMBDA = 0.01:0.01:0.30;
```

To create whole-brain turbulence violin/statistics plots:

```matlab
whole_brain_turbulence_violin_plots
```

Main parameters in `whole_brain_turbulence_violin_plots.m`:

```matlab
target_lambdas = [0.01 0.12 0.30];
alpha = 0.05;
```

## Step 5: Network/Node-Level Turbulence

To compute node-level turbulence directly from Hilbert phases:

```matlab
network_node_turbulence
```

This creates:

```text
node_turbulence_results.mat
```

The fixed node-level lambda is set in `network_node_turbulence_helper.m`:

```matlab
lambda = 0.12;
```

To compute network-level spider/violin plots from `turbulence_results.mat`:

```matlab
network_turbulence_spider_plots
```

The network plotting lambda is set in `network_turbulence_spider_plots.m`:

```matlab
target_lambda = 0.12;
```

## Step 6: Information Cascade

```matlab
information_cascade
```

This computes information cascade from `turbulence_results.mat`, creates violin plots, and saves:

```text
information_cascade_results.mat
```

The cascade calculation uses all available lambdas in `output.LAMBDA` and tests adjacent lambda transitions:

```matlab
lambda(1) -> lambda(2)
lambda(2) -> lambda(3)
...
lambda(end-1) -> lambda(end)
```

The current significance threshold is:

```matlab
sig = pp(:) < 0.05;
```

To recreate only the information-cascade violin plots later:

```matlab
information_cascade_violin_plots
```

## Model-Free Run Order

```matlab
BOLD2hilbertupdated
metastability_kuramoto
whole_brain_turbulence
whole_brain_turbulence_violin_plots
network_node_turbulence
network_turbulence_spider_plots
information_cascade
information_cascade_violin_plots
```
