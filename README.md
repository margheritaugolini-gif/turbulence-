# Metastability And Turbulence Analysis

This repository contains MATLAB/Python scripts for Schaefer-1000 BOLD time-series analyses using both model-free measures and model-based Hopf simulations.

The workflows are intentionally separated:

- Model-free measures: see [README_model_free.md](README_model_free.md)
- Model-based Hopf analyses: see [README_model_based.md](README_model_based.md)

## Repository Contents

Core preprocessing:

- `create_timeseries_data_schaefer1000.py`
- `check_schaefer1000_alignment.py`
- `BOLD2hilbertupdated.m`
- `load_schaefer_coords.m`
- `load_schaefer1000_sc.m`
- `readNPY.m`
- `Violin.m`

Model-free analysis scripts:

- `metastability_kuramoto.m`
- `metastability_kuramoto_helper.m`
- `whole_brain_turbulence.m`
- `whole_brain_turbulence_helper.m`
- `whole_brain_turbulence_violin_plots.m`
- `network_node_turbulence.m`
- `network_node_turbulence_helper.m`
- `network_turbulence_spider_plots.m`
- `information_cascade.m`
- `information_cascade_violin_plots.m`

Model-based Hopf scripts:

- `hopf_model_fit.m`
- `hopf_model_fit_plots.m`
- `hopf_capacity_susceptibility.m`
- `hopf_capacity_susceptibility_violin_plots.m`

## Data Policy

Raw time series and generated MATLAB result files are not intended to be committed to GitHub. They are ignored by `.gitignore`.

Expected local input data:

- `time_series/`
- `schaefer1000_7N_coords_decolab.csv`
- `SC_schaefer1000_decolab.mat` or `sc_schaefer1000_decolab.npy`

Generated local files include:

- `timeseries_data_all.mat`
- `preprocessed_data.mat`
- `phase_for_turbulence.mat`
- `turbulence_results.mat`
- `node_turbulence_results.mat`
- `information_cascade_results.mat`
- `hopf_sc_fit_results_*.mat`
- `hopf_subjectwise_IC_SUS.mat`

Included partial result:

- `hopf_checkpoint_schaefer1000_G0_to_4p5.mat`

## Software

MATLAB is required for the main analyses. The scripts use functions from standard MATLAB toolboxes, including signal processing and statistics functions such as `butter`, `filtfilt`, `hilbert`, `signrank`, `friedman`, and `tiedrank`.

Python is used only to assemble/check the input `.mat` file:

```bash
python3 create_timeseries_data_schaefer1000.py
python3 check_schaefer1000_alignment.py
```

Python dependencies:

```bash
pip install numpy scipy pandas
```

## Quick Start

From MATLAB, run scripts from the repository root:

```matlab
cd('/path/to/metastability')
```

For model-free analyses, follow [README_model_free.md](README_model_free.md).

For model-based Hopf analyses, follow [README_model_based.md](README_model_based.md).
