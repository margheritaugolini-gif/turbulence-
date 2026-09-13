# Model-Based Hopf Analysis

This workflow fits an SC-coupled Hopf model to empirical BOLD functional connectivity and computes model-based information capacity and susceptibility.

## Inputs

Required local files:

```text
timeseries_data_all.mat
schaefer1000_7N_coords_decolab.csv
SC_schaefer1000_decolab.mat
```

If `SC_schaefer1000_decolab.mat` is absent, `load_schaefer1000_sc.m` can load:

```text
sc_schaefer1000_decolab.npy
```

`timeseries_data_all.mat` is created with:

```bash
python3 create_timeseries_data_schaefer1000.py
```

## Step 1: Fit The Hopf Model

In MATLAB:

```matlab
hopf_model_fit
```

This saves the main fit file:

```text
hopf_sc_fit_results_schaefer1000_G0_to_4p5.mat
```

and checkpoint files such as:

```text
hopf_checkpoint_schaefer1000_G0_to_4p5.mat
```

Main parameters in `hopf_model_fit.m`:

```matlab
Cfg.filt.lb = 0.008;
Cfg.filt.ub = 0.08;

Cfg.NR = 20;
Cfg.fit_r_range_mm = [8.13 33.82];

Cfg.nNodes = 1000;
Cfg.TR = 2.9;
Cfg.min_T = 10;

Cfg.Glower = 0;
Cfg.Gstep = 0.25;
Cfg.Gupper = 4.5;

Cfg.NSIM_FIT = 10;
Cfg.NSIM_FINAL = 100;
Cfg.run_final_simulations = false;

Cfg.a = -0.02;
Cfg.sig = 0.01;
Cfg.excTp = 0;
```

Output filenames are controlled by:

```matlab
Cfg.save_file = 'hopf_checkpoint_schaefer1000_G0_to_4p5.mat';
fit_results_file = 'hopf_sc_fit_results_schaefer1000_G0_to_4p5.mat';
```

## Step 2: Plot Hopf Fit Results

```matlab
hopf_model_fit_plots
```

This reads:

```matlab
fit_results_file = 'hopf_sc_fit_results_schaefer1000_G0_to_4p5.mat';
```

and creates:

```text
Figure1_OptimalG_Violins.png
Figure1B_OptimalG_Control_vs_Meditation.png
Figure2_Empirical_vs_Simulated_FC.png
Figure3_GroupFit_Error_vs_G.png
```

If you want to plot a different saved fit, edit:

```matlab
fit_results_file = 'your_fit_results_file.mat';
```

## Step 3: Compute Information Capacity And Susceptibility

```matlab
hopf_capacity_susceptibility
```

This reads the fitted Hopf file:

```matlab
fit_results_file = 'hopf_sc_fit_results_schaefer1000_G0_to_4p5.mat';
```

and saves:

```text
hopf_subjectwise_IC_SUS.mat
```

Main parameters in `hopf_capacity_susceptibility.m`:

```matlab
Cfg.TR = 2.9;
Cfg.filt.lb = 0.008;
Cfg.filt.ub = 0.08;
Cfg.NSIM = 10;
Cfg.a = -0.02;
Cfg.sig = 0.01;
Cfg.burn_in_seconds = 2000;
```

## Step 4: Plot Capacity And Susceptibility

```matlab
hopf_capacity_susceptibility_violin_plots
```

This reads:

```text
hopf_subjectwise_IC_SUS.mat
```

and creates violin plots for information capacity and susceptibility across conditions.

## Model-Based Run Order

```matlab
hopf_model_fit
hopf_model_fit_plots
hopf_capacity_susceptibility
hopf_capacity_susceptibility_violin_plots
```

## Notes

The Hopf fitting step is computationally expensive. Increase or decrease `Cfg.NSIM_FIT`, `Cfg.Gstep`, and the G range depending on available compute time and required precision.

The published default pipeline uses:

```matlab
Cfg.Glower = 0;
Cfg.Gstep = 0.25;
Cfg.Gupper = 4.5;
```

If you use an extended or merged G search, update `fit_results_file` in the plotting and capacity/susceptibility scripts.
