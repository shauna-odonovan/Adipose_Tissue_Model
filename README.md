# Adipose Tissue Model

MATLAB implementation of a physiology-based mathematical model of postprandial adipose tissue metabolism.

## Overview

This model simulates the metabolic fluxes through human adipose tissue in the postprandial (fed) state. It captures the key processes of:

- Triglyceride (TG) lipolysis by lipoprotein lipase (LPL)
- Glucose uptake via GLUT1 and GLUT4 transporters
- Intracellular lipolysis and re-esterification
- Glycerol and NEFA transport

The model is described in detail in:

> O'Donovan, S.D., et al. *A physiology-based mathematical model of whole-body postprandial adipose tissue metabolism.* PLOS Computational Biology (2019) 15(10): e1007400.

## Repository Contents

| File | Description |
|------|-------------|
| `Parameter_Fit_Adipose_Model.m` | Main parameter fitting script |
| `Estimating_Model_Parameters_Script.m` | Parameter estimation with initial values and bounds |
| `Flux_Predictions.m` | Flux prediction functions |
| `Flux_Predictions_Script.m` | Script to run flux predictions with default parameters |
| `Glucose_G3P.m` | Glucose and G-3-P delay model |
| `Glucose_G3P_NEFA_ODE.m` | ODE system for glucose, G-3-P, and NEFA dynamics |
| `Glycerol_ODE.m` | ODE for glycerol dynamics |
| `insulin_delay.m` | Three-compartment insulin delay model |
| `Model_Error.m` | Objective function for parameter fitting |
| `sample_data.mat` | Sample experimental data |
| `model_description.pdf` | Detailed model description |

## Model Parameters

The model uses **14 parameters** that capture key physiological processes in adipose tissue metabolism. For complete biological descriptions, default values, units, and estimation bounds for each parameter, see **[PARAMETERS.md](PARAMETERS.md)**.

A summary table is provided below:

| # | Parameter | Default | Category |
|---|-----------|---------|----------|
| 1 | `K_ad` | 0.01 | LPL lipolysis |
| 2 | `LPL_delay` | 158.13 | LPL lipolysis |
| 3 | `D_spill` | 41.53 | Spillover |
| 4 | `GLUT_1` | 0.0159 | Glucose flux |
| 5 | `GLUT_4` | 0.00082 | Glucose flux |
| 6 | `AT_delay` | 13.85 | Glucose flux |
| 7 | `pgly` | 0.2705 | Glycerol flux |
| 8 | `B_ATL` | 0.0013 | Glycerol flux |
| 9 | `ATL_max` | 1.4373 | Glycerol flux |
| 10 | `K_ATL` | 1.49185 | Glycerol flux |
| 11 | `re_ester` | 0.000295 | NEFA |
| 12 | `pnefa` | 0.0478 | NEFA |
| 13 | `G3P_delay` | 31.74 | NEFA |
| 14 | `glucose_use` | 0.8803 | NEFA |

## Further Reading

- **[PARAMETERS.md](PARAMETERS.md)** — Full biological explanations for all 14 model parameters
- **[model_description.pdf](model_description.pdf)** — Detailed mathematical description of the model, parameter estimation methodology, and flux prediction framework

## License

This project is licensed under the terms found in the [LICENSE](LICENSE) file.
