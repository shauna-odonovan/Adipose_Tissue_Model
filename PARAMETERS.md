# Adipose Tissue Model Parameters

This document describes the parameters used in the MATLAB implementation of the adipose tissue model from O'Donovan et al. PLOS Comp. Bio. (2019) 15(10): e1007400.

## Parameter List

The model includes 14 parameters, indexed in the order they appear in the `parameters` vector:

| Index | Name       | Description | Units | Equation in paper |
|-------|------------|-------------|-------|-------------------|
| 1 | `K_ad` | Adipose tissue lipoprotein lipase (LPL) activity coefficient | mL/100mL tissue/min | Eq. X |
| 2 | `LPL_delay` | Delay time for insulin effect on LPL activity | min | Eq. Y |
| 3 | `D_spill` | Fractional spillover parameter | % | Eq. Z |
| 4 | `GLUT_1` | Basal glucose transporter activity | mL/100mL tissue/min | Eq. A |
| 5 | `GLUT_4` | Insulin-stimulated glucose transporter activity | mL/100mL tissue/min | Eq. B |
| 6 | `AT_delay` | Delay time for insulin effect on adipose tissue glucose uptake | min | Eq. C |
| 7 | `pgly` | Glycerol permeability coefficient | mL/100mL tissue/min | Eq. D |
| 8 | `B_ATL` | Basal adipose tissue lipolysis rate | µmol/100mL tissue/min | Eq. E |
| 9 | `ATL_max` | Maximum adipose tissue lipolysis rate | µmol/100mL tissue/min | Eq. F |
|10 | `K_ATL` | Michaelis constant for adipose tissue lipolysis | µmol/mL | Eq. G |
|11 | `re_ester` | Re-esterification rate coefficient | dimensionless | Eq. H |
|12 | `pnefa` | NEFA permeability coefficient | mL/100mL tissue/min | Eq. I |
|13 | `G3P_delay` | Delay time for G-3-P production from glucose | min | Eq. J |
|14 | `glucose_use` | Fraction of glucose used for G-3-P synthesis | dimensionless | Eq. K |

*Note: The exact equation references and units may need verification from the original publication.*

## Default Values and Ranges

Default parameter values used in the model fitting are not explicitly listed in the code. Refer to the `sample_data.mat` and the model fitting scripts for typical values.

## References

- O'Donovan, S. D., et al. (2019). A computational model of postprandial adipose tissue lipid metabolism derived using human arteriovenous stable isotope tracer data. *PLOS Computational Biology*, 15(10), e1007400. https://doi.org/10.1371/journal.pcbi.1007400