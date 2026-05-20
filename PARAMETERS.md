# Model Parameters

This document provides biological explanations for all 14 model parameters used in the adipose tissue metabolic model. These parameters are estimated from arteriovenous difference data using the `lsqnonlin` algorithm (see `Estimating_Model_Parameters_Script.m` and `Parameter_Fit_Adipose_Model.m`).

## Parameter Overview

| # | Symbol | Description | Group | Initial Value | Lower Bound | Upper Bound |
|---|--------|-------------|-------|---------------|-------------|-------------|
| 1 | K_ad | LPL lipolysis rate constant | LPL lipolysis | 0.00958 | 0.001 | 0.1 |
| 2 | LPL_delay | Insulin delay for LPL activation (mins) | LPL lipolysis | 155 | 100 | 300 |
| 3 | D_spill | Fractional spillover scaling factor | Spillover | 43.15 | 38 | 50 |
| 4 | GLUT_1 | GLUT1-mediated glucose uptake rate | Glucose flux | 0.01809 | 0.001 | 0.1 |
| 5 | GLUT_4 | GLUT4-mediated insulin-sensitive glucose uptake rate | Glucose flux | 0.000552 | 0 | 0.1 |
| 6 | AT_delay | Insulin delay for adipose tissue (mins) | Glucose flux | 5.384 | 1 | 45 |
| 7 | pgly | Glycerol transport permeability coefficient | Glycerol flux | 0.2675 | 0.001 | 1 |
| 8 | B_ATL | Basal adipose tissue lipolysis rate | Glycerol flux | 0.02781 | 0.001 | 1 |
| 9 | ATL_max | Maximum insulin-suppressible lipolysis rate | Glycerol flux | 1.329 | 0.01 | 5 |
| 10 | K_ATL | Insulin sensitivity constant for lipolysis suppression | Glycerol flux | 1.713 | 0.01 | 5 |
| 11 | re_ester | NEFA re-esterification rate | NEFA flux | 0.001 | 0 | 0.05 |
| 12 | pnefa | NEFA transport permeability coefficient | NEFA flux | 0.05319 | 0 | 1.5 |
| 13 | G3P_delay | G-3-P delay time constant (mins) | NEFA flux | 60 | 1 | 100 |
| 14 | glucose_use | Proportion of glucose used for G-3-P production | NEFA flux | 0.9 | 0 | 1 |

---

## Detailed Parameter Descriptions

### LPL Lipolysis Parameters

#### 1. K_ad — LPL Lipolysis Rate Constant

The rate constant governing the hydrolysis of circulating triglyceride (TG) by lipoprotein lipase (LPL) on the endothelial surface of adipose tissue capillaries. LPL is the primary enzyme responsible for the breakdown of chylomicron and VLDL triglycerides into free fatty acids and glycerol for uptake by adipocytes. In the model, the triglyceride arteriovenous difference is computed as `D_TG = -K_ad × I_LPL × TG_art`, where `I_LPL` is the delayed insulin signal that activates LPL. This parameter thus scales the overall efficiency of LPL-mediated TG clearance from plasma by adipose tissue.

#### 2. LPL_delay — Insulin Delay for LPL Activation

The time delay (in minutes) representing the lag between changes in plasma insulin concentration and the corresponding activation of LPL on the capillary endothelium. Insulin stimulates LPL synthesis and translocation to the endothelial surface, but this process is not instantaneous. The model uses a three-compartment delay (implemented in `insulin_delay.m`) to produce a smooth, physiologically realistic delayed insulin signal `I_LPL` from the measured arterial insulin. A longer delay implies slower postprandial LPL activation and hence delayed TG clearance.

---

### Spillover Parameter

#### 3. D_spill — Fractional Spillover Scaling Factor

A scaling factor that determines the fractional spillover of NEFA released by LPL-mediated lipolysis of plasma triglycerides. Not all fatty acids released by LPL are retained by the adipose tissue; a fraction "spills over" into the venous circulation without being taken up by adipocytes. The model computes fractional spillover as `spill = (D_spill × I_b / I_PL) / 100`, where `I_b` is basal insulin and `I_PL` is the current plasma insulin. This inverse relationship with insulin reflects the physiological observation that insulin promotes fatty acid entrapment within adipose tissue, so spillover is higher in the fasted (low insulin) state and lower postprandially.

---

### Glucose Flux Parameters

#### 4. GLUT_1 — GLUT1-Mediated Glucose Uptake Rate

The rate constant for insulin-independent (basal) glucose uptake into adipocytes via the GLUT1 transporter. GLUT1 is constitutively present on the adipocyte plasma membrane and facilitates basal glucose transport regardless of insulin levels. In the model, the insulin-independent component of glucose flux is `GLUT_1 × G_art`, where `G_art` is the arterial glucose concentration. This parameter represents the basal, always-active pathway for glucose entry into adipose tissue.

#### 5. GLUT_4 — GLUT4-Mediated Insulin-Sensitive Glucose Uptake Rate

The rate constant for insulin-stimulated glucose uptake into adipocytes via the GLUT4 transporter. GLUT4 is stored in intracellular vesicles and translocates to the cell membrane in response to insulin signaling, dramatically increasing glucose uptake. In the model, the insulin-dependent glucose flux component is `GLUT_4 × I_AT × G_art`, where `I_AT` is the delayed insulin signal in the adipose tissue compartment. This parameter captures the postprandial increase in glucose uptake driven by the insulin-mediated recruitment of GLUT4 transporters.

#### 6. AT_delay — Insulin Delay for Adipose Tissue (mins)

The time delay (in minutes) representing the lag between changes in plasma insulin concentration and the insulin signal reaching the intracellular adipose tissue compartment. This delay accounts for the time required for insulin to cross the endothelium, bind to its receptor, and initiate downstream signaling (including GLUT4 translocation and suppression of intracellular lipolysis). Like `LPL_delay`, it is implemented as a three-compartment delay in `insulin_delay.m`. The resulting delayed insulin signal `I_AT` modulates GLUT4-mediated glucose uptake, intracellular lipolysis suppression, and NEFA re-esterification.

---

### Glycerol Flux Parameters

#### 7. pgly — Glycerol Transport Permeability Coefficient

The permeability-surface area product governing the passive transport of glycerol between the adipose tissue interstitial space and the arterial plasma. Glycerol is a product of both LPL-mediated lipolysis (extracellular) and intracellular lipolysis of stored TG. Since adipocytes lack significant glycerol kinase activity, glycerol produced intracellularly must exit the cell and enter the circulation. In the model, glycerol efflux is driven by the concentration gradient between adipose tissue glycerol and the effective arterial glycerol concentration, scaled by `pgly`.

#### 8. B_ATL — Basal Adipose Tissue Lipolysis Rate

The basal (insulin-independent) rate of intracellular lipolysis of stored triglycerides within adipocytes. Even in the absence of insulin, hormone-sensitive lipase (HSL) and adipose triglyceride lipase (ATGL) continuously hydrolyze stored TG, releasing glycerol and fatty acids. `B_ATL` represents this constitutive lipolytic activity, which contributes to the fasting glycerol flux. In the model, the total intracellular lipolysis rate is `B_ATL + ATL_max / (1 + I_AT / K_ATL)`, where the second term is the insulin-suppressible component.

#### 9. ATL_max — Maximum Insulin-Suppressible Lipolysis Rate

The maximum rate of intracellular lipolysis that can be suppressed by insulin. This parameter represents the dynamic range of insulin-mediated suppression of adipose tissue lipolysis. In the fed state, rising insulin levels suppress HSL and ATGL activity, reducing the release of glycerol and NEFA from stored TG. The insulin-suppressible lipolysis component is modeled as `ATL_max / (1 + I_AT / K_ATL)`, which decreases toward zero as insulin rises. `ATL_max` thus quantifies the total capacity of insulin-suppressible lipolysis.

#### 10. K_ATL — Insulin Sensitivity Constant for Lipolysis Suppression

The insulin concentration at which the insulin-suppressible component of intracellular lipolysis is reduced to half of its maximum value (`ATL_max / 2`). This parameter reflects the sensitivity of adipose tissue lipolysis to insulin: a lower `K_ATL` means that less insulin is needed to suppress lipolysis (greater insulin sensitivity), while a higher `K_ATL` indicates insulin resistance at the level of intracellular lipolysis. The Michaelis-Menten form `ATL_max / (1 + I_AT / K_ATL)` ensures that suppression is graded and saturable.

---

### NEFA Flux Parameters

#### 11. re_ester — NEFA Re-esterification Rate

The rate constant governing the re-esterification of intracellular NEFA back into stored triglycerides. Re-esterification is a critical metabolic pathway in adipose tissue that cycles fatty acids back into TG storage, and it is stimulated by insulin. The process requires glycerol-3-phosphate (G-3-P), which is derived from glucose metabolism. In the model, the re-esterification rate is `re_ester × I_AT × NEFA_at × G_3_P`, where `I_AT` is the delayed adipose insulin signal, `NEFA_at` is the adipose tissue NEFA concentration, and `G_3_P` is the adipose G-3-P concentration. This ensures that re-esterification is enhanced by both insulin and the availability of G-3-P from glucose uptake.

#### 12. pnefa — NEFA Transport Permeability Coefficient

The permeability-surface area product governing the passive transport (diffusion) of NEFA between the adipose tissue interstitial space and the plasma. NEFA can cross the adipocyte membrane and endothelium along their concentration gradient. In the model, NEFA efflux from adipose tissue is driven by the gradient between adipose NEFA and the effective plasma NEFA concentration (which includes NEFA from spillover), scaled by `pnefa`. This parameter captures the net passive exchange of NEFA between tissue and plasma compartments.

#### 13. G3P_delay — G-3-P Delay Time Constant (mins)

The time constant (in minutes) for the delay between glucose uptake by adipose tissue and the appearance of glycerol-3-phosphate (G-3-P) available for re-esterification. After glucose enters the adipocyte, it must be phosphorylated to glucose-6-phosphate and then metabolized through glycolysis to produce G-3-P. This metabolic processing introduces a time lag, which the model captures using a two-compartment delay (implemented in `Glucose_G3P.m`). A longer `G3P_delay` means that G-3-P becomes available more slowly after glucose uptake, which in turn delays the onset of NEFA re-esterification.

#### 14. glucose_use — Proportion of Glucose Used for G-3-P Production

The fraction of glucose taken up by adipose tissue that is directed toward the production of glycerol-3-phosphate (G-3-P) for triglyceride synthesis (re-esterification). Not all glucose taken up by adipocytes is used for G-3-P production; some is oxidized for energy via the TCA cycle or used in other metabolic pathways. This parameter scales the glucose flux entering the G-3-P delay compartment: `glucose_use × (GLUT_1 × G_art + GLUT_4 × I_AT × G_art)`. A value of 1.0 would indicate that all glucose is used for G-3-P production, while lower values indicate diversion of glucose to other metabolic fates.

---

## Parameter Groups in the Code

The 14 parameters are organized into functional groups within the MATLAB code:

- **LPL lipolysis** (`parameters_LPL`): `K_ad` (1), `LPL_delay` (2)
- **Spillover** (`parameters_spill`): `D_spill` (3)
- **Glucose flux** (`parameters_glu`): `GLUT_1` (4), `GLUT_4` (5), `AT_delay` (6)
- **Glycerol flux / Adipose tissue lipolysis** (`parameters_ATL`): `pgly` (7), `B_ATL` (8), `ATL_max` (9), `K_ATL` (10)
- **NEFA flux** (`parameters_nefa`): `re_ester` (11), `pnefa` (12), `G3P_delay` (13), `glucose_use` (14)

## References

O'Donovan, S. D., et al. A computational model of postprandial adipose tissue lipid metabolism derived using human arteriovenous stable isotope tracer data. *PLOS Computational Biology* (2019) 15(10): e1007400.