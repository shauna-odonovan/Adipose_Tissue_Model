# Model Parameters

This document provides biological descriptions for each of the 14 parameters used in the adipose tissue metabolism model. Parameters are listed in the order they appear in the parameter vector used throughout the MATLAB code (e.g., in `Estimating_Model_Parameters_Script.m`, `Model_Error.m`, and `Flux_Predictions.m`).

---

## LPL Lipolysis Parameters

### 1. K_ad — LPL Lipolysis Rate Constant
- **Index in parameter vector:** 1
- **Default initial value:** 0.0096
- **Bounds:** [0.001, 0.1]

**Biological explanation:** K_ad is the rate constant governing lipoprotein lipase (LPL)-mediated hydrolysis of circulating plasma triglycerides (TG) within the capillary endothelium of adipose tissue. LPL is the primary enzyme responsible for cleaving triglycerides into free fatty acids (NEFA) and glycerol, enabling their uptake by adipocytes. The rate of LPL action is modulated by insulin via a delayed signaling pathway (see LPL_delay). K_ad directly scales the triglyceride flux equation: `D_TG = -K_ad * I_LPL * TG_art`, where higher K_ad values indicate greater efficiency of LPL-mediated TG hydrolysis per unit of arterial triglyceride and insulin-delayed signal.

---

### 2. LPL_delay — Insulin Delay for LPL Activation
- **Index in parameter vector:** 2
- **Default initial value:** 155
- **Bounds:** [100, 300]

**Biological explanation:** LPL_delay represents the time delay (in minutes) for insulin to activate lipoprotein lipase on the capillary endothelium. Insulin does not immediately increase LPL activity; rather, its stimulatory effect is transmitted through a cascade of intracellular signaling events (including transcription, translation, and translocation of LPL to the endothelial surface). This delay is modeled using a three-compartment chain (`insulin_delay.m`), which produces a smoothed, delayed insulin signal `I_LPL` that modulates LPL-mediated TG lipolysis. A larger LPL_delay corresponds to a slower insulin-mediated upregulation of LPL activity.

---

## Spillover Parameter

### 3. D_spill — Fractional Spillover Scaling Factor
- **Index in parameter vector:** 3
- **Default initial value:** 43.15
- **Bounds:** [38, 50]

**Biological explanation:** D_spill governs the fraction of NEFA released by LPL-mediated hydrolysis of plasma triglycerides that escapes immediate uptake by adipocytes and instead "spills over" into the venous circulation. During LPL action on circulating TG, not all released NEFA are captured by the adipose tissue; some diffuse directly into the venous outflow. The fractional spillover is modeled as inversely proportional to plasma insulin: `spill = (D_spill * I_b / I_art) / 100`, where I_b is basal insulin. At higher insulin concentrations, spillover is reduced because insulin promotes NEFA re-esterification and retention within adipocytes. D_spill scales the magnitude of this spillover phenomenon.

---

## Glucose Uptake Parameters

### 4. GLUT_1 — GLUT1 Transporter Rate
- **Index in parameter vector:** 4
- **Default initial value:** 0.0181
- **Bounds:** [0.001, 0.1]

**Biological explanation:** GLUT_1 represents the rate constant for insulin-independent (basal) glucose uptake into adipocytes via the GLUT1 transporter. GLUT1 is constitutively expressed on the adipocyte membrane and facilitates passive diffusion of glucose along its concentration gradient. This basal uptake is independent of insulin and accounts for the component of glucose flux that persists even in the fasting state. The contribution to glucose flux is: `-GLUT_1 * G_art`, where G_art is the arterial glucose concentration.

---

### 5. GLUT_4 — GLUT4 Transporter Rate
- **Index in parameter vector:** 5
- **Default initial value:** 0.00055
- **Bounds:** [0, 0.1]

**Biological explanation:** GLUT_4 represents the rate constant for insulin-stimulated glucose uptake into adipocytes via the GLUT4 transporter. Unlike GLUT1, GLUT4 is stored in intracellular vesicles and is translocated to the cell membrane in response to insulin signaling. This insulin-dependent pathway is the primary mechanism by which postprandial glucose uptake into adipose tissue is enhanced. The contribution to glucose flux is: `-GLUT_4 * I_AT * G_art`, where I_AT is the delayed insulin signal in the adipose tissue compartment and G_art is the arterial glucose concentration.

---

### 6. AT_delay — Insulin Delay for Adipose Tissue
- **Index in parameter vector:** 6
- **Default initial value:** 5.38
- **Bounds:** [1, 45]

**Biological explanation:** AT_delay represents the time delay (in minutes) for insulin to exert its effects on adipose tissue, including stimulation of GLUT4 translocation, activation of re-esterification, and inhibition of intracellular lipolysis. Similar to LPL_delay, this is modeled using a three-compartment chain (`insulin_delay.m`), producing a delayed insulin signal `I_AT` that modulates insulin-dependent processes within the adipocyte. The delay reflects the time required for insulin receptor binding, intracellular signaling cascades (e.g., PI3K/Akt pathway), and downstream metabolic effects. AT_delay also serves as the glycerol delay parameter (GLY_delay) in some model configurations.

---

## Glycerol and Adipose Tissue Lipolysis Parameters

### 7. pgly — Glycerol Permeability Coefficient
- **Index in parameter vector:** 7
- **Default initial value:** 0.267
- **Bounds:** [0.001, 1]

**Biological explanation:** pgly is the permeability coefficient governing the passive diffusion of glycerol between the adipose tissue interstitial space and the plasma. Glycerol is a product of both LPL-mediated hydrolysis of circulating TG and intracellular lipolysis of stored TG. Since adipocytes lack glycerol kinase, glycerol produced within the adipocyte cannot be re-utilized and must be released into the circulation. The glycerol flux equation includes the term `pgly * (gly_at - GLY_art)`, representing diffusion-driven transport along the glycerol concentration gradient between adipose tissue and arterial plasma.

---

### 8. B_ATL — Basal Adipose Tissue Lipolysis Rate
- **Index in parameter vector:** 8
- **Default initial value:** 0.0278
- **Bounds:** [0.001, 1]

**Biological explanation:** B_ATL represents the basal rate of intracellular lipolysis within adipocytes, independent of insulin regulation. Adipose tissue lipolysis (ATL) is the process by which stored triglycerides are hydrolyzed into glycerol and NEFA by hormone-sensitive lipase (HSL) and adipose triglyceride lipase (ATGL). Even in the presence of insulin, there is a constitutive level of lipolysis that occurs. B_ATL defines this insulin-insensitive, constitutive component of intracellular TG breakdown, contributing a constant flux of glycerol and NEFA from the adipocyte.

---

### 9. ATL_max — Maximum Insulin-Suppressible Lipolysis Rate
- **Index in parameter vector:** 9
- **Default initial value:** 1.33
- **Bounds:** [0.01, 5]

**Biological explanation:** ATL_max represents the maximum rate of insulin-suppressible intracellular lipolysis in adipose tissue. In addition to the basal lipolysis (B_ATL), there is an insulin-sensitive component of lipolysis that is maximally active in the fasted state (low insulin) and suppressed when insulin rises. The total lipolysis rate is modeled as: `B_ATL + ATL_max / (1 + I_AT / K_ATL)`, where the second term represents the insulin-suppressible component that decreases as insulin increases. ATL_max defines the upper bound of this suppressible lipolysis.

---

### 10. K_ATL — Insulin Sensitivity for Lipolysis Inhibition
- **Index in parameter vector:** 10
- **Default initial value:** 1.71
- **Bounds:** [0.01, 5]

**Biological explanation:** K_ATL is the insulin concentration at which the insulin-suppressible component of adipose tissue lipolysis is reduced to half of its maximum value (ATL_max / 2). It represents the sensitivity of the lipolytic pathway to insulin-mediated inhibition. In the model, lipolysis is given by: `B_ATL + ATL_max / (1 + I_AT / K_ATL)`. A lower K_ATL indicates that lipolysis is more sensitive to insulin suppression (i.e., half-maximal suppression occurs at lower insulin concentrations), while a higher K_ATL indicates insulin resistance of the lipolytic pathway.

---

## NEFA and Re-esterification Parameters

### 11. re_ester — Re-esterification Rate Constant
- **Index in parameter vector:** 11
- **Default initial value:** 0.001
- **Bounds:** [0, 0.05]

**Biological explanation:** re_ester is the rate constant for the re-esterification of NEFA back into triglycerides within the adipocyte. Re-esterification is the process by which NEFA (either taken up from the plasma or released by intracellular lipolysis) are combined with glycerol-3-phosphate (G-3-P) to form stored triglycerides. This process is promoted by insulin, which enhances the availability of G-3-P through stimulation of glucose uptake and glycolysis. The re-esterification rate in the model is: `re_ester * I_AT * NEFA_at * G_3_P`, where I_AT is the delayed insulin signal, NEFA_at is the adipose tissue NEFA concentration, and G_3_P is the glycerol-3-phosphate concentration. This term couples NEFA dynamics to glucose metabolism through the shared intermediate G-3-P.

---

### 12. pnefa — NEFA Permeability Coefficient
- **Index in parameter vector:** 12
- **Default initial value:** 0.0532
- **Bounds:** [0, 1.5]

**Biological explanation:** pnefa is the permeability coefficient governing the transport of NEFA between the adipose tissue interstitial/cellular space and the plasma. NEFA can diffuse across the adipocyte membrane along their concentration gradient, facilitated by transport proteins. The NEFA flux includes the term `pnefa * (NEFA_at - NEFA_art_eff)`, where NEFA_art_eff is the effective arterial NEFA concentration adjusted for spillover NEFA from LPL hydrolysis. This represents the passive, concentration-gradient-driven exchange of NEFA between the adipose tissue and plasma compartments.

---

### 13. G3P_delay — G-3-P Production Delay
- **Index in parameter vector:** 13
- **Default initial value:** 60
- **Bounds:** [1, 100]

**Biological explanation:** G3P_delay represents the time delay (in minutes) for the conversion of glucose to glycerol-3-phosphate (G-3-P) within adipocytes. After glucose is taken up by the cell via GLUT1 and GLUT4 transporters, it is metabolized through glycolysis to produce G-3-P, which is a critical substrate for the re-esterification of NEFA into triglycerides. This delay is modeled using a two-compartment chain (`Glucose_G3P.m`), where `tau = G3P_delay` governs the time constant of the delay. A larger G3P_delay means that G-3-P becomes available for re-esterification more slowly after glucose uptake, reflecting the time required for glycolytic processing.

---

### 14. glucose_use — Fraction of Glucose Used for G-3-P Production
- **Index in parameter vector:** 14
- **Default initial value:** 0.9
- **Bounds:** [0, 1]

**Biological explanation:** glucose_use represents the fraction of glucose taken up by the adipocyte that is directed toward the production of glycerol-3-phosphate (G-3-P) for triglyceride synthesis. Not all glucose taken up by adipocytes is used for G-3-P production; some is oxidized for energy or converted to lactate. This parameter scales the glucose uptake signal before it enters the G-3-P delay compartments: `usage * (GLUT_1 * G_art + GLUT_4 * I_AT * G_art)`, where usage = glucose_use. A value of 1.0 would mean that all glucose uptake contributes to G-3-P production, while lower values indicate that a proportion of glucose is diverted to other metabolic fates.

---

## Summary Table

| # | Parameter | Symbol | Description | Typical Range |
|---|-----------|--------|-------------|---------------|
| 1 | K_ad | K_ad | LPL lipolysis rate constant | 0.001 - 0.1 |
| 2 | LPL_delay | tau_LPL | Insulin delay for LPL activation (min) | 100 - 300 |
| 3 | D_spill | D_spill | Fractional spillover scaling factor | 38 - 50 |
| 4 | GLUT_1 | GLUT_1 | GLUT1 basal glucose uptake rate | 0.001 - 0.1 |
| 5 | GLUT_4 | GLUT_4 | GLUT4 insulin-stimulated glucose uptake rate | 0 - 0.1 |
| 6 | AT_delay | tau_AT | Insulin delay for adipose tissue (min) | 1 - 45 |
| 7 | pgly | p_gly | Glycerol permeability coefficient | 0.001 - 1 |
| 8 | B_ATL | B_ATL | Basal adipose tissue lipolysis rate | 0.001 - 1 |
| 9 | ATL_max | ATL_max | Max insulin-suppressible lipolysis rate | 0.01 - 5 |
| 10 | K_ATL | K_ATL | Insulin sensitivity for lipolysis inhibition | 0.01 - 5 |
| 11 | re_ester | r_ester | Re-esterification rate constant | 0 - 0.05 |
| 12 | pnefa | p_NEFA | NEFA permeability coefficient | 0 - 1.5 |
| 13 | G3P_delay | tau_G3P | G-3-P production delay (min) | 1 - 100 |
| 14 | glucose_use | f_glu | Fraction of glucose used for G-3-P production | 0 - 1 |

---

*For further details on the model structure and equations, see `Model_description.pdf` and the associated publication: O'Donovan et al., PLOS Comp. Bio. (2019) 15(10): e1007400.*
