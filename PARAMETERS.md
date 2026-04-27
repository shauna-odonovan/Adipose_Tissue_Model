# Model Parameters

This document provides biological explanations for all 14 parameters used in the adipose tissue metabolic model implemented in the MATLAB code. These parameters are estimated from human arteriovenous stable isotope tracer data, as described in O'Donovan et al. *A computational model of postprandial adipose tissue lipid metabolism derived using human arteriovenous stable isotope tracer data.* PLOS Comp. Bio. (2019) 15(10): e1007400.

The parameters are defined in `Estimating_Model_Parameters_Script.m` and used throughout the model ODE files (`Glucose_G3P_NEFA_ODE.m`, `Glycerol_ODE.m`, `Parameter_Fit_Adipose_Model.m`).

---

## 1. K_ad — LPL-Mediated Triglyceride Hydrolysis Rate Constant

**Code reference:** `p_opt(1)`, `parameters(1)`  
**Typical range:** 0.001–0.1  
**Unit context:** min⁻¹ (normalised by plasma TG concentration)

K_ad is the rate constant governing lipoprotein lipase (LPL)–mediated hydrolysis of circulating triglyceride-rich lipoproteins (chylomicrons and VLDL) at the adipose tissue endothelial surface. In the model, the rate of triglyceride extraction from arterial plasma is computed as `D_TG = -K_ad × I_LPL × TG_art`, where `I_LPL` is the delayed insulin signal acting on LPL and `TG_art` is the arterial triglyceride concentration. Biologically, LPL is anchored to the luminal surface of capillary endothelial cells and is the primary enzyme responsible for hydrolysing circulating triglycerides into fatty acids for tissue uptake. Its activity is upregulated by insulin through transcriptional and post-translational mechanisms, which is captured by the insulin delay term `I_LPL`.

---

## 2. LPL_delay — Insulin Action Delay on LPL

**Code reference:** `p_opt(2)`, `parameters(2)`  
**Typical range:** 100–300 min  
**Unit:** minutes

LPL_delay represents the time delay between a rise in plasma insulin and the subsequent activation of LPL on the adipose tissue endothelium. This delay is implemented as a three-compartment insulin delay system (see `insulin_delay.m`), where the delay time τ = LPL_delay governs the rate at which the delayed insulin signal `I_LPL` approaches the plasma insulin level `I_pl`. Biologically, this reflects the time required for insulin to stimulate LPL synthesis in adipocytes, secretion of the enzyme, translocation to the capillary endothelium, and activation — processes that are known to operate on a timescale of several hours rather than minutes.

---

## 3. D_spill — Fatty Acid Spillover Fraction

**Code reference:** `p_opt(3)`, `parameters(3)`  
**Typical range:** 38–50  
**Unit:** percentage (%)

D_spill represents the fraction of fatty acids released by LPL-mediated triglyceride hydrolysis that are not immediately taken up by adipocytes but instead "spill over" directly into the venous plasma. In the model, the spillover fraction is computed as `spill = (D_spill × (I_b / I_PL)) / 100`, where `I_b` is basal insulin and `I_PL` is the current plasma insulin level. Biologically, spillover occurs because not all fatty acids liberated by LPL at the endothelial surface are captured by adjacent adipocytes; a proportion escapes into the circulation. This phenomenon is well documented in human tracer studies and is influenced by the local concentration gradient and the rate of fatty acid uptake by adipocytes.

---

## 4. GLUT_1 — Insulin-Independent Glucose Transport Rate

**Code reference:** `p_opt(4)`, `parameters(4)`  
**Typical range:** 0.001–0.1  
**Unit context:** min⁻¹ (normalised by arterial glucose)

GLUT_1 is the rate constant for insulin-independent (basal) glucose uptake into adipocytes, mediated primarily by GLUT1 glucose transporters constitutively expressed on the adipocyte plasma membrane. In the model, basal glucose uptake is computed as `GLUT_1 × G_art`, where `G_art` is the arterial glucose concentration. Biologically, GLUT1 transporters facilitate a low level of glucose entry into adipocytes regardless of insulin status, ensuring a basal supply of glucose for essential metabolic processes such as glycolysis and glycerol-3-phosphate production.

---

## 5. GLUT_4 — Insulin-Dependent Glucose Transport Rate

**Code reference:** `p_opt(5)`, `parameters(5)`  
**Typical range:** 0–0.1  
**Unit context:** min⁻¹ (normalised by arterial glucose and insulin)

GLUT_4 is the rate constant for insulin-stimulated glucose uptake into adipocytes, mediated by GLUT4 glucose transporters that translocate from intracellular vesicles to the plasma membrane upon insulin stimulation. In the model, insulin-dependent glucose uptake is computed as `GLUT_4 × I_AT × G_art`, where `I_AT` is the delayed insulin signal in adipose tissue. Biologically, GLUT4 translocation is the primary mechanism by which insulin acutely increases glucose uptake in adipose tissue; this process accounts for the majority of postprandial glucose disposal by adipocytes.

---

## 6. AT_delay — Adipose Tissue Insulin Action Delay

**Code reference:** `p_opt(6)`, `parameters(6)`  
**Typical range:** 1–45 min  
**Unit:** minutes

AT_delay represents the time delay between a rise in plasma insulin and the physiological response of adipose tissue (specifically, the insulin-stimulated glucose uptake via GLUT4). This delay is implemented using the same three-compartment insulin delay system as `LPL_delay` (see `insulin_delay.m`), with τ = AT_delay. Biologically, this delay captures the time required for insulin receptor binding, intracellular signalling cascade activation, GLUT4 vesicle translocation to the cell surface, and the resultant increase in glucose transport. The AT_delay is typically much shorter than LPL_delay because GLUT4 translocation occurs on a timescale of minutes rather than hours.

---

## 7. pgly — Glycerol Efflux Rate Constant

**Code reference:** `p_opt(7)`, `parameters(7)`  
**Typical range:** 0.001–1  
**Unit context:** min⁻¹ (permeability coefficient)

pgly is the permeability/transport coefficient governing the passive efflux of glycerol from adipocytes into the venous plasma. In the model (`Glycerol_ODE.m`), glycerol transport is driven by the concentration gradient between intracellular glycerol (`gly_at`) and arterial glycerol (`GLY_art`), plus the glycerol generated from LPL-mediated triglyceride hydrolysis. Biologically, adipocytes lack significant glycerol kinase activity and therefore cannot re-utilise the glycerol produced during lipolysis; instead, glycerol exits the cell via aquaporin-7 and other membrane channels down its concentration gradient, making it a key marker of lipolytic activity in the tissue.

---

## 8. B_ATL — Basal Adipose Tissue Lipolysis Rate

**Code reference:** `p_opt(8)`, `parameters(8)`  
**Typical range:** 0.001–1  
**Unit context:** μmol·kg⁻¹·min⁻¹

B_ATL represents the insulin-independent (basal) rate of intracellular triglyceride lipolysis in adipocytes. In the model, total lipolysis is given by `B_ATL + ATL_max / (1 + I_AT / K_ATL)`, where the first term `B_ATL` is the basal component that persists even when insulin is elevated. Biologically, basal lipolysis reflects constitutive activity of hormone-sensitive lipase (HSL) and adipose triglyceride lipase (ATGL) that is not fully suppressed by insulin. This basal rate ensures a continuous supply of fatty acids for energy metabolism even under conditions of low lipolytic stimulation.

---

## 9. ATL_max — Maximum Insulin-Suppressible Lipolysis Rate

**Code reference:** `p_opt(9)`, `parameters(9)`  
**Typical range:** 0.01–5  
**Unit context:** μmol·kg⁻¹·min⁻¹

ATL_max is the maximum rate of lipolysis that can be suppressed by insulin action in adipose tissue. In the model, the insulin-suppressible component of lipolysis is `ATL_max / (1 + I_AT / K_ATL)`, which decreases towards zero as insulin levels rise. Biologically, ATL_max represents the maximal capacity of insulin-regulated lipolysis driven by catecholamines and other lipolytic hormones. Under fasting or stress conditions with low insulin, lipolysis approaches `B_ATL + ATL_max`, while in the postprandial state with elevated insulin, the suppressible fraction is largely inhibited.

---

## 10. K_ATL — Insulin Sensitivity of Adipose Tissue Lipolysis

**Code reference:** `p_opt(10)`, `parameters(10)`  
**Typical range:** 0.01–5  
**Unit context:** insulin concentration units (mU/L)

K_ATL is the half-saturation constant for insulin's suppression of adipose tissue lipolysis. In the model, the insulin-suppressible lipolysis term `ATL_max / (1 + I_AT / K_ATL)` follows Michaelis-Menten kinetics, where K_ATL represents the insulin concentration at which lipolysis is half-maximally suppressed. Biologically, a lower K_ATL indicates greater insulin sensitivity (i.e., less insulin is required to suppress lipolysis), while a higher K_ATL indicates insulin resistance of the anti-lipolytic pathway. This parameter captures inter-individual variation in the responsiveness of adipocytes to insulin's anti-lipolytic effects.

---

## 11. re_ester — Re-esterification Rate Constant

**Code reference:** `p_opt(11)`, `parameters(11)`  
**Typical range:** 0–0.05  
**Unit context:** min⁻¹ (normalised by insulin, NEFA, and G-3-P)

re_ester is the rate constant for the re-esterification of intracellular non-esterified fatty acids (NEFA) back into triglycerides within adipocytes. In the model, re-esterification is computed as `re_ester × I_AT × NEFA_at × G_3_P`, where the rate depends on the delayed insulin signal `I_AT`, intracellular NEFA concentration `NEFA_at`, and glycerol-3-phosphate availability `G_3_P`. Biologically, re-esterification is the process by which fatty acids are esterified with glycerol-3-phosphate to reform triglycerides for storage in the adipocyte lipid droplet. Insulin promotes re-esterification by increasing glucose uptake and thus glycerol-3-phosphate production, and by activating key enzymes such as glycerol-3-phosphate acyltransferase (GPAT). This parameter is critical for the balance between fatty acid storage and release.

---

## 12. pnefa — NEFA Transport Coefficient

**Code reference:** `p_opt(12)`, `parameters(12)`  
**Typical range:** 0–1.5  
**Unit context:** min⁻¹ (permeability coefficient)

pnefa is the permeability/transport coefficient for the movement of non-esterified fatty acids (NEFA) across the adipocyte plasma membrane. In the model (`Glucose_G3P_NEFA_ODE.m`), NEFA transport is driven by the concentration gradient: `pnefa × (NEFA_at - NEFA_art_effective)`. Biologically, NEFA cross the adipocyte membrane via both passive diffusion (flip-flop across the lipid bilayer) and protein-facilitated transport (e.g., CD36/FAT, FATP). This parameter captures the net effective permeability of the membrane to fatty acids, which determines how rapidly intracellular fatty acid levels equilibrate with plasma concentrations.

---

## 13. G3P_delay — Glycerol-3-Phosphate Production Delay

**Code reference:** `p_opt(13)`, `parameters(13)`  
**Typical range:** 1–100 min  
**Unit:** minutes

G3P_delay is the time delay (τ) in the compartmental model linking glucose uptake to glycerol-3-phosphate (G-3-P) production (see `Glucose_G3P.m`). The model uses a two-compartment delay system where glucose is first converted to intermediate metabolites (`G_11`) and then to G-3-P, with the rate of each transition governed by `1/τ`. Biologically, this delay represents the time required for glucose to be taken up by the adipocyte, metabolised through glycolysis, and converted to glycerol-3-phosphate — a multi-step enzymatic process involving hexokinase, phosphofructokinase, and other glycolytic enzymes. G-3-P is essential for re-esterification of fatty acids into triglycerides.

---

## 14. glucose_use — Fraction of Glucose Directed to G-3-P Production

**Code reference:** `p_opt(14)`, `parameters(14)`  
**Typical range:** 0–1  
**Unit:** dimensionless (fraction)

glucose_use (referred to as `usage` in `Glucose_G3P.m`) represents the proportion of glucose taken up by adipocytes that is directed toward glycerol-3-phosphate (G-3-P) production via glycolysis. In the model, G-3-P production is computed as `usage × (GLUT_1 × G_art + GLUT_4 × I_AT × G_art)`. Biologically, not all glucose entering the adipocyte is used for G-3-P synthesis; a fraction is oxidised via the TCA cycle, used for lactate production, or directed to other biosynthetic pathways. The remaining fraction (glucose_use) that enters glycolysis generates dihydroxyacetone phosphate, which is reduced to G-3-P — the backbone molecule required for triglyceride re-esterification. This parameter is therefore central to coupling glucose metabolism with lipid storage.

---

## Parameter Summary Table

| # | Parameter | Symbol | Group | Initial Value | Lower Bound | Upper Bound |
|---|-----------|--------|-------|--------------|-------------|-------------|
| 1 | LPL hydrolysis rate | K_ad | LPL lipolysis | 0.00958 | 0.001 | 0.1 |
| 2 | LPL insulin delay | LPL_delay | LPL lipolysis | 155 | 100 | 300 |
| 3 | Fatty acid spillover | D_spill | Spillover | 43.15 | 38 | 50 |
| 4 | GLUT1 transport rate | GLUT_1 | Glucose flux | 0.01809 | 0.001 | 0.1 |
| 5 | GLUT4 transport rate | GLUT_4 | Glucose flux | 0.000552 | 0 | 0.1 |
| 6 | Adipose insulin delay | AT_delay | Glucose flux | 5.384 | 1 | 45 |
| 7 | Glycerol efflux rate | pgly | Glycerol flux | 0.2675 | 0.001 | 1 |
| 8 | Basal lipolysis rate | B_ATL | ATL lipolysis | 0.02781 | 0.001 | 1 |
| 9 | Max suppressible lipolysis | ATL_max | ATL lipolysis | 1.329 | 0.01 | 5 |
| 10 | Lipolysis insulin sensitivity | K_ATL | ATL lipolysis | 1.713 | 0.01 | 5 |
| 11 | Re-esterification rate | re_ester | NEFA flux | 0.001 | 0 | 0.05 |
| 12 | NEFA transport coefficient | pnefa | NEFA flux | 0.0532 | 0 | 1.5 |
| 13 | G-3-P production delay | G3P_delay | NEFA flux | 60 | 1 | 100 |
| 14 | Glucose fraction to G-3-P | glucose_use | NEFA flux | 0.9 | 0 | 1 |
