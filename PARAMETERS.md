# Model Parameters

This document provides biological explanations for all 14 model parameters used in the adipose tissue metabolic model implemented in the MATLAB code. The parameters are estimated in `Parameter_Fit_Adipose_Model.m` and referenced as `p_opt(1)` through `p_opt(14)`.

## 1. K_ad — LPL Lipolysis Rate Constant

- **Symbol in code:** `p_opt(1)`
- **Typical value:** 0.0096
- **Biological explanation:** This parameter represents the rate constant for lipoprotein lipase (LPL)-mediated hydrolysis of circulating triglycerides (TG) in adipose tissue. LPL is an enzyme bound to the capillary endothelium that hydrolyzes chylomicron and VLDL triglycerides into free fatty acids and glycerol, enabling their uptake by adipocytes. `K_ad` scales the insulin-activated LPL activity against the arterial triglyceride concentration, determining the net TG extraction flux: `D_TG = -K_ad * I_LPL * TG_art`.

## 2. LPL_delay — Insulin Delay for LPL Activation

- **Symbol in code:** `p_opt(2)`
- **Typical value:** 155 (minutes)
- **Biological explanation:** This parameter represents the time delay (in minutes) between changes in plasma insulin concentration and the corresponding activation of lipoprotein lipase (LPL) in adipose tissue capillaries. LPL activation by insulin is not instantaneous; it involves transcriptional regulation, enzyme synthesis, and translocation to the capillary lumen. The model implements this as a three-fold insulin delay chain in `insulin_delay.m`, where `tau = LPL_delay` governs the lag between plasma insulin (`I_pl`) and the delayed insulin signal (`I_LPL`) that modulates LPL activity.

## 3. D_spill — Fractional Spillover

- **Symbol in code:** `p_opt(3)`
- **Typical value:** 43.15
- **Biological explanation:** This parameter quantifies the fraction of fatty acids released by LPL-mediated lipolysis of triglycerides that "spill over" into the plasma rather than being taken up by the adipose tissue. In postprandial metabolism, not all fatty acids liberated from chylomicrons and VLDL by LPL are captured by adipocytes — a proportion escapes into the venous circulation. The spillover fraction depends on the ratio of interstitial to plasma insulin (`I_b/I_PL`), reflecting insulin's role in promoting fatty acid trapping: `spill = (D_spill * (I_b / I_PL)) / 100`.

## 4. GLUT_1 — Insulin-Independent Glucose Transport

- **Symbol in code:** `p_opt(4)`
- **Typical value:** 0.0181
- **Biological explanation:** This parameter represents the rate of basal (insulin-independent) glucose uptake into adipose tissue mediated by GLUT1 transporters. GLUT1 is constitutively present on the adipocyte membrane and facilitates passive glucose diffusion down its concentration gradient. This basal uptake operates regardless of insulin levels and is combined with insulin-dependent uptake to determine total glucose extraction: `D_G = -GLUT_1 * G_art - GLUT_4 * I_AT * G_art`.

## 5. GLUT_4 — Insulin-Dependent Glucose Transport

- **Symbol in code:** `p_opt(5)`
- **Typical value:** 0.000552
- **Biological explanation:** This parameter represents the rate of insulin-stimulated glucose uptake into adipose tissue mediated by GLUT4 transporters. In the basal state, GLUT4 is sequestered in intracellular vesicles. Upon insulin stimulation, these vesicles translocate to the cell membrane, increasing glucose transport capacity. `GLUT_4` scales the product of delayed adipose tissue insulin (`I_AT`) and arterial glucose concentration to determine insulin-stimulated glucose uptake. This is the primary mechanism by which insulin promotes postprandial glucose disposal in adipose tissue.

## 6. AT_delay — Adipose Tissue Insulin Delay

- **Symbol in code:** `p_opt(6)`
- **Typical value:** 5.38 (minutes)
- **Biological explanation:** This parameter represents the time delay (in minutes) between changes in plasma insulin concentration and the corresponding insulin action within the adipose tissue interstitial compartment. Insulin must traverse the capillary endothelium and diffuse through the interstitial space before reaching adipocyte insulin receptors. The model implements this as a three-fold delay chain in `insulin_delay.m`, where `tau = AT_delay` governs the lag between plasma insulin (`I_pl`) and the interstitial insulin signal (`I_AT`) that modulates GLUT4 translocation, lipolysis suppression, and re-esterification.

## 7. pgly — Glycerol Transport Rate

- **Symbol in code:** `p_opt(7)`
- **Typical value:** 0.267
- **Biological explanation:** This parameter represents the rate of glycerol exchange between adipose tissue and plasma. Glycerol is produced intracellularly when triglycerides are hydrolyzed by adipose tissue lipolysis (via ATGL and HSL). Since adipocytes lack glycerol kinase, glycerol cannot be re-esterified within the same cell and must diffuse out. `pgly` governs the arteriovenous glycerol flux across the adipose tissue bed, accounting for both the passive diffusion of glycerol out of adipocytes and the venous glycerol concentration: `d(gly_at) = -pgly * (gly_at - GLY_art + K_ad * TG_art * I_LPL) + B_ATL + ATL_max / (1 + I_AT / K_ATL)`.

## 8. B_ATL — Basal Adipose Tissue Lipolysis Rate

- **Symbol in code:** `p_opt(8)`
- **Typical value:** 0.0278
- **Biological explanation:** This parameter represents the basal (insulin-independent) rate of intracellular triglyceride lipolysis within adipocytes. Even in the absence of hormonal stimulation, adipose triglyceride lipase (ATGL) and hormone-sensitive lipase (HSL) maintain a low rate of triglyceride hydrolysis, releasing glycerol and fatty acids. `B_ATL` is a constant component of lipolysis that persists regardless of insulin levels, contributing to the basal glycerol and NEFA release from adipose tissue.

## 9. ATL_max — Maximum Insulin-Suppressible Lipolysis Rate

- **Symbol in code:** `p_opt(9)`
- **Typical value:** 1.33
- **Biological explanation:** This parameter represents the maximum rate of insulin-suppressible lipolysis in adipose tissue. In the fasted state (low insulin), intracellular lipolysis is active and driven by ATGL and HSL. As insulin levels rise postprandially, insulin signaling through the cAMP/PKA pathway suppresses lipolysis. `ATL_max` quantifies the insulin-sensitive component of lipolysis that can be fully suppressed at high insulin concentrations, contributing to the total lipolysis equation: `ATL_max / (1 + I_AT / K_ATL)`.

## 10. K_ATL — Insulin Sensitivity of Adipose Tissue Lipolysis

- **Symbol in code:** `p_opt(10)`
- **Typical value:** 1.71
- **Biological explanation:** This parameter represents the insulin concentration at which insulin-suppressible lipolysis is reduced to half of its maximum rate (analogous to a Michaelis-Menten constant). It determines how sensitive adipose tissue lipolysis is to insulin-mediated suppression. A lower `K_ATL` means lipolysis is more sensitive to insulin (suppressed at lower insulin concentrations), while a higher value means greater insulin concentrations are needed for equivalent suppression. The insulin-suppressible lipolysis follows the equation: `ATL_max / (1 + I_AT / K_ATL)`.

## 11. re_ester — Re-esterification Rate

- **Symbol in code:** `p_opt(11)`
- **Typical value:** 0.001
- **Biological explanation:** This parameter represents the rate at which free fatty acids (NEFA) within adipocytes are re-esterified back into triglycerides. Re-esterification is a critical metabolic cycle in adipose tissue where fatty acids are recombined with glycerol-3-phosphate (G-3-P) to form triglycerides. This process is insulin-stimulated and requires G-3-P, which is primarily derived from glucose metabolism. The re-esterification flux is modeled as: `re_ester * I_AT * NEFA_at * G_3_P`, reflecting its dependence on interstitial insulin, intracellular NEFA availability, and G-3-P supply.

## 12. pnefa — NEFA Transport Rate

- **Symbol in code:** `p_opt(12)`
- **Typical value:** 0.0532
- **Biological explanation:** This parameter represents the rate of NEFA (non-esterified fatty acid) exchange between adipose tissue and plasma. NEFA are transported across the adipocyte membrane and the capillary endothelium via both passive diffusion and protein-mediated transport (e.g., CD36, FABPpm). `pnefa` governs the net flux of NEFA between the intracellular and plasma compartments, contributing to the arteriovenous NEFA difference. It determines how rapidly NEFA released by lipolysis enter the circulation and how quickly plasma NEFA can be taken up by adipocytes.

## 13. G3P_delay — G-3-P Production Delay

- **Symbol in code:** `p_opt(13)`
- **Typical value:** 60 (minutes)
- **Biological explanation:** This parameter represents the time delay (in minutes) between glucose uptake by adipocytes and the appearance of glycerol-3-phosphate (G-3-P) from glucose metabolism. After glucose enters the adipocyte via GLUT transporters, it must undergo glycolysis to produce dihydroxyacetone phosphate (DHAP), which is then reduced to G-3-P. This metabolic processing, including insulin-activated glycolytic enzyme regulation, introduces a time lag. The model implements this as a two-compartment delay chain in `Glucose_G3P.m`, where `tau = G3P_delay` governs the lag between glucose uptake and G-3-P availability for re-esterification.

## 14. glucose_use — Proportion of Glucose Converted to G-3-P

- **Symbol in code:** `p_opt(14)`
- **Typical value:** 0.9
- **Biological explanation:** This parameter represents the fraction of glucose taken up by adipose tissue that is metabolized to glycerol-3-phosphate (G-3-P). Not all glucose taken up by adipocytes is directed toward G-3-P production — some may be oxidized via the TCA cycle, converted to lactate, or used for other biosynthetic pathways. `glucose_use` (referred to as `usage` in the code) scales the total glucose uptake to determine the G-3-P production flux: `usage * (GLUT_1 * G_art + GLUT_4 * I_AT * G_art)`. A value of 0.9 indicates that approximately 90% of glucose uptake is directed toward G-3-P generation, which is essential for triglyceride re-esterification.
