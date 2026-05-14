# Model Parameters

This document provides biological explanations for all 14 model parameters used in the Adipose Tissue Model MATLAB code. The model simulates postprandial adipose tissue metabolism, as described in O'Donovan et al., *PLOS Computational Biology* (2019) 15(10): e1007400.

For further mathematical detail on parameter estimation and flux predictions, see `model_description.pdf`.

---

## 1. `K_ad` — LPL-Mediated TG Lipolysis Rate Constant

**Default value:** 0.01 (1/min)

**Biological explanation:** This parameter represents the rate constant for lipoprotein lipase (LPL)-mediated hydrolysis of circulating triglycerides (TG) within the adipose tissue capillary bed. LPL is an enzyme anchored to the luminal surface of endothelial cells that catalyses the breakdown of TG-rich lipoproteins (chylomicrons and VLDL) into free fatty acids and glycerol. The rate of this lipolytic reaction is proportional to `K_ad`, modulated by the local insulin-delayed LPL activity signal. An increase in `K_ad` corresponds to enhanced LPL activity and greater TG extraction from the bloodstream.

---

## 2. `LPL_delay` — Insulin Delay for LPL Activation

**Default value:** 158.13 (min)

**Biological explanation:** This parameter is the time constant of the three-compartment insulin delay model specifically governing LPL activation. After a meal, insulin levels rise in the plasma, but the activation of LPL on the endothelial surface is not instantaneous — it involves transcriptional and translational regulation, intracellular processing, and translocation of the enzyme to the capillary lumen. `LPL_delay` captures this lag, with a larger value indicating a slower insulin-mediated activation of LPL.

---

## 3. `D_spill` — Fractional Spillover Parameter

**Default value:** 41.53

**Biological explanation:** During LPL-mediated TG lipolysis in the capillary bed, not all liberated fatty acids are taken up by the adipocyte; a fraction "spills over" into the venous circulation and contributes to the systemic NEFA pool. `D_spill` controls the partitioning of LPL-released NEFA between adipocyte uptake and spillover. A higher value means a greater proportion of LPL-released fatty acids escape entrapment by the adipocyte, increasing the spillover fraction observed in the venous effluent.

---

## 4. `GLUT_1` — GLUT1 Transporter-Mediated Glucose Uptake Rate

**Default value:** 0.0159 (1/min)

**Biological explanation:** This parameter governs the basal (insulin-independent) uptake of glucose into adipocytes via the GLUT1 transporter. GLUT1 is constitutively present on the adipocyte plasma membrane and facilitates passive glucose diffusion along its concentration gradient. `GLUT_1` represents the rate at which glucose is transported into the cell independent of insulin stimulation, providing a basal glucose supply for glycolysis and glycerol-3-phosphate (G-3-P) production.

---

## 5. `GLUT_4` — GLUT4 Transporter-Mediated Insulin-Stimulated Glucose Uptake Rate

**Default value:** 0.00082 (1/min per μU/mL insulin)

**Biological explanation:** This parameter scales the insulin-dependent component of glucose uptake mediated by the GLUT4 transporter. Upon insulin stimulation, GLUT4 vesicles are translocated from intracellular storage compartments to the plasma membrane, dramatically increasing the cell's glucose permeability. `GLUT_4` represents the sensitivity of this translocation process to the insulin signal arriving at the adipocyte (after the AT_delay), thus controlling the postprandial augmentation of glucose uptake.

---

## 6. `AT_delay` — Adipose Tissue Insulin Delay

**Default value:** 13.85 (min)

**Biological explanation:** This parameter is the time constant of the three-compartment insulin delay model governing insulin action on adipose tissue glucose uptake and intracellular lipolysis. After plasma insulin rises, there is a physiological lag before insulin exerts its effects on the adipocyte — due to receptor binding, intracellular signalling cascades (e.g., PI3K/Akt pathway), and GLUT4 translocation for glucose uptake or hormone-sensitive lipase inhibition for lipolysis suppression. `AT_delay` captures this transduction lag for adipocyte-specific insulin actions.

---

## 7. `pgly` — Glycerol Permeability / Transport Rate

**Default value:** 0.2705 (1/min)

**Biological explanation:** This parameter represents the rate constant for glycerol transport across the adipocyte plasma membrane. Adipocytes express aquaglyceroporins (primarily AQP7) that facilitate glycerol efflux. During lipolysis, intracellular glycerol is generated from TG breakdown and can either be re-utilised for re-esterification or released into the circulation. `pgly` governs the net flux of glycerol out of the adipocyte, driven by the concentration gradient between intracellular and arterial glycerol levels.

---

## 8. `B_ATL` — Basal Adipose Tissue Lipolysis Rate

**Default value:** 0.0013 (mmol/L/min)

**Biological explanation:** This parameter represents the basal (insulin-independent) rate of intracellular triglyceride lipolysis in adipose tissue. Even in the absence of lipolytic stimulation (e.g., by catecholamines), there is a constitutive rate of TG breakdown driven by the basal activity of adipose triglyceride lipase (ATGL) and hormone-sensitive lipase (HSL). `B_ATL` sets this minimum lipolytic flux, which contributes to the background NEFA and glycerol release from adipocytes.

---

## 9. `ATL_max` — Maximum Adipose Tissue Lipolysis Rate

**Default value:** 1.4373 (mmol/L/min)

**Biological explanation:** This parameter represents the maximum achievable rate of intracellular TG lipolysis in adipose tissue under strong lipolytic stimulation. In the model, the insulin-suppressible component of lipolysis is represented as `ATL_max / (1 + I_AT / K_ATL)`, where `I_AT` is the delayed intracellular insulin signal. When insulin is absent, lipolysis approaches `B_ATL + ATL_max`; when insulin is high, the lipolytic rate is suppressed toward `B_ATL`. `ATL_max` thus captures the dynamic range of lipolytic regulation.

---

## 10. `K_ATL` — ATL Insulin Sensitivity (Half-Saturation Constant)

**Default value:** 1.49185 (μU/mL equivalent)

**Biological explanation:** This parameter is the half-saturation constant in the Michaelis-Menten-type equation governing insulin suppression of adipose tissue lipolysis. It determines the insulin concentration at which the insulin-suppressible component of lipolysis is reduced by 50%. A lower `K_ATL` means the adipocyte is more sensitive to insulin's anti-lipolytic effect (half-maximal suppression at a lower insulin level), whereas a higher `K_ATL` indicates insulin resistance of the lipolytic pathway.

---

## 11. `re_ester` — Re-Esterification Rate of Adipose NEFA

**Default value:** 0.000295 (1/min per mmol/L G-3-P)

**Biological explanation:** This parameter governs the rate at which intracellular NEFA are re-esterified with glycerol-3-phosphate (G-3-P) to form triglycerides. Re-esterification is a critical cycle in adipose tissue metabolism: NEFA released by lipolysis can be re-incorporated into TG, and this process is highly dependent on the availability of G-3-P (itself derived from glucose metabolism). `re_ester` scales the G-3-P-dependent re-esterification flux, linking glucose availability to fatty acid trapping within the adipocyte.

---

## 12. `pnefa` — NEFA Transport Permeability / Rate

**Default value:** 0.0478 (1/min)

**Biological explanation:** This parameter represents the rate constant for NEFA transport across the adipocyte plasma membrane. NEFA enter and exit adipocytes via both passive diffusion and protein-mediated transport (e.g., CD36, FATP family). `pnefa` governs the net efflux of NEFA from the intracellular pool into the venous circulation, driven by the NEFA concentration gradient. It also affects the rate at which arterial NEFA are taken up by the adipocyte. A higher `pnefa` facilitates greater NEFA exchange between the cell and the bloodstream.

---

## 13. `G3P_delay` — G-3-P Delay Time Constant

**Default value:** 31.74 (min)

**Biological explanation:** This parameter represents the time constant for the two-compartment delay model governing the production of glycerol-3-phosphate (G-3-P) from glucose. After glucose enters the adipocyte and is metabolised via glycolysis to produce G-3-P, there is a metabolic lag due to the time required for glycolytic processing and enzymatic conversion. `G3P_delay` captures this lag, ensuring that the G-3-P concentration (which drives re-esterification) responds to changes in glucose uptake with a physiologically realistic time course.

---

## 14. `glucose_use` — Proportion of Glucose Used for G-3-P Production

**Default value:** 0.8803

**Biological explanation:** This parameter represents the fraction of glucose taken up by the adipocyte that is directed toward the production of glycerol-3-phosphate (G-3-P) via glycolysis. Not all internalised glucose is converted to G-3-P; a portion is oxidised via the TCA cycle, converted to lactate, or used for other metabolic purposes. `glucose_use` determines the efficiency of the glyceroneogenic pathway, directly affecting the availability of G-3-P for NEFA re-esterification. A value close to 1 indicates that nearly all glucose is used for G-3-P generation, maximising the capacity for fatty acid trapping.

---

## Parameter Summary Table

| # | Parameter | Default Value | Units | Category |
|---|-----------|---------------|-------|----------|
| 1 | `K_ad` | 0.01 | 1/min | LPL lipolysis |
| 2 | `LPL_delay` | 158.13 | min | LPL lipolysis |
| 3 | `D_spill` | 41.53 | — | Spillover |
| 4 | `GLUT_1` | 0.0159 | 1/min | Glucose flux |
| 5 | `GLUT_4` | 0.00082 | 1/min per μU/mL | Glucose flux |
| 6 | `AT_delay` | 13.85 | min | Glucose flux |
| 7 | `pgly` | 0.2705 | 1/min | Glycerol flux |
| 8 | `B_ATL` | 0.0013 | mmol/L/min | Glycerol flux |
| 9 | `ATL_max` | 1.4373 | mmol/L/min | Glycerol flux |
| 10 | `K_ATL` | 1.49185 | μU/mL equiv. | Glycerol flux |
| 11 | `re_ester` | 0.000295 | 1/min per mmol/L G-3-P | NEFA |
| 12 | `pnefa` | 0.0478 | 1/min | NEFA |
| 13 | `G3P_delay` | 31.74 | min | NEFA |
| 14 | `glucose_use` | 0.8803 | — | NEFA |

## Parameter Estimation Bounds

| Parameter | Lower Bound | Upper Bound |
|-----------|-------------|-------------|
| `K_ad` | 0.001 | 0.1 |
| `LPL_delay` | 100 | 300 |
| `D_spill` | 38 | 50 |
| `GLUT_1` | 0.001 | 0.1 |
| `GLUT_4` | 0 | 0.1 |
| `AT_delay` | 1 | 45 |
| `pgly` | 0.001 | 1 |
| `B_ATL` | 0.001 | 1 |
| `ATL_max` | 0.01 | 5 |
| `K_ATL` | 0.01 | 5 |
| `re_ester` | 0 | 0.05 |
| `pnefa` | 0 | 1.5 |
| `G3P_delay` | 1 | 100 |
| `glucose_use` | 0 | 1 |
