# Model Parameters

This document provides biological explanations for all 14 model parameters used in the adipose tissue metabolism model (O'Donovan et al., *PLOS Comp. Bio.*, 2019). These parameters are estimated in `Estimating_Model_Parameters_Script.m` and used throughout the ODE solver functions (`Glucose_G3P.m`, `Glucose_G3P_NEFA_ODE.m`, `Glycerol_ODE.m`).

---

## LPL-Mediated Lipolysis Group

### 1. K_ad (`p_opt(1)`)
**Rate constant for LPL-mediated lipolysis of plasma triglycerides.**

Lipoprotein lipase (LPL) is an enzyme anchored to the capillary endothelium that hydrolyzes circulating triglyceride-rich lipoproteins (chylomicrons and VLDL) into free fatty acids and glycerol. K_ad quantifies the efficiency with which LPL clears arterial triglycerides in adipose tissue. A higher K_ad reflects greater LPL activity, leading to more rapid uptake of circulating fatty acids into adipocytes for storage. In the model, arterial TG uptake is proportional to `K_ad * TG_art * I_LPL`, where `I_LPL` is the insulin-delayed signal activating LPL.

- **Typical range:** 0.001–0.1
- **Fitted value:** ~0.0096

### 2. LPL_delay (`p_opt(2)`)
**Time delay (minutes) for insulin activation of lipoprotein lipase.**

Insulin stimulates LPL synthesis and translocation to the capillary lumen, but this process is not instantaneous. LPL_delay represents the lag between rising plasma insulin levels and the corresponding increase in LPL activity within adipose tissue capillaries. This delay captures the time required for insulin signaling, LPL gene transcription, translation, and secretion to the endothelial surface.

- **Typical range:** 100–300 minutes
- **Fitted value:** ~155 minutes

---

## Spillover Group

### 3. D_spill (`p_opt(3)`)
**Fractional spillover parameter for fatty acids escaping adipose tissue uptake.**

Not all fatty acids released by LPL-mediated hydrolysis of triglycerides are taken up by the underlying adipocytes. A fraction "spills over" into the venous circulation and becomes available to other tissues. D_spill determines the proportion of LPL-released fatty acids that are not captured by adipocytes, modulated by the insulin ratio between systemic plasma and the local LPL compartment. Spillover is calculated as `(D_spill * (I_b / I_PL)) / 100`.

- **Typical range:** 38–50
- **Fitted value:** ~43.15

---

## Glucose Flux Group

### 4. GLUT_1 (`p_opt(4)`)
**Rate constant for insulin-independent (basal) glucose uptake via GLUT1 transporters.**

GLUT1 transporters are constitutively present on the adipocyte membrane and mediate basal glucose uptake that is not regulated by insulin. This parameter represents the rate at which glucose enters adipocytes through GLUT1, driven by the arterial glucose concentration gradient. It ensures a baseline level of glucose supply for glycolysis and glycerol-3-phosphate (G-3-P) production even in the fasted state.

- **Typical range:** 0.001–0.1
- **Fitted value:** ~0.018

### 5. GLUT_4 (`p_opt(5)`)
**Rate constant for insulin-dependent glucose uptake via GLUT4 transporters.**

GLUT4 transporters are stored in intracellular vesicles and translocate to the adipocyte membrane in response to insulin signaling. This parameter quantifies the insulin-stimulated glucose uptake rate, proportional to both the local insulin concentration in adipose tissue (`I_AT`) and the arterial glucose level. After a meal, rising insulin causes GLUT4 translocation, dramatically increasing glucose flux into adipocytes. Total glucose uptake: `GLUT_1 * G_art + GLUT_4 * I_AT * G_art`.

- **Typical range:** 0–0.1
- **Fitted value:** ~0.00055

### 6. AT_delay (`p_opt(6)`)
**Time delay (minutes) for insulin action on adipose tissue glucose uptake and lipolysis.**

After insulin enters the interstitial space of adipose tissue, there is a physiological delay before it exerts its metabolic effects (GLUT4 translocation, suppression of lipolysis). AT_delay represents this lag, which is modeled using a three-compartment delay system. The delay accounts for insulin receptor binding, signal transduction cascade, and the time needed for cellular responses (e.g., vesicle trafficking for GLUT4).

- **Typical range:** 1–45 minutes
- **Fitted value:** ~5.4 minutes

---

## Glycerol Flux / Adipose Tissue Lipolysis Group

### 7. pgly (`p_opt(7)`)
**Glycerol transport/diffusion parameter across the adipocyte membrane.**

Glycerol is produced intracellularly by the hydrolysis of stored triglycerides (lipolysis) and can exit adipocytes through aquaporin-7 channels into the interstitial space, from where it enters the circulation. pgly governs the rate of glycerol efflux from adipose tissue, influencing the arteriovenous glycerol difference. It appears in the glycerol mass balance as a proportional transport term.

- **Typical range:** 0.001–1
- **Fitted value:** ~0.267

### 8. B_ATL (`p_opt(8)`)
**Basal (insulin-independent) adipose tissue lipolysis rate.**

Even in the presence of insulin, there is a basal rate of triglyceride hydrolysis within adipocytes, mediated by constitutive activity of adipose triglyceride lipase (ATGL) and hormone-sensitive lipase (HSL). B_ATL represents this minimum lipolytic flux, ensuring a constant release of glycerol and fatty acids regardless of hormonal status. This basal lipolysis is important for providing fuel during fasting.

- **Typical range:** 0.001–1
- **Fitted value:** ~0.028

### 9. ATL_max (`p_opt(9)`)
**Maximum adipose tissue lipolysis rate under insulin-suppressed conditions.**

When insulin levels are very low (e.g., during prolonged fasting), lipolysis is maximally activated by catecholamines and other lipolytic hormones. ATL_max represents the upper limit of the lipolytic rate. The actual lipolytic rate is modeled as `B_ATL + ATL_max / (1 + I_AT / K_ATL)`, where insulin suppresses lipolysis in a dose-dependent manner through the insulin-sensitivity constant K_ATL.

- **Typical range:** 0.01–5
- **Fitted value:** ~1.33

### 10. K_ATL (`p_opt(10)`)
**Insulin sensitivity constant for suppression of adipose tissue lipolysis.**

K_ATL determines the insulin concentration at which lipolysis is half-maximally suppressed. It appears in the denominator of the lipolysis inhibition term `ATL_max / (1 + I_AT / K_ATL)`. A lower K_ATL means adipose tissue is more sensitive to insulin's anti-lipolytic effect (lipolysis is suppressed at lower insulin concentrations). This parameter captures inter-individual variation in insulin sensitivity at the level of the adipocyte.

- **Typical range:** 0.01–5
- **Fitted value:** ~1.71

---

## NEFA (Non-Esterified Fatty Acid) Group

### 11. re_ester (`p_opt(11)`)
**Rate constant for insulin-stimulated re-esterification of adipose NEFA.**

Re-esterification is the process by which free fatty acids (released by lipolysis) are re-incorporated into triglycerides within adipocytes. This process requires glycerol-3-phosphate (G-3-P), which is derived from glucose metabolism. Insulin promotes re-esterification by enhancing glucose uptake (providing G-3-P) and by stimulating the esterification enzymes. The re-esterification rate is modeled as `re_ester * I_AT * NEFA_at * G_3_P`, linking it to both insulin and G-3-P availability.

- **Typical range:** 0–0.05
- **Fitted value:** ~0.001

### 12. pnefa (`p_opt(12)`)
**NEFA transport/diffusion parameter across the adipocyte membrane.**

Non-esterified fatty acids released by lipolysis can either be re-esterified or released into the circulation. pnefa governs the rate of NEFA efflux from the adipose tissue interstitial space into venous blood, influencing the arteriovenous NEFA difference. A higher pnefa value means fatty acids more readily leave adipose tissue, contributing to higher circulating NEFA levels.

- **Typical range:** 0–1.5
- **Fitted value:** ~0.053

### 13. G3P_delay (`p_opt(13)`)
**Time delay (minutes) for glycerol-3-phosphate (G-3-P) production from glucose.**

After glucose enters the adipocyte, it must be metabolized through glycolysis to produce glycerol-3-phosphate, the backbone required for triglyceride synthesis. G3P_delay represents the lag between glucose uptake and the availability of G-3-P for re-esterification. This delay is modeled using a three-compartment system and accounts for the time required for glycolytic processing and G-3-P accumulation.

- **Typical range:** 1–100 minutes
- **Fitted value:** ~60 minutes

### 14. glucose_use (`p_opt(14)`)
**Proportion of intracellular glucose directed toward G-3-P production.**

Not all glucose taken up by adipocytes is used for glycerol-3-phosphate synthesis; some is fully oxidized for energy or used in other metabolic pathways. glucose_use represents the fraction of glucose flux that is diverted to G-3-P production, which is essential for triglyceride re-esterification. A value near 1 indicates nearly all glucose is used for G-3-P, while lower values indicate greater glucose diversion to alternative pathways (e.g., oxidation, lactate production).

- **Typical range:** 0–1
- **Fitted value:** ~0.9

---

## Parameter Estimation

All 14 parameters are estimated simultaneously by fitting the model predictions to arteriovenous difference measurements across adipose tissue, using stable isotope tracer data from human subjects. The estimation is performed in `Estimating_Model_Parameters_Script.m` using a constrained optimization approach (`fmincon`) that minimizes the sum of squared errors between model predictions and observed data. Parameter bounds are set based on physiological plausibility.

For details on the fitting procedure, see `Parameter_Fit_Adipose_Model.m` and the accompanying publication.
