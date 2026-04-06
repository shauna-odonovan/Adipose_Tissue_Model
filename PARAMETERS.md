# Model Parameters

The adipose tissue model uses 14 physiological parameters that describe lipid and glucose metabolism in adipocytes. Below is a brief biological explanation of each parameter.

| # | Symbol in code | Biological meaning |
|---|----------------|--------------------|
| 1 | `K_ad` | Lipoprotein lipase (LPL) catalytic constant. Controls the rate at which arterial triglyceride is hydrolyzed to release fatty acids for uptake by adipose tissue. |
| 2 | `LPL_delay` | Time delay (minutes) representing the latency between LPL activation and observable lipolysis in the tissue. |
| 3 | `D_spill` | Fractional fatty‑acid spillover from peripheral tissues. Specifies the proportion of released NEFA that bypasses adipose uptake and appears in systemic circulation. |
| 4 | `GLUT_1` | Insulin‑independent glucose transport rate (GLUT1) into adipocytes. Governs basal glucose uptake from blood. |
| 5 | `GLUT_4` | Insulin‑stimulated glucose transport rate (GLUT4) into adipocytes. Modulated by insulin concentration in the adipose compartment. |
| 6 | `AT_delay` | Delay (minutes) for insulin signaling to affect adipocyte glucose uptake. |
| 7 | `pgly` | Basal glycerol production rate from triglyceride lipolysis. |
| 8 | `B_ATL` | Basal glycerol clearance (or transport) rate from adipose tissue. |
| 9 | `ATL_max` | Maximum (saturable) glycerol transport capacity in adipocytes. |
|10 | `K_ATL` | Affinity constant for glycerol transport, controlling how quickly transport approaches `ATL_max`. |
|11 | `re_ester` | Rate constant for re‑esterification of non‑esterified fatty acids (NEFA) back into triglyceride within adipocytes. |
|12 | `pnefa` | Basal NEFA production rate from lipolysis, independent of insulin signaling. |
|13 | `G3P_delay` | Delay (minutes) for conversion of glucose to glycerol‑3‑phosphate, the glycerol backbone for triglyceride synthesis. |
|14 | `glucose_use` | Overall glucose consumption rate that fuels glycerol‑3‑phosphate production and re‑esterification processes. |
