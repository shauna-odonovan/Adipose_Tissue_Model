# Model Parameters

The adipose tissue model uses 14 fitted parameters (p_opt(1)–p_opt(14)) that describe key metabolic processes. The biological meaning of each parameter is summarized below.

| Index | Parameter (MATLAB variable) | Biological description |
|-------|-----------------------------|------------------------|
| 1 | **K_ad** | LPL (lipoprotein lipase) activity coefficient that scales insulin‑stimulated triglyceride hydrolysis in adipose tissue. It reflects the sensitivity of LPL to insulin. |
| 2 | **LPL_delay** | Time delay (minutes) for insulin to activate LPL activity. |
| 3 | **D_spill** | Fractional spill‑over coefficient describing the proportion of fatty acids released from adipose tissue that appear in the systemic circulation as free NEFA, modulated by insulin. |
| 4 | **GLUT_1** | Basal (insulin‑independent) glucose uptake rate constant via the GLUT1 transporter in adipocytes. |
| 5 | **GLUT_4** | Insulin‑dependent glucose uptake rate constant via the GLUT4 transporter in adipocytes. |
| 6 | **AT_delay** | Delay (minutes) for insulin to affect adipose‑tissue glucose uptake and associated metabolic pathways. |
| 7 | **pgly** | Basal glycerol production rate from lipolysis, independent of insulin signalling. |
| 8 | **B_ATL** | Basal term governing insulin‑stimulated glycerol release (A_TL) from adipose tissue. |
| 9 | **ATL_max** | Maximal rate of insulin‑stimulated glycerol release (A_TL) from adipocytes. |
|10 | **K_ATL** | Half‑maximal insulin concentration for the glycerol release response (A_TL). |
|11 | **re_ester** | Rate constant for re‑esterification of NEFA back into triglycerides within the adipocyte. |
|12 | **pnefa** | Scaling factor for the net NEFA release component that is not re‑esterified; determines the magnitude of NEFA appearance in plasma. |
|13 | **G3P_delay** | Delay (minutes) for glucose‑derived glycerol‑3‑phosphate to become available for NEFA synthesis. |
|14 | **glucose_use** | Fraction of intracellular glucose directed toward glycerol‑3‑phosphate formation for NEFA re‑esterification. |

These parameters are estimated by fitting the model to arteriovenous tracer data (see `Model_description.pdf`).
