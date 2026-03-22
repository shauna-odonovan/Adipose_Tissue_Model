# Model Parameters

The adipose tissue model uses 14 calibrated parameters. Below are brief biological explanations for each parameter.

| Parameter | Description |
|-----------|-------------|
| K_ad | Maximum rate constant for insulin‑dependent activation of adipose tissue lipoprotein lipase (LPL)–mediated triglyceride hydrolysis. |
| LPL_delay | Time delay (minutes) representing the lag between changes in plasma insulin and the activation of LPL activity in the adipose microvasculature. |
| D_spill | Fractional spillover coefficient describing the proportion of fatty acids released by LPL that spill over into the circulating free fatty acid (NEFA) pool, modulated by insulin. |
| GLUT_1 | Basal glucose uptake rate via insulin‑independent GLUT‑1 transporters in adipocytes. |
| GLUT_4 | Insulin‑stimulated glucose uptake rate via GLUT‑4 transporters in adipocytes. |
| AT_delay | Delay (minutes) capturing the intracellular signaling lag for insulin to affect glucose uptake in adipose tissue. |
| pgly | Basal glycerol release rate constant from triglyceride hydrolysis, independent of insulin. |
| B_ATL | Maximum rate of glycerol‑3‑phosphate conversion to triacylglycerol (TAG) synthesis (ATL pathway). |
| ATL_max | Saturation level (Vmax) for glycerol‑3‑phosphate utilization in TAG synthesis. |
| K_ATL | Michaelis‑Menten constant (Km) for the ATL conversion step. |
| re_ester | Fraction of free fatty acids that are re‑esterified into intracellular triglycerides within the adipocyte. |
| pnefa | Basal non‑esterified fatty acid (NEFA) release rate constant representing insulin‑independent lipolysis. |
| G3P_delay | Delay (minutes) describing the time lag between glucose uptake and the production of glycerol‑3‑phosphate in the adipocyte. |
| glucose_use | Rate constant for intracellular glucose utilization (glycolysis) that supplies glycerol‑3‑phosphate for TAG synthesis. |
