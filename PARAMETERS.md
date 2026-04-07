# Model Parameters

The adipose tissue model uses 14 key parameters. Below is a brief biological description of each parameter and typical units.

| Parameter | Description | Typical Units |
|-----------|-------------|---------------|
| **K_ad** | Rate constant for adipose tissue lipoprotein lipase (LPL) mediated triglyceride hydrolysis. Controls the speed at which LPL clears circulating TG in the presence of insulin. | min⁻¹ |
| **LPL_delay** | Time delay (minutes) between the insulin signal and activation of LPL activity in adipose tissue. | min |
| **D_spill** | Fractional spillover factor representing the proportion of fatty acids released from adipose tissue that appear in the arterial circulation, scaled by insulin concentration. | % |
| **GLUT_1** | Basal glucose transport rate constant via the GLUT1 transporter, independent of insulin. | min⁻¹ |
| **GLUT_4** | Insulin‑dependent glucose transport rate constant via the GLUT4 transporter. | min⁻¹ |
| **AT_delay** | Delay (minutes) for insulin to exert its effect on adipose tissue glucose uptake. | min |
| **pgly** | Rate constant for glycerol production from triglyceride hydrolysis (or other glycerol‑generating pathways). | min⁻¹ |
| **B_ATL** | Basal activity level of adipose triglyceride lipase (ATGL) before insulin stimulation. | min⁻¹ |
| **ATL_max** | Maximum (Vmax) activity of ATGL under saturating substrate conditions. | min⁻¹ |
| **K_ATL** | Michaelis–Menten constant for ATGL; substrate concentration at which ATGL operates at half its maximum rate. | mM |
| **re_ester** | Rate constant for re‑esterification of non‑esterified fatty acids (NEFA) back into triglycerides within adipose tissue. | min⁻¹ |
| **pnefa** | Rate constant governing insulin‑independent NEFA release from adipose tissue. | min⁻¹ |
| **G3P_delay** | Delay (minutes) between glucose uptake and formation of glycerol‑3‑phosphate, which is required for triglyceride synthesis. | min |
| **glucose_use** | Rate constant for utilization of glucose for glycerol‑3‑phosphate production and other metabolic pathways in adipose cells. | min⁻¹ |
