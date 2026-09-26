# Table 13. External public 2bRAD-M route-B classification

All runs used route B, BcgI, mismatch 0, no residual GC correction, and LOOCV. Growth PTR features used robust estimates only; missingness was represented explicitly. Bootstrap intervals resample samples (5,000 times). Taxonomic support is the route-B coverage/support layer, not an independently published abundance table.

| Dataset | Feature arm | AUROC | AUROC bootstrap 95% CI | AUPRC |
|---|---|---:|---:|---:|
| Ovarian cancer vs benign | Taxonomic support | 0.780 | 0.515–1.000 | 0.730 |
| Ovarian cancer vs benign | Growth PTR | **0.930** | 0.786–1.000 | 0.929 |
| Ovarian cancer vs benign | Combined | 0.810 | 0.560–1.000 | 0.739 |
| Bladder NMIBC vs MIBC | Taxonomic support | **0.990** | 0.946–1.000 | 0.982 |
| Bladder NMIBC vs MIBC | Growth PTR | 0.957 | 0.857–1.000 | 0.913 |
| Bladder NMIBC vs MIBC | Combined | 0.981 | 0.914–1.000 | 0.968 |
| Maternal faeces vs meconium | Taxonomic support | 0.617 | 0.476–0.750 | 0.665 |
| Maternal faeces vs meconium | Growth PTR | 0.538 | 0.393–0.683 | 0.581 |
| Maternal faeces vs meconium | Combined | 0.594 | 0.454–0.734 | 0.649 |

Point estimates are exploratory: no cohort had a paired growth-versus-support improvement that passed permutation testing at 0.05 after accounting for the direction specified in the text.
