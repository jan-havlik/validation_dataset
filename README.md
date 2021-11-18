# Repository for data validaton between QmRLFS finder and R-loop tracker tools

Each folder contains source data for the validation:
- `gene_data` folder contains raw sequences in `fasta` format
- `experimental_data` folder contains DRIPc sequenced R-loops in `fasta` format
- `rlooptracker_data` folder contains R-loop tracker analyses in `bedgraph` format

The script provides the following metrics for validation:

- `Accuracy`
- `Sensitivity`
- `Specificity`
- `Precision`
- `Matthews Correlation Coefficient`

Given dataset should provide following stats:


# Validating R-loop forming signal detecton with experimental data

## Plus strand

| Gene | RLFS predicted by DRIPc sequencing | RLFS predicted by R-loop tracker |
|:---:|:---:|:---:|
| Immunoglobulin | - | + |
| MYC | + | + |
| RHOH | - | + |
| ACTB | + | + |
| FMR1 | + | + |
| SNRPN | + | + |
| HK2 | + | + |
| CIRH1A | + | + |
| APOE | + | + |
| FHIT | + | + |
| PPM1D | + | + |
| TP53 | - | + |
| JTB | - | - |
| PBX1 | + | + |

```
Accuracy: 78.57 %
Specificity: 25.00 %
Sensitivity: 100.00 %
Precision: 76.92 %
Matthews Correlation Coefficient: 0.44
```

The detection algorithm produces more False positives compared to experimental methods which can be seen in Specificity, but overall performance and accuracy is adequate.

## Minus strand

| Gene | RLFS predicted by DRIPc sequencing | RLFS predicted by R-loop tracker |
|:---:|:---:|:---:|
| Immunoglobulin | - | - |
| MYC | - | + |
| RHOH | + | + |
| ACTB | + | + |
| FMR1 | - | - |
| SNRPN | + | - |
| HK2 | + | + |
| CIRH1A | + | + |
| APOE | + | + |
| FHIT | + | + |
| PPM1D | - | + |
| TP53 | + | + |
| JTB | + | - |
| PBX1 | - | + |

```
Accuracy: 64.29 %
Specificity: 40.00 %
Sensitivity: 77.78 %
Precision: 70.00 %
Matthews Correlation Coefficient: 0.19
```

Results for negative strand are bit worse than for the positive strand and interestingly the algorithmic detection produces more False negatives.
