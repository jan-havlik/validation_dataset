# Repository for data validaton between QmRLFS finder and R-loop tracker tools

Each folder contains source data for the validation:
- `source_data` folder contains sequences in `fasta` format
- `qmrlfs_data` folder contains QmRLFS finder analyses in `bedgraph` format
- `rlooptracker_data` folder contains R-loop tracker analyses in `bedgraph` format

There is also script included to obtain validation statistics for given datasets. Function `download_seq` was used for downloading sequences from the genome browser via API and the function `compare_tools` runs the validation. The output metrics are as follows:

- `Accuracy`
- `Sensitivity`
- `Specificity`
- `Precision`
- `Matthews Correlation Coefficient`

Given dataset should provide following stats:
```
Accuracy: 90.91 %
Sensitivity: 95.50 %
Specificity: 40.00 %
Precision: 94.64 %
Matthews Correlation Coefficient: 0.37
```

