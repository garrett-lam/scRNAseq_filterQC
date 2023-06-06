# Benchmark: scRNAseqFilterQC (Lab 6 Custom) vs Lab 6 Reference

To benchmark our tool (scRNAseqFilterQC) against the Lab 6 code: we used the `time` and `memusg` packages
- The `memusg` was downloaded from https://github.com/jhclark/memusg 

## Basic Usage
To run the benchmarking code:
`./benchmark.sh`


## Summary of Results 


|                      | Lab 6 Reference | scRNAFilterQC (Lab 6 Custom) |
|----------------------|-----------------|----------------------------|
| real (Runtime)                 | 2m0.955s        | 3m37.270s                  |
| user (Runtime)            | 4m10.671s       | 8m8.331s                   |
| sys  (Runtime)                | 2m43.736s       | 0m37.429s                  |
| Memory Usage         | 5.943 GB        | 5.397 GB                   |
