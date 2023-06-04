## LAB 6 Benchmarks

# lab6 Code:
time python benchmark.py 

../memusg/memusg python ../lab6-spring23/benchmark.py

# scRNAseqFilterQC Code:
time python scRNAseqFilterQC.py counts/ -n 0.01 -t 0.01 -p 0.125 \
    -g GCG TTR IAPP GHRL PPY COL3A1 CPA1 CLPS REG1A CTRB1 CTRB2 PRSS2 CPA2 KRT19 INS SST CELA3A VTCN1

../memusg/memusg python scRNAseqFilterQC.py counts/ -n 0.01 -t 0.01 -p 0.125 \
    -g GCG TTR IAPP GHRL PPY COL3A1 CPA1 CLPS REG1A CTRB1 CTRB2 PRSS2 CPA2 KRT19 INS SST CELA3A VTCN1


