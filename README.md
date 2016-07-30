## Usage

1. Process data into input.Rdata with input.R or similar tool. The output is Rdata file that contains objects with pre-processed input data.

2. Modify parameters in param.R. Use input.R:inspect_dataset() to estimate co-occurrence values.

3. Run algorithm:

```
source('klink-2.R') # compiles C code as well
klink2('input.Rdata')
```

triples object will contain output semantic relations.
