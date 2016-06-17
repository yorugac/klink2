## Usage

1. Process data into input.Rdata with input.R or similar tool.

2. Modify parameters in param.R.

3. Run algorithm:

```
source('klink-2.R')
load('input.Rdata') # objects with pre-processed input data
klink2()
```

triples object will contain output semantic relations.
