# Description

Run the following steps for the atlas integration:
- Highly variable genes search
- Integration with scVI for 36 different hyper-parameter sets
- scIB to compare those 36 models

# Workflow

Create all necessary directories
```
mkdir -p logs results data
```

Run highly variable genes search. You will need `src`, `config` and `data` packages installed and accessible to python script
```
bsub < ./scripts/hvg_script.bsub
```

Run scvi hyper-parameter search
```
bsub < ./scripts/run_integration.bsub
```