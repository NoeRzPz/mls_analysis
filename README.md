# Longevity analysis

## Exploratory Descriptive Analysis of phenotipcic traits

We create eda_r conda environment , to work with traits_eda.Rmd script

```
mamba env create --file envs/eda_r.yaml
```
## CAAStools

We have followed the instructions on `https://github.com/linudz/caastools`  to intall it.

We create conda environment with dependencies needed by CAAStools

```
mamba env create --file envs/caast_test.yaml
```

We can run localy the example provided with: 

```
bash scripts/local_caas_test.sh
```
