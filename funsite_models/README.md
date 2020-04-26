### Funsite predictor
The three jupyter notebooks were the used to generate the Funsite models and predictions for the three types of functional sites - CS, LIG and PPI. For getting LIG-metal predictions, use the LIG notebook and use the filter for restricting analysis to metal-binding sites.
The datasets used by all the jupyter notebooks can be found in `../datasets`. Some large dataset files could not be uploaded to this git repository, they can be downloaded from ftp://orengoftp.biochem.ucl.ac.uk/cath/supplementary-materials/2020-cath-funsite-predictor/

### Python 3 environment for Funsite predictor notebooks
```
virtualenv -p python3.6 FunsiteEnv
source FunsiteEnv/bin/activate

# install dependencies
pip install -r scripts/FunsiteEnv_requirements.txt
```
