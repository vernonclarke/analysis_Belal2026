# Analysis instructions for ACHGRAB and RNAscope data

## Index
- [Setup Instructions](#setup-instructions)
- [Create .env for token](#create-env-for-token)
- [Quick Checks](#quick-checks)
- [Export environment](#export-environment)
- [Run Jupyter Notebook](#run-jupyter-notebook)
- [Download DANDI dataset](#download-dandi-dataset)
- [Ensure Git ignore](#ensure-git-ignore)
- [Create html](#create-html)



## Setup Instructions
```bash
cd "/Users/euo9382/Documents/Repositories/analysis_Belal2026"
conda create -n image_analysis python=3.12 -y
conda activate image_analysis
pip install matplotlib seaborn python-dotenv remfile
pip install notebook jupyterlab
pip install plotly
pip install -U kaleido
pip install ipykernel
pip install ipympl
pip install pytz


```

---

## Create .env for token
```bash
cd /Users/euo9382/Documents/Repositories/CatalystNeuro
echo 'DANDI_API_TOKEN=your_actual_token_here' > .env
head .env
echo ".env" >> .gitignore
```

---

## Quick Checks
```bash
python -c "import pynwb, hdmf; print('pynwb:', pynwb.__version__, 'hdmf:', hdmf.__version__)"
python -c "import importlib.metadata; print('neuroconv:', importlib.metadata.version('neuroconv'))"
python -c "import surmeier_lab_to_nwb as sln; print('pkg OK')"
```

---

## Export environment
```bash
conda env export --from-history > environment.yml
conda env export > environment.lock.yml
pip freeze > requirements.txt
```

---

## Run Jupyter Notebook
```bash
cd "/Users/euo9382/Documents/Repositories/analysis_Belal2026"
conda activate image_analysis
jupyter notebook
```

---

## Download DANDI dataset
```bash
export $(grep DANDI_API_TOKEN /Users/euo9382/Documents/Repositories/CatalystNeuro/.env | sed 's/DANDI_API_TOKEN/DANDI_API_KEY/')
dandi download DANDI:001538/draft -o '/Users/euo9382/Documents/Repositories/CatalystNeuro/Zhai Paper'
```

---

## Ensure Git ignore
```bash
cd /Users/euo9382/Documents/Repositories/CatalystNeuro
echo "Zhai paper/" >> .gitignore
git status --ignored
```


## Create html
```bash
cd /Users/euo9382/Documents/Repositories/CatalystNeuro/notebooks
jupyter nbconvert figure_5F_acetylcholine_biosensor_VC.ipynb --to html --output-dir=docs
```
