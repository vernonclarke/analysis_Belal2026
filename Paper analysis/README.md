## Processed datasets and figure outputs

**Cholinergic interneuron control of intrastriatal GABAergic circuits targeting spiny
projection neurons is disrupted in Parkinson's disease models**

Belal M, Perez-Rosello T, Guven E.B, Kocaturk S, Xie Z, Ilijic E, Tkatch T, Li J, Dauer W,
Assous M, Tepper J.M, Clarke V.R.J, Surmeier D.J

This folder contains processed data spreadsheets and SVG figure panels generated from the raw
electrophysiology, imaging and in-situ-hybridisation recordings

## Raw data and analysis code availability

- Raw NWB data are available from the DANDI Archive as `DANDI:001832`  
  `https://doi.org/10.48324/dandi.001832/0.260611.2102`

- Analysis code is available on GitHub  
  `https://github.com/vernonclarke/analysis_Belal2026`

- The analysis code is archived on Zenodo  
 `https://doi.org/10.5281/zenodo.20658500`

- Computational modeling figures, Figure 4, Figure 5, Figure 6 and Figure S1-S4, are not included in this processed data folder.
  See the GitHub repository and Zenodo analysis-code archive listed below for these workflow
  
- Computational modeling code is available on GitHub
`https://github.com/vernonclarke/msNEURON_Belal2026`

- The computational modeling code is archived on Zenodo  
`https://doi.org/10.5281/zenodo.20705696`


## Licence

This processed dataset is released under **CC BY 4.0**

You may share and adapt the material for any purpose, including commercially, provided you give
appropriate credit to the publication above

## Figure folders

Each figure folder has its own README with the scripts, 
processed data files, SVG panels and statistics for that figure

`Figure 1`:   gabazine sensitivity of ChI-evoked PSCs in dSPNs and iSPNs    
`Figure 2`:   ChI-evoked PSCs compared between dSPNs and iSPNs  
`Figure 3`:   GABAA delta-subunit CRISPR knockdown and the slow PSC component  
`Figure 3B`:  NGF-interneuron-evoked GABAergic currents in SPNs  
`Figure 7`:   GRAB-ACh3.0 imaging of ChI-evoked ACh release in MCI-Park and control mice  
`Figure 8`:   ChI-evoked GABAergic PSCs in MCI-Park and control mice  
`Figure 9`:   ChI-evoked GABAergic input after 6-OHDA MFB lesion  
`Figure 10`:  NPY-interneuron-evoked input after 6-OHDA MFB lesion  
`Figure 11`:  RNAscope and ChI to NGF PSC analysis, including Figure S6 RNAscope data  
`Figure S5`:  independent 6-OHDA ChI-evoked PSC study  

## Folder layout

Most figure folders contain:

```text
Figure X/
  *.R or *.ipynb analysis scripts
  xlsx/ processed spreadsheets and csv tables
  svg/ figure panels and fit diagnostics
```

## Workflow

For electrophysiology figures, the usual workflow is:

1. Download the raw NWB files from DANDI Archive `DANDI:001832`
2. Run the figure data-processing script to create averaged traces in `xlsx/`
3. Run the figure analysis script to fit traces and create summary spreadsheets `xlsx` and `csv` formats
4. Run the figure graphs script to plot SVG panels and write statistics tables

Figure 7 imaging and Figure 11 RNAscope use related workflows described in their figure READMEs
