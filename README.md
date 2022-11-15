# Blue Genes Manuscript

Code and data for *Schoenoplectus americanus* Blue Genes manuscript

All code was run in R version 4.1.2 (2021-11-01). Necessary packages to run script are identified at the top of each script.

All R scripts are written assuming the below organizational structure:
```
BlueGenes
│   README.md
└─── code
└─── supp_code
└─── supp_data
└─── figs
```

To generate figures using the files "Fig4.R", "Fig5.R", "Fig6.R", "Fig.7", or "SuppFigs.R", the respective (generalized) linear model objects must first be fit by running the code "TraitAnalysis.R", "CompetitionAnalysis.R", and/or "ExtinctAnalysis.R". Thus, we recommend running those three scripts in full prior to generating figures. 


