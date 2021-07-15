Aanalyses and figures to be published in Ecology Letters. "Hydraulic prediction of drought-induced plant dieback and top-kill depends on leaf habit and growth form" by YJ. Chen, B. Choat, F. Sterck, P. Maenpuen, M. Katabuchi, SB. Zhang, K. W. Tomlinson, R. S. Oliveira, YJ. Zhang, JX. Shen, KF. Cao, and S. Jansen .

## Data

- `data/traits.csv`: the raw dataset used in the paper
- `data/PCA_Log.csv`: log-transfomred dataset

## Code

- `table.Rmd`: R markdown to generate Table 1. `Rscript -e "rmarkdown::render('table.Rmd')"`

- `PCA.R`: R code to generate Fig. 2 `Rscript PCA.R`

- `Fig_regression.Rmd`: R markdown to generate Fig. 4-5 and Fig. S9-11 `Rscript -e "rmarkdown::render('Fig_regression.Rmd')"`
