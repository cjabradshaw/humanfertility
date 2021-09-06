# Human fertility models

R code to reproduce models testing for cross-national human fertility patterns relative to underlying drivers

code updated September 2021

R code accompanies the following paper:

<strong>Perry, C, <a href="https://globalecologyflinders.com/people/#DIRECTOR">CJA Bradshaw</a>, CM Saraswati, M Judge, <a href="https://research-repository.uwa.edu.au/en/persons/jane-heyworth">J Heyworth</a>, <a href="https://research-repository.uwa.edu.au/en/persons/peter-le-souef">PN Le Souëf</a></strong>. In review. Lower infant mortality and access to contraception reduce fertility in low- and middle-income nations. <em>Lancet Global Health<em>

Prof Corey J. A. Bradshaw <br>
<a href="http://globalecologyflinders.com" target="_blank">Global Ecology</a>, <a href="http://flinders.edu.au" target="_blank">Flinders University</a>, Adelaide, Australia <br>
September 2021 <br>
<a href=mailto:corey.bradshaw@flinders.edu.au>e-mail</a> <br>


## Data:
- basedata.csv = main input data
- matmort.csv = maternal mortality data
- pop.yr.csv = national population data
- continent.country2.csv = naming file necessary for data merging

## Requires the following libraries:
- lme4
- Hmisc
- ggplot2
- plotly
- nlme
- car
- dismo
- gbm
- rgeos
- rworldmap
- rworldxtra
- rcompanion
- SpatialEpi
- ggridges
- dplyr
- ggpubr
- plyr
- fields
- ncf
- AICcmodavg
- modEvA

## and following source-code files:
- new_lmer_AIC_tables3.R
- r.squared.R
