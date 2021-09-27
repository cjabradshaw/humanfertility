# Human fertility models

<img align="right" src="contraception.png" alt="contraception" width="200" style="margin-top: 20px">

<a href="https://cran.r-project.org">R</a> code to reproduce models testing for cross-national human fertility patterns relative to underlying drivers

code updated September 2021

R code accompanies the following paper:

<strong><a href="https://globalecologyflinders.com/people/#DIRECTOR">Bradshaw, CJA</a>, C Perry, <a href="https://www.linkedin.com/in/chitra-maharani-saraswati-6bab3510b?originalSubdomain=au">CM Saraswati</a>, M Judge, <a href="https://research-repository.uwa.edu.au/en/persons/jane-heyworth">J Heyworth</a>, <a href="https://research-repository.uwa.edu.au/en/persons/peter-le-souef">PN Le Souëf</a></strong>. To be submitted. Lower infant mortality and access to contraception reduce fertility in low- and middle-income nations. <em>Lancet Global Health</em>

Prof Corey J. A. Bradshaw <br>
<a href="http://globalecologyflinders.com" target="_blank">Global Ecology</a>, <a href="http://flinders.edu.au" target="_blank">Flinders University</a>, Adelaide, Australia <br>
September 2021 <br>
<a href=mailto:corey.bradshaw@flinders.edu.au>e-mail</a> <br>


## Data
- basedata.csv = <em>main input data</em> (sourced from the following online databases:)
- matmort.csv = <em>maternal mortality data</em>
- pop.yr.csv = <em>national population data</em>
- continent.country2.csv = <em>naming file necessary for data merging</em>

## Requires the following libraries
- <code>lme4</code>
- <code>Hmisc</code>
- <code>ggplot2</code>
- <code>plotly</code>
- <code>nlme</code>
- <code>car</code>
- <code>dismo</code>
- <code>gbm</code>
- <code>rgeos</code>
- <code>rworldmap</code>
- <code>rworldxtra</code>
- <code>rcompanion</code>
- <code>SpatialEpi</code>
- <code>ggridges</code>
- <code>dplyr</code>
- <code>ggpubr</code>
- <code>plyr</code>
- <code>fields</code>
- <code>ncf</code>
- <code>AICcmodavg</code>
- <code>modEvA</code>

## and the following source-code files
- <code>new_lmer_AIC_tables3.R</code>
- <code>r.squared.R</code>
