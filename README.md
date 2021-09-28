# Human fertility models

<img align="right" src="contraception.png" alt="contraception" width="200" style="margin-top: 20px">

<a href="https://cran.r-project.org">R</a> code to reproduce models testing for cross-national human fertility patterns relative to underlying drivers

code updated September 2021

R code accompanies the following paper:

<strong><a href="https://globalecologyflinders.com/people/#DIRECTOR">Bradshaw, CJA</a>, C Perry, <a href="https://www.linkedin.com/in/chitra-maharani-saraswati-6bab3510b?originalSubdomain=au">CM Saraswati</a>, M Judge, <a href="https://research-repository.uwa.edu.au/en/persons/jane-heyworth">J Heyworth</a>, <a href="https://research-repository.uwa.edu.au/en/persons/peter-le-souef">PN Le Souëf</a></strong>. To be submitted. Lower infant mortality and access to contraception reduce fertility in low- and middle-income nations. <em>Lancet Global Health</em>

## Summary
# Background
Family planning has been associated with reduced fertility and improved maternal and child survival. Although fertility is declining in many nations, progress on improving maternal and child-health outcomes is slower. To assess why, we tested hypotheses regarding the potential predictors of variation in fertility among 46 low- and middle-income countries. We predicted that higher nation-level availability/quality of family-planning services, maternal education, and wealth reduce average fertility, whereas higher infant mortality and religious adherence increase it.
# Methods
Using publicly available data, we collated test and control variables in six themes: (<em>i</em>) availability of family planning, (<em>ii</em>) quality of family planning, (<em>iii</em>) maternal education, (<em>iv</em>) religion, (<em>v</em>) mortality, and (<em>vi</em>) socio-economic conditions. Given the sample size, we first constructed general linear models to test for relationships between fertility and the variables from each theme, from which we retained those with the highest explanatory power within a final general linear model set to determine the partial correlation of dominant test variables. We also applied boosted regression trees, generalised least-squares models, and a generalised linear mixed-effects models to account for non-linearity and spatial autocorrelation.
# Findings
On average among all countries, we found an association between all main variables and fertility, with reduced infant mortality having the strongest relationship with reduced fertility. Access to contraception was the next-highest correlate with reduced fertility, with female secondary education, home health visitations, and adherence to Catholicism having weak, if any, explanatory power.
# Interpretation
Our models suggest that decreasing infant mortality and increasing access to contraception will have the greatest effect on decreasing global fertility. We thus provide new evidence that progressing the United Nation’s Sustainable Development Goals for reducing infant mortality can be accelerated by increasing access to any form of family planning.

Prof Corey J. A. Bradshaw <br>
<a href="http://globalecologyflinders.com" target="_blank">Global Ecology</a> | <em><a href="https://globalecologyflinders.com/partuyarta-ngadluku-wardli-kuu/" target="_blank">Partuyarta Ngadluku Wardli Kuu</a></em>, <a href="http://flinders.edu.au" target="_blank">Flinders University</a>, Adelaide, Australia <br>
September 2021 <br>
<a href=mailto:corey.bradshaw@flinders.edu.au>e-mail</a> <br>

## Data
All data sourced from the following online databases: (<em>i</em>) <a href="http://dhsprogram.com">Demographic and Health Surveys</a>, (<em>ii</em>) <a href="http://track20.org/pages/data_analysis/policy/FPE.php">Family Planning Effort Index</a>, (<em>iii</em>) <a href="http://mics.unicef.org">Multiple Indicator Cluster Surveys</a>, (<em>iv</em>) <a href="http://track20.org/pages/data_analysis/policy/NCIFP.php">National Composite Index on Family Planning</a>, (<em>v</em>) <a href="http://data.worldbank.org">World Bank</a>, (<em>vi</em>) <a href="http://cia.gov/the-world-factbook">World Factbook</a>, and (<em>vii</em>) World Health Organization <a href="http://who.int/data/gho">Global Health Observatory</a>.
- <strong>basedata.csv</strong> = <em>main input data</em>
- <strong>matmort.csv</strong> = <em>maternal mortality data</em>
- <strong>pop.yr.csv</strong> = <em>national population data</em>
- <strong>continent.country2.csv</strong> = <em>naming file necessary for data merging</em>

## Requires the following R libraries
- <code><a href="https://cran.r-project.org/web/packages/lme4/index.html">lme4</a></code>
- <code><a href="https://cran.r-project.org/web/packages/Hmisc/index.html">Hmisc</a></code>
- <code><a href="https://ggplot2.tidyverse.org/">ggplot2</a></code>
- <code><a href="https://plotly.com/r/">plotly</a></code>
- <code><a href="https://cran.r-project.org/web/packages/nlme/index.html">nlme</a></code>
- <code><a href="https://cran.r-project.org/web/packages/car/index.html">car</a></code>
- <code><a href="https://cran.r-project.org/web/packages/dismo/index.html">dismo</a></code>
- <code><a href="https://cran.r-project.org/web/packages/gbm/index.html">gbm</a></code>
- <code><a href="https://cran.r-project.org/web/packages/rgeos/index.html">rgeos</a></code>
- <code><a href="https://cran.r-project.org/web/packages/rworldmap/index.html">rworldmap</a></code>
- <code><a href="https://cran.r-project.org/web/packages/rworldxtra/index.html">rworldxtra</a></code>
- <code><a href="https://cran.r-project.org/web/packages/rcompanion/index.html">rcompanion</a></code>
- <code><a href="https://cran.r-project.org/web/packages/SpatialEpi/index.html">SpatialEpi</a></code>
- <code><a href="https://cran.r-project.org/web/packages/ggridges/vignettes/introduction.html">ggridges</a></code>
- <code><a href="https://cran.r-project.org/web/packages/dplyr/vignettes/dplyr.html">dplyr</a></code>
- <code><a href="https://cran.r-project.org/web/packages/ggpubr/index.html">ggpubr</a></code>
- <code><a href="https://cran.r-project.org/web/packages/plyr/index.html">plyr</a></code>
- <code><a href="https://cran.r-project.org/web/packages/fields/index.html">fields</a></code>
- <code><a href="https://cran.r-project.org/web/packages/ncf/index.html">ncf</a></code>
- <code><a href="https://cran.r-project.org/web/packages/AICcmodavg/index.html">AICcmodavg</a></code>
- <code><a href="https://cran.r-project.org/web/packages/modEvA/index.html">modEvA</a></code>

## and the following source-code files
- <code>new_lmer_AIC_tables3.R</code>
- <code>r.squared.R</code>
