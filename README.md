# Human fertility models

<a href="https://doi.org/10.5281/zenodo.7496142"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.7496142.svg" alt="DOI"></a>

<img align="right" src="www/contraception.png" alt="contraception" width="200" style="margin-top: 20px">

<a href="https://cran.r-project.org"><em>R</em></a> code to reproduce models testing for cross-national human fertility patterns relative to underlying drivers

code updated December 2022

<em>R</em> code accompanies the following paper:

<strong><a href="https://globalecologyflinders.com/people/#DIRECTOR">Bradshaw, CJA</a>, C Perry, <a href="https://orcid.org/0000-0002-9948-1865">M Judge</a>, <a href="https://www.linkedin.com/in/chitra-maharani-saraswati-6bab3510b?originalSubdomain=au">CM Saraswati</a>, <a href="https://research-repository.uwa.edu.au/en/persons/jane-heyworth">J Heyworth</a>, <a href="https://research-repository.uwa.edu.au/en/persons/peter-le-souef">PN Le Souëf</a></strong>. 2023. <a href="http://doi.org/10.1371/journal.pone.0280260">Lower infant mortality, lower household size, and more access to contraception reduce fertility in low- and middle-income nations</a>. <em>PLoS One</em> doi:10.1371/journal.pone.0280260  

Pre-print (out-of-date) also available <a href="http://doi.org/10.1101/2021.12.16.21267946">here</a>.

## Abstract
Although average contraceptive use has increased globally in recent decades, an estimated 222 million (26%) of women of child-bearing age worldwide face an unmet need for family planning — defined as a discrepancy between fertility preferences and contraception practice, or failing to translate desires to avoid pregnancy into preventative behaviours and practices. While many studies have reported relationships between availability/quality of contraception and family planning, infant mortality, and fertility, these relationships have not been evaluated quantitatively across a broad range of low- and middle-income countries. Using publicly available data from 64 low- and middle-income countries, we collated test and control variables in six themes: (<em>i</em>) availability of family planning, (<em>ii</em>) quality of family planning, (<em>iii</em>) female education, (<em>iv</em>) religion, (<em>v</em>) mortality, and (<em>vi</em>) socio-economic conditions. We predicted that higher nation-level availability/quality of family-planning services and female education reduce average fertility, whereas higher infant mortality, great household size (a proxy for population density), and religious adherence increase it. Given the sample size, we first constructed general linear models to test for relationships between fertility and the variables from each theme, from which we retained those with the highest explanatory power within a final general linear model set to determine the partial correlation of dominant test variables. We also applied boosted regression trees, generalised least-squares models, and generalised linear mixed-effects models to account for non-linearity and spatial autocorrelation. On average among all countries, we found the strongest associations between fertility and infant mortality, household size, and access to any form of contraception. Higher infant mortality and household size increased fertility, whereas greater access to any form of contraception decreased it fertility. Female education, home visitations by health workers, quality of family planning, and religious adherence all had weak, if any, explanatory power. Our models suggest that decreasing infant mortality, ensuring sufficient housing to reduce household size, and increasing access to contraception will have the greatest effect on decreasing global fertility. We thus provide new evidence that progressing the United Nation’s Sustainable Development Goals for reducing infant mortality can be accelerated by increasing access to family planning.  

<br>
Prof Corey J. A. Bradshaw <br>
<a href="http://globalecologyflinders.com" target="_blank">Global Ecology</a> | <em><a href="https://globalecologyflinders.com/partuyarta-ngadluku-wardli-kuu/" target="_blank">Partuyarta Ngadluku Wardli Kuu</a></em>, <a href="http://flinders.edu.au" target="_blank">Flinders University</a>, Adelaide, Australia <br>
August 2022 <br>
<a href=mailto:corey.bradshaw@flinders.edu.au>e-mail</a> <br>

## <a href="https://github.com/cjabradshaw/humanfertility/tree/main/scripts">Scripts</a>
- main script <code>humanFertilityGitHubV3.R</code> includes all data preparation and modelling steps
- <code>new_lmer_AIC_tables3.R</code> (source file)
- <code>r.squared.R</code> (source file)

## <a href="https://github.com/cjabradshaw/humanfertility/tree/main/data">Data</a>
All data sourced from the following online databases: (<em>i</em>) <a href="http://dhsprogram.com">Demographic and Health Surveys</a>, (<em>ii</em>) <a href="http://track20.org/pages/data_analysis/policy/FPE.php">Family Planning Effort Index</a>, (<em>iii</em>) <a href="http://mics.unicef.org">Multiple Indicator Cluster Surveys</a>, (<em>iv</em>) <a href="http://track20.org/pages/data_analysis/policy/NCIFP.php">National Composite Index on Family Planning</a>, (<em>v</em>) <a href="http://data.worldbank.org">World Bank</a>, (<em>vi</em>) <a href="http://wid.world/data">World Inequality Database</a>, and (<em>vii</em>) World Health Organization <a href="http://who.int/data/gho">Global Health Observatory</a>, <a href="http://www.thearda.com">Association of Religion Data Archives</a>.
- <strong>basedata.update.csv</strong> = <em>main input data</em>
- <strong>matmort.update.csv</strong> = <em>maternal mortality data</em>
- <strong>pop.yr.update.csv</strong> = <em>national population data</em>
- <strong>continent.country2.csv</strong> = <em>naming file necessary for data merging</em>

## Requires the following <em>R</em> libraries
<code><a href="https://cran.r-project.org/web/packages/lme4/index.html">lme4</a></code>, <code><a href="https://cran.r-project.org/web/packages/Hmisc/index.html">Hmisc</a></code>, <code><a href="https://ggplot2.tidyverse.org/">ggplot2</a></code>, <code><a href="https://plotly.com/r/">plotly</a></code>, <code><a href="https://cran.r-project.org/web/packages/nlme/index.html">nlme</a></code>, <code><a href="https://cran.r-project.org/web/packages/car/index.html">car</a></code>, <code><a href="https://cran.r-project.org/web/packages/dismo/index.html">dismo</a></code>, <code><a href="https://cran.r-project.org/web/packages/gbm/index.html">gbm</a></code>, <code><a href="https://cran.r-project.org/web/packages/rgeos/index.html">rgeos</a></code>, <code><a href="https://cran.r-project.org/web/packages/rworldmap/index.html">rworldmap</a></code>, <code><a href="https://cran.r-project.org/web/packages/rworldxtra/index.html">rworldxtra</a></code>, <code><a href="https://cran.r-project.org/web/packages/rcompanion/index.html">rcompanion</a></code>, <code><a href="https://cran.r-project.org/web/packages/SpatialEpi/index.html">SpatialEpi</a></code>, <code><a href="https://cran.r-project.org/web/packages/ggridges/vignettes/introduction.html">ggridges</a></code>, <code><a href="https://cran.r-project.org/web/packages/dplyr/vignettes/dplyr.html">dplyr</a></code>, <code><a href="https://cran.r-project.org/web/packages/ggpubr/index.html">ggpubr</a></code>, <code><a href="https://cran.r-project.org/web/packages/plyr/index.html">plyr</a></code>, <code><a href="https://cran.r-project.org/web/packages/fields/index.html">fields</a></code>, <code><a href="https://cran.r-project.org/web/packages/ncf/index.html">ncf</a></code>, <code><a href="https://cran.r-project.org/web/packages/AICcmodavg/index.html">AICcmodavg</a></code>, <code><a href="https://cran.r-project.org/web/packages/modEvA/index.html">modEvA</a></code>, <code><a href="https://cran.r-project.org/web/packages/VIM/index.html">VIM</a></code>, <code><a href="https://www.rdocumentation.org/packages/mice/versions/3.14.0/topics/mice">mice</a></code>. <code><a href="https://indrajeetpatil.github.io/ggstatsplot/">ggstatsplot</a></code>
<br>
<br>
<a href="https://www.flinders.edu.au"><img align="bottom-left" src="www/Flinders_University_Logo_Horizontal_RGB_Master.png" alt="Flinders University logo" width="200" style="margin-top: 20px"></a><a href="https://github.com/FutureChildHealth"><img align="bottom-left" src="https://github.com/cjabradshaw/humanfertility/blob/main/www/FCHlogoFinaltransp.png" alt="Future Child Health logo" width="200" style="margin-top: 20px"></a>
<a href="https://globalecologyflinders.com"><img align="bottom-left" src="www/GEL Logo Kaurna New Transp.png" alt="GEL logo" width="200" style="margin-top: 20px"></a>
<a href="https://www.uwa.edu.au"><img align="bottom-left" src="www/uwa2.png" alt="UWA logo" width="150" style="margin-top: 20px"></a>

