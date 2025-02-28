# Code related to the paper "A neural mechanism for learning from delayed postingestive feedback"


This repository contains:
- Code to generate the figures in the paper "A neural mechanism for learning from delayed postingestive feedback"
- Code to fit generalized linear mixed models (GLMMs) to the data in the paper
- Code to run the brainwide cell detection pipeline (designed for Princeton computing clusters)


The cell detection pipeline incorporates a number of existing resources:
- ClearMap2 (https://github.com/ChristophKirst/ClearMap2)
- Cellfinder (https://github.com/brainglobe/cellfinder)
- BrainPipe (https://github.com/PrincetonUniversity/BrainPipe)
- Elastix (https://github.com/SuperElastix/elastix)
- Pystripe (https://github.com/chunglabmit/pystripe)
- Allen Mouse Brain Common Coordinate Framework (https://atlas.brain-map.org)


The GLMM code requires the following R packages:
- glmmTMB (https://github.com/glmmTMB/glmmTMB)
- marginaleffects (https://github.com/vincentarelbundock/marginaleffects)
- stringr (https://github.com/tidyverse/stringr)


Data: https://doi.org/10.6084/m9.figshare.28327118

bioRxiv: https://doi.org/10.1101/2023.10.06.561214