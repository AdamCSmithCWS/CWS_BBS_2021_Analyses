# CWS 2021 data-version Trend Analysis

Canadian Wildlife Service, BBS status and trend analyses for 2021 data version.

This analysis was done using bbsBayes (<https://github.com/bbsBayes/bbsBayes>), and some additional functions and code to fit the models in Stan, instead of JAGS.

This repo includes some git commits that tracked initial work to support the new version of the R-package bbsBayes2 (<https://github.com/bbsBayes/bbsBayes2>).

These analyses were conducted using Stan, and versions of the models that are included in the new bbsBayes2 package.

Scripts that are numbered 1, 1a, 2, 3, 4, 5, and 6, include the full workflow from raw BBS data to the trend estimates, trend-maps, population trajectories, and rolling short-term trends necessary to recreate the 2021 CWS trend analyses. THe summarized output from this analysis is available on a shareable [Google Drive Folder](https://drive.google.com/drive/folders/134W6w4yvMhVrlrXyG17uwE96d4u_flSf?usp=share_link). 

If you are interested in a customized summary of trends (e.g., different regions or time-periods), please contact me directly and I can share the full model output for use with the bbsBayes2 package. You're of course, welcome to re-fit the model yourself, but the model fitting process takes significant computing resources (time, memory, and processors).

