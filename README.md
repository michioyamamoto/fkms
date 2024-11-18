# fkms: Functional K-means Clustering
R functions for the functional k-means clustering analysis developed by Yamamoto and Terada (2024), which estimates an optimal cluster structure of longitudinal data. This package efficiently estimates cluster structures even for irregular data with measurement points varying across individuals or sparse data that include subjects with extremely limited observations.

To install a package on GitHub, run the following code on the R console. Note that the following code may need to be run with administrative privileges.

```
install.packages("devtools") ## if the devtools package is not installed
library(devtools)
install_github("michioyamamoto/fkms").
```
