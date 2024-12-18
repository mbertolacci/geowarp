# GeoWarp example

Within R, run the following commands to install the required dependencies and the geowarp R package:

```{r}
install.packages(c(
  'dplyr', 'patchwork', 'ggplot2', 'patchwork', 'matrixStats', 'devtools'
))

devtools::install_github('mbertolacci/geowarp')
```

You can then run the two R scripts in this folder.

The first R script, `01-fit-vertical-only.R`, fits a model with no horizontal dependencies to simulated data. This can be useful for exploring the vertical structure of the data: how the mean and variance vary with depth, and whether there is vertical nonstationarity. This can also be used to tune some of the choices such as the number of vertical knots for the mean. These quantities are documented with comments in the file.

The second R script, '02-fit-3d.R`, generates the values of a full 3-D field and some soundings within that field. Then it fits GeoWarp to the soundings and attempts to predict the full 3-D field.
