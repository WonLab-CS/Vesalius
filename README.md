[![DOI](https://zenodo.org/badge/306332649.svg)](https://zenodo.org/badge/latestdoi/306332649)

# Vesalius

Welcome to the Vesalius GitHub page!
<img src="man/figures/banner.png" />

UNDER CONSTRUCTION - NEW VERSION COMMING SOON


## What is Vesalius ?
Vesalius is an R package to decipher tissue anatomy by embracing various
image analysis techniques for high-resolution ST data. Vesalius identifies
spatially expressed genes linked to the morphology of tissue structures.

## How to install Vesalius?

If you do not have `devtools` already installed, please do so beforehand.

```
install.packages("devtools")
```  
Ensure that the library has been loaded
```
library(devtools)
```
Install Vesalius via GitHub
```
install_github("patrickCNMartin/Vesalius")
```

## How to use Vesalius?
Vesalius provides an internal data set taken from real Spatial Transcriptomic
data. This can be used as a dummy data set to get a feel for the Vesalius
workflow.

The "Quick Start" guide can be found [here](https://github.com/WonLab-CS/Vesalius/blob/main/vignettes/vesalius.Rmd)

An in depth view of the Vesalius workflow can be found [here](https://github.com/WonLab-CS/Vesalius/blob/main/vignettes/Vesalius_Analysis/Vesalius_MSB_analysis.Rmd). This contains the entire analysis related to
the Vesalius [Molecular Systems Biology](https://www.embopress.org/doi/full/10.15252/msb.202211080).

NOTE: The in depth analysis file, contains path to files that should be changed
accordingly.


## What's next?
The Vesalius package is in its early stage of development. We would ask you to
share with us any bugs, concerns, or features you wish to see implemented.

Please open a GitHub issue or send an email to Patrick Martin (Patrick.Martin@cshs.org - pcnmartin@gmail.com)
