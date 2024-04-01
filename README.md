# insideFlow v1.1.0

## Overview

InsideFlow is an R toolkit for managing plugin processes from a FlowJo workspace,
developed and maintained by Joshua Luthy. 

## Installation 

You can install the development version from GitHub with `devtools`:

``` r
install.packages("devtools")
devtools::install_github("jluthy/insideFlow")
```

## Usage

InsideFlow is designed to work with FlowJo workspace data as it is exported from FlowJo during a plugin process. You can create a new insideFlow object and load a csv file of expression data exported from a FlowJo plugin process like this:

``` r
library(insideFlow)

# Create an insideFlow object
obj <- new("insideFlow", runID = "aloha_aina")

# Load the expression data from a FlowJo workspace
obj <- loadFlowJoCSV(obj, data_path = '/path/toFlowJo/data.csv', eventNumber = 'EventNumberDP')

```

Use the insideFlow object to perform various analyses on the data, such as clustering with FlowSOM, or batch normalization with cyCombine. Review the documentation for each function for more information on how to use it and what other slots are available in the insideFlow object.

## Contributing
Improvements and new features will be added on a regular basis, please post on
the [github page](https://github.com/jluthy/insideFlow) with any questions or if
you would like to contribute. 
