# HLB-ABM

This repo contains the code required to run the agent-based model used in "Belief in Neighbor Behavior and Confidence in Scientific Information as Barriers to Cooperative Disease Control", published in the American Journal of Agricultural Economics.

This repository contains the files for both the model which generates the data, as well as the R code which was used to produce the results and figures used in the paper. The analysis code is confined to the `Analysis` folder, the rest of the repository contains source files for the agent-based model.

## Analysis

To run the analysis, the individual files can be run and the results coalesced. Many paths are hardcoded into the files, these must be changed in order for the scripts to run properly. The data can be acquired [here](https://hdl.handle.net/20.500.12741/rep:2267). The file names will have an equivalent to those in the scripts but the paths will need to be adjusted.

## Libraries

To use this model, please download boost and cereal and place them in the headers folder as "boost" and "cereal"

Boost: https://www.boost.org/users/download/
Cereal: https://uscilab.github.io/cereal/index.html

The file structure should look like the following:

```
headers
  -> boost
    -> accumulators
    -> algorithms
    -> etc...
  -> cereal
    -> archives
    -> details
    -> etc...
```

Once these folders are in place, the model can be compiled by running `make`, producing an executable called `serial.out`. To run the model use the following:

`<executable-name> <path/to/econConfig> <path/to/bioConfig>`

If no config parameters are provided, the defaults located at configs/econConfig.json and configs/bioConfigs.json will be used.
