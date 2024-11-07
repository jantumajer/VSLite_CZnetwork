# VSLite_CZnetwork
Data supporting the manuscript using the VS-Lite model to forecast intra-annual growth kinetics and phenology under climate change for temperate forests in Central Europe

### Rscripts 
The key function implementing the VS-Lite model is `RUN.R`. This script can be used to load predefined functions from additional .R files, load the input data, calibrate the model and save the results for the baseline period. To run a model in a forecasting model, use `Forecast_modeling.R`. We acknowledge [Suzan Tolwinski-Ward](https://github.com/suztolwinskiward/VSLiteR) who defined the model and published its original scripts. We modified them according to recommendations of the recent studies including *Campelo et al.* [(2018)](https://doi.org/10.1016/j.dendro.2018.03.001), *Seftigen et al.* [(2018)](https://onlinelibrary.wiley.com/doi/10.1111/geb.12802), *Tumajer et al.* [(2017)](http://linkinghub.elsevier.com/retrieve/pii/S0168192317302290), and *Tumajer et al.* [(2023)](https://linkinghub.elsevier.com/retrieve/pii/S0048969723057807).

### Climatic data
Monthly mean temperature and precipitation totals for each site from January 1961 to December 2020. 

### Climate projections
Decadal anomalies of mean monthly temperature and precipitation totals concerning the normal period 1995-2014.

### NDVI
Normalized Difference Vegetation Index for each site for the 2000-2023 period.

### Outputs
A standard set of outputs of the VS-Lite model for a baseline calibration period and four scenarios of future climate change. For each site and simulation, the `.Rda` file stores a list with the following tables and variable:

- **mod.trw** Modeled tree-ring width chronology
- **gT** Matrix of simulate partial growth rates to temperature
- **gM** Matrix of simulated partial growth rates to soil moisture
- **gINT** Matrix of simulated integral growth rates
- **gE** Ration of mean monthly daylength to daylength of summer solstice (=partial growth rate to photoperiod)
- **Moist** Simulated monthly soil moisture [v/v]
- **Temp** Mean monthly air temperature [°C]
- **syear** First year of simulations
- **eyear** Last year of simulations
- **par** Model parameters
  - **T1, T2, T3, T4** Parameters defining the shape of the response function to temperature
  - **M1, M2** Parameters defining the shape of the response function to temperature
  - **M3, M4** Not-applied parameters defining a shape of the response function to temperature
  - **Acor, I_0** Parameters of previous-year autocorrelation effects
  - **corr** Correlation coefficient between observed and simulated chronology (only available for outputs based on baseline period)

