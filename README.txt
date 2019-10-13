This repository contains the code used for the analysis and visualisation of IUCN range maps and threat assessment classifications to generate maps of impact probability for the publication Harfoot et al., “Towards a global map of threats to species”

Data availability:
Data on range maps are freely available at https://www.iucnredlist.org/resources/spatial-data-download and http://datazone.birdlife.org/species/requestdis. IUCN threat classification assessment data can be downloaded using the Red List API (http://apiv3.iucnredlist.org/api/v3/docs) or on request from IUCN’s Global Species Programme Red List Unit (redlist@iucn.org). Other data are freely available using citations in the manuscript.

“SimulatingSpatialThreatMap.R”:
Contains code for simulating maps of impact with specified degrees of spatial autocorrelation. It also simulates the threat assessment process for real species ranges assuming different assessment protocols and degrees of uncertainty.

“Compare_Methods_RMSE-tidy.R”:
Runs a set of different models against the simulated impact probability maps to compare there predictive ability against those underlying impact maps.

“RegionalLogisticRegressionModels.R”
Runs logistic regression models for each pixel of a global 50km x 50km grid with separate models for the three taxonomic groups amphibians, birds and mammals, and for the six threats: agriculture, hunting & trapping, logging, invasive species, pollution, and climate change. The code saves the prediction probabilities for each grid cell to a ecoregion based directory structure.

“ActivitiesThreatsCorrelation-tidy.R”
Performs an evaluation between the predictions of the logistic regression models for mammals and for logging against the Global Forest Watch deforestation dataset.

PlotGlobalMapsFromRegionalAllTaxaMettByRealm-tidy.R
Reads the predicted impact probabilities per pixel, taxaonomic group and threat and (1) visualises these as maps, (2) calculates impact risks and hotspots, (3) compares the predictions within a set of protected areas against the METT database of expert-based assessment of threats facing the same protected areas, and (4) compares the spatial distribution of probability of impact against current best estimates of human impact as measured by the human footprint index.

