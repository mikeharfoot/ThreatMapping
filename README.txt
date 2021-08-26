This repository contains the code used for the analysis and visualisation of IUCN range maps and threat assessment classifications to generate maps of impact probability for the publication Harfoot et al., “Using the IUCN Red List to map threats to terrestrial vertebrates at global scale”

Data availability:
Data on range maps are freely available at https://www.iucnredlist.org/resources/spatial-data-download and http://datazone.birdlife.org/species/requestdis. IUCN threat classification assessment data can be downloaded using the Red List API (http://apiv3.iucnredlist.org/api/v3/docs) or on request from IUCN’s Global Species Programme Red List Unit (redlist@iucn.org). Other data are freely available using citations in the manuscript.

“SimulatingSpatialThreatMap.R”:
Contains code for simulating maps of impact with specified degrees of spatial autocorrelation. It also simulates the threat assessment process for real species ranges assuming different assessment protocols and degrees of uncertainty.

“Compare_Methods_RMSE-tidy.R”:
Runs a set of different models against the simulated impact probability maps to compare there predictive ability against those underlying impact maps.

“RegionalLogisticRegressionModels.R”
Runs logistic regression models for each pixel of a global 50km x 50km grid with separate models for the three taxonomic groups amphibians, birds and mammals, and for the six threats: agriculture, hunting & trapping, logging, invasive species, pollution, and climate change. The code saves the prediction probabilities for each grid cell to a ecoregion based directory structure.

“RegionalLogisticRegressionModelsDDThreatCodes.R”
Runs logistic regression models for each pixel of a global 50km x 50km grid with separate models for the three taxonomic groups amphibians, birds and mammals, and for the six threats: agriculture, hunting & trapping, logging, invasive species, pollution, and climate change. This code uses the drawn threat codes informed by the proportion of data deficient species within each group to assess the implications of this uncertainty on predicted threat probability. The prediction probabilities are saved for each grid cell to a ecoregion based directory structure.

“ActivitiesThreatsCorrelation.R”
Performs an evaluation between the predictions of the logistic regression models for mammals and for logging against the Global Forest Watch deforestation dataset.

"PlotGlobalMapsFromRegionalAllTaxaUncertainty.R"
Reads the predicted impact probabilities per pixel, taxaonomic group and threat and (1) visualises these as maps (including presenting the proportion of data deficient species within the group as a measure of uncertainty), (2) calculates impact risks and hotspots, (3) compares the predictions within a set of protected areas against the KBA database assessments of threats facing each KBA, and (4) compares the spatial distribution of probability of impact against current best estimates of human impact as measured by the human footprint index.

"CreateDataDeficientLayer.R"
Counts the number of data deficient species within each taxonomic group occuring in each grid cell and expresses this as the proportion of the total number of species occurring in that cell. These outputs are used to inform the uncertainty associated with predictions.

"ThreatMapFunctions.R"
Contains a set of functions used within other scripts for analysing and visualising threat maps.

"KBA_evaluation_Dissolve.R"
Reads the predicted impact probabilities per pixel, taxaonomic group and threat and compares the predicted impact probabilities with the assessed threat severity category for each human activity within KBAs using (1) Kendall's rank correlation and (2) comparing the predicted threat for high severity KBAs (Very rapid or moderate) with low severity (no or slow) KBAs. 

"DrawThreatCodes.R"
Draws a simulated set of threat codes for species using the proportion of data deficient species as a hypothetical indication of how uncertain Red List assessors might be of the threat status in each grid cell.