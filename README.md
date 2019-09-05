# rotavirusNYC
MATLAB code for model describing rotavirus transmission dynamics and impact of vaccination in New York City

WORK IN PROGRESS: Model code is now posted. Data will be uploaded prior to publication.

CONTENTS:
rotamodelNYmo.m - Simulates rotavirus transmission model and compared to NYC hospitalization and laboratory-reported data.
rasisAM.m - Function file containing the differential equations for the model. Called by 'rotamodelNYmo'.
monthly_vcov_nyc.csv - Data file containing monthly rotavirus vaccine coverage with 1 dose (first column) or full 2- or 3-dose series (second column)
monthly_RVhosp_nyc.csv - Data file containing monthly rotavirus hospitalizations for all ages (first column) and by age group (columns 2-5: <2 yrs, 2-4 yrs, 5-12 yrs, 13+ yrs)
monthly_RVlab_nyc.csv - Data file containing monthly laboratory-reported rotavirus cases by age group (columns 1-4: <2 yrs, 2-4 yrs, 5-12 yrs, 13+ yrs) and for all ages (last column)
