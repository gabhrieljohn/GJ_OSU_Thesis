---
editor_options: 
  markdown: 
    wrap: 72
---

# FES524 Final Project Data README
NOTE: Viewing this document in Visual mode in R Studio will automatically alter 
formatting, please only view in Source. 

In "Analyses/data_prep.R", the object "data_prep_422_230" is loaded.
This object is so named because it contains select relevant columns from
two datasets (dataset DSCMET422 and PRIMET230) derived from the HJ
Andrews Provisional Data Portal, URL for which is:
<https://andrewsforest.oregonstate.edu/data/streaming/provisional-data-portal>

These variables have been selected for their potential value in analyses
exploring relationships between tree growth and various climate
variables such as precipitation, soil moisture, and soil temperature.
This README may be updated to include additional variables of
exploratory interest.

Below is a table describing each variable. 

| Name                            | Definition                                                                              | Data Type | Missing Value | Unit        | Precision/Range                                                                  |
| ------------------------------- | --------------------------------------------------------------------------------------- | --------- | ------------- | ----------- | -------------------------------------------------------------------------------- |
| Site                            | Site from which data are derived                                                        | chr       | –             | –           | –                                                                                |
| Date                            | Timestamp of corresponding data                                                         | POSIXct   | –             | –           | 5 minutes                                                                        |
| DEND309\_AVG\_4200\_0\_01       | Band dendrometer reading of old-growth Douglas-fir ID "309" 42 meters above the ground  | num       | NaN           | micron (μm) | Resolution on CR1000: 3.3 μm. Accuracy: ±0.12%. Recorded to nearest whole number |
| Flag\_DEND309\_AVG\_4200\_0\_01 | QA/QC flag for the corresponding dendrometer                                            | chr       | NA            | –           | –                                                                                |
| DEND310\_AVG\_150\_0\_01        | Band dendrometer reading of old-growth Douglas-fir ID "310" 1.5 meters above the ground | num       | NaN           | micron (μm) | Resolution on CR1000: 3.3 μm. Accuracy: ±0.12%. Recorded to nearest whole number |
| Flag\_DEND310\_AVG\_150\_0\_01  | QA/QC flag for the corresponding dendrometer                                            | chr       | NA            | –           | –                                                                                |
| DEND311\_AVG\_150\_0\_01        | Band dendrometer reading of old-growth Douglas-fir ID "311" 1.5 meters above the ground | num       | NaN           | micron (μm) | Resolution on CR1000: 3.3 μm. Accuracy: ±0.12%. Recorded to nearest whole number |
| Flag\_DEND311\_AVG\_150\_0\_01  | QA/QC flag for the corresponding dendrometer                                            | chr       | NA            | –           | –                                                                                |
| DEND646\_AVG\_150\_0\_01        | Band dendrometer reading of old-growth western hemlock ID "646" 1.5 meters above ground | num       | NaN           | micron (μm) | Resolution on CR1000: 3.3 μm. Accuracy: ±0.12%. Recorded to nearest whole number |
| Flag\_DEND646\_AVG\_150\_0\_01  | QA/QC flag for the corresponding dendrometer                                            | chr       | NA            | –           | –                                                                                |
| DEND648\_AVG\_150\_0\_01        | Band dendrometer reading of old-growth western hemlock ID "648" 1.5 meters above ground | num       | NaN           | micron (μm) | Resolution on CR1000: 3.3 μm. Accuracy: ±0.12%. Recorded to nearest whole number |
| Flag\_DEND648\_AVG\_150\_0\_01  | QA/QC flag for the corresponding dendrometer                                            | chr       | NA            | –           | –                                                                                |
| DEND652\_AVG\_150\_0\_01        | Band dendrometer reading of old-growth western hemlock ID "652" 1.5 meters above ground | num       | NaN           | micron (μm) | Resolution on CR1000: 3.3 μm. Accuracy: ±0.12%. Recorded to nearest whole number |
| Flag\_DEND652\_AVG\_150\_0\_01  | QA/QC flag for the corresponding dendrometer                                            | chr       | NA            | –           | –                                                                                |
| DEND654\_AVG\_150\_0\_01        | Band dendrometer reading of second-growth Douglas-fir ID "654" 1.5 meters above ground  | num       | NaN           | micron (μm) | Resolution on CR1000: 3.3 μm. Accuracy: ±0.12%. Recorded to nearest whole number |
| Flag\_DEND654\_AVG\_150\_0\_01  | QA/QC flag for the corresponding dendrometer                                            | chr       | NA            | –           | –                                                                                |
| DEND655\_AVG\_150\_0\_01        | Band dendrometer reading of second-growth Douglas-fir ID "655" 1.5 meters above ground  | num       | NaN           | micron (μm) | Resolution on CR1000: 3.3 μm. Accuracy: ±0.12%. Recorded to nearest whole number |
| Flag\_DEND655\_AVG\_150\_0\_01  | QA/QC flag for the corresponding dendrometer                                            | chr       | NA            | –           | –                                                                                |
| DEND656\_AVG\_150\_0\_01        | Band dendrometer reading of second-growth Douglas-fir ID "656" 1.5 meters above ground  | num       | NaN           | micron (μm) | Resolution on CR1000: 3.3 μm. Accuracy: ±0.12%. Recorded to nearest whole number |
| Flag\_DEND656\_AVG\_150\_0\_01  | QA/QC flag for the corresponding dendrometer                                            | chr       | NA            | –           | –                                                                                |
| DEND860\_AVG\_150\_0\_01        | Band dendrometer reading of second-growth Douglas-fir ID "860" 1.5 meters above ground  | num       | NaN           | micron (μm) | Resolution on CR1000: 3.3 μm. Accuracy: ±0.12%. Recorded to nearest whole number |
| Flag\_DEND860\_AVG\_150\_0\_01  | QA/QC flag for the corresponding dendrometer                                            | chr       | NA            | –           | –                                                                                |
| ts                              | Same as "Date" above, renamed for consistent coding                                     | POSIXct   | –             | micron (μm) | 5 minutes                                                                        |
| yr                              | Timestamp column with year extracted                                                    | factor    | –             | year        | years                                                                            |
| DOY                             | Day of Year                                                                             | num       | –             | –           | days                                                                             |
| SOILTEMP\_MEAN\_010\_01         | Mean soil temperature 10 cm belowground                                                 | num       | NA            | °C          | To nearest hundredth                                                             |
| Flag\_SOILTEMP\_MEAN\_010\_01   | QA/QC flag for soil temperature probe                                                   | chr       | NA            | –           | –                                                                                |
| SOILTEMP\_MEAN\_020\_02         | Mean soil temperature 20 cm belowground                                                 | num       | NA            | °C          | To nearest hundredth                                                             |
| Flag\_SOILTEMP\_MEAN\_020\_02   | QA/QC flag for soil temperature probe                                                   | chr       | NA            | –           | –                                                                                |
| SOILTEMP\_MEAN\_050\_03         | Mean soil temperature 50 cm belowground                                                 | num       | NA            | °C          | To nearest hundredth                                                             |
| Flag\_SOILTEMP\_MEAN\_050\_03   | QA/QC flag for soil temperature probe                                                   | chr       | NA            | –           | –                                                                                |
| SOILTEMP\_MEAN\_100\_04         | Mean soil temperature 100 cm belowground                                                | num       | NA            | °C          | To nearest hundredth                                                             |
| Flag\_SOILTEMP\_MEAN\_100\_04   | QA/QC flag for soil temperature probe                                                   | chr       | NA            | –           | –                                                                                |
| SWC\_MEAN\_010\_01              | Mean soil water content 10 cm belowground                                               | num       | NA            | fraction    | To nearest thousandth                                                            |
| Flag\_SWC\_MEAN\_010\_01        | QA/QC flag for soil moisture probe                                                      | chr       | NA            | –           | –                                                                                |
| SWC\_MEAN\_020\_02              | Mean soil water content 20 cm belowground                                               | num       | NA            | fraction    | To nearest thousandth                                                            |
| Flag\_SWC\_MEAN\_020\_02        | QA/QC flag for soil moisture probe                                                      | chr       | NA            | –           | –                                                                                |
| SWC\_MEAN\_050\_03              | Mean soil water content 50 cm belowground                                               | num       | NA            | fraction    | To nearest thousandth                                                            |
| Flag\_SWC\_MEAN\_050\_03        | QA/QC flag for soil moisture probe                                                      | chr       | NA            | –           | –                                                                                |
| SWC\_MEAN\_100\_04              | Mean soil water content 100 cm belowground                                              | num       | NA            | fraction    | To nearest thousandth                                                            |
| Flag\_SWC\_MEAN\_100\_04        | QA/QC flag for soil moisture probe                                                      | chr       | NA            | –           | –                                                                                |
| PRECIP\_TOT\_100\_0\_01         | Total precipitation (TE525M rain gauge, 9.7in; TR525 retired 2023-06-07)                | num       | NA            | mm          | To nearest thousandth                                                            |
| Flag\_PRECIP_TOT_100\_0\_01     | QA/QC flag for the corresponding precipitation probe                                    | chr       | NA            | -           | -                                                                                |
| PRECIP\_ACC\_100\_0\_01         | Cumulative sum of total precipitation above                                             | num       | NA            | mm          | To the nearest hundredth                                                         |
| Flag\_PRECIP\_ACC\_100\_0\_01   | QA/QC flag for the corresponding precipitation probe                                    | chr       | NA            | -           | -                                                                                |

: Metadata for "data_prep_422_230"
