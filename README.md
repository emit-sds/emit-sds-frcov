# EMIT Fractional Cover Product Delivery - Algorithm Theoretical Basis Document (ATBD)

*K.D. Chadwick**<sup>1</sup>, *Red Willow Coleman*<sup>1</sup>

<sup>1</sup>Jet Propulsion Laboratory, California Institute of Technology

Corresponding author: K. Dana Chadwick (dana.chadwick@jpl.nasa.gov)

**Key Points:**
   1. Please note that this ATBD will be updated on an ongoing basis as the EMIT extended mission progresses. This is intended to be a place where the community can find the most up-to-date information on the current plans for algorithm development and offer contributions.
   2. 
   3. 

**Version:** 1.0

**Release Date:** TBD

**DOI:** TBD

## Abstract

## Plain Language Summary

### Keywords: hyperspectral imaging, imaging spectroscopy, fractional cover, terrestrial

## 1 Version Description

This is Version 1.0 of the EMIT fractional cover and fractional cover quality assessment (QA) product.

## 2 Introduction

## 3 Context/Background

### 3.1 Historical Perspective

### 3.2 Additional Information

## 4 Algorithm Description

### 4.1 Scientific Theory

### 4.2 Mathematical Theory

### 4.3 Fractional Cover Algorithm Input Variables

The required input files for fractional cover production are in Table 1.

**Table 1.** _Input variables_
| Name | Description | Unit | Required |
| --- | --- | --- | --- |
| HDRF (Surface Reflectance) | hemispherical-directional reflectance factor per wavelength | unitless | true |
| Endmember Spectral Library | collection of spectral surface endmembers (HDRF per wavelength) | unitless | true |
| Observation Geometry | solar zenith angle, view zenith angle, relative azimuth angle | degree | false |

### 4.4 Fractional Cover Algorithm Output Variables
The EMIT output data products delivered to the DAAC use their formatting conventions, the system operates internally on data products stored as binary data cubes with detatched human-readable ASCII header files.

### 4.5 Fractional Cover QA Product Input Variables
The required input files for fractional cover QA production are in Table 2.

**Table 2.** _QA input variables_
| Name | Description | Spatial Resolution | Required | Type |
| --- | --- | --- | --- | --- |
| EMIT L2A Surface Reflectance | hemispherical-directional reflectance factor (HDRF) per wavelength | 60 m | true | raster
| EMIT L2A Mask | estimated surface reflectance uncertainty and masks | 60 m | true | raster
| ESA WorldCover | global landcover classification derived from Sentinel-1 and Sentinel-2 data | 10 m | true | raster
| GSHHG | the Global Self-consistent, Hierarchical, High-resolution Geography Database for coastlines and rivers | <60 m | true | vector

#### 4.5.1 Clouds

Flag all pixels identified as "cirrus" or "cloud" by the EMIT L2A Mask product (Green, 2022a) as **cloud**.

#### 4.5.2 Urban

Flag all pixels with a ESA WorldCover raster value of 50 as **urban**, which corresponds to the "built-up" class. ESA WoldCover documentation defines the built-up class as: "Land covered by buildings, roads and other man-made structures such as railroads. Buildings include both residential and industrial building. Urban green (parks, sport facilities) is not included in this class. Waste dump deposits and extraction sites are considered as bare" (Zanaga et al., 2021). 


#### 4.5.3 Water

Flag all pixels identified as "water" by the EMIT L2A Mask product (Green, 2022a) and all pixels that intersect with GSHHG global database of coastlines and rivers (Wessel and Smith, 1996) as **water**. 

#### 4.5.4 Snow/Ice

The Normalized Difference Snow Index (NDSI) is a index-based metric used to identify pixels that are most likely to be snow and/or ice (Hall et al., 1995). The following equation is used to calculate NDSI from the EMIT L2A surface reflectance product (Green, 2022b): 

$$
NDSI = \frac{Green - SWIR}{Green + SWIR}
$$

Where 560 nm is the "Green" wavelength band and 1600 nm is the "SWIR (shortwave-infrared)" wavelength band. Following recommendations in the literature that suggest a global NDSI threshold of 0.4 (Hall et al., 2015), flag all pixels with an NDSI value greater than 0.4 as **snow/ice**.

### 4.6 Fractional Cover QA Product Output Variables
The EMIT output data products delivered to the DAAC use their formatting conventions, the system operates internally on data products stored as binary data cubes with detatched human-readable ASCII header files.

The QA product is a single band cloud-optimized GeoTIFF (COG), where each flagged QA pixel is assigned one of the following values:  
 * 1 = Cloud 
 * 2 = Urban 
 * 3 = Water
 * 4 = Snow/Ice

 For pixels that contain multiple QA flags (e.g., a water pixel covered by clouds), the following hierarchy is employed: 
 1. If the pixel contains clouds, QA = 1 
 2. If the pixel contains built-up material, QA = 2 
 3. If the pixel contains water or is a coastal pixel, QA = 3
 4. If the pixel is classified as snow/ice, QA = 4

This hierarchy order minimizes incorrect classification of pixels with NDSI thresholding, which is known to over-identify liquid water as snow/ice. 

## 5 Algorithm Usage Constraints

## 6 Performance Assessment

### 6.1 Validation Methods

### 6.2 Uncertainties

### 6.3 Known Issues

#### 6.3.1 Fractional Cover 

#### 6.3.2 Fractional Cover QA 
* The GSHHG dataset does not include major inland lakes, including the Great Lakes, which should be classified as QA water pixels
* EMIT L2A mask product may fail to identify liquid water for acquisitions with substantial sun glint, leading to misclassified water pixels in the QA product
* NDSI threhsold value of 0.4 may over-identify liquid water pixels as snow/ice pixels in the QA product

## 7 Algorithm Implementation

### 7.1 Algorithm Availability

### 7.2 Input Data Access
All input data required to run the code not found in the above repositories is available at the NASA Digital Archive. 

### 7.3 Output Data Access
All output data *will be* available at the NASA Digital Archive and searchable from NASA Earthdata Search (https://search.earthdata.nasa.gov/search). 

## 8 Significance Discussion
This ATBD describes the implementation of a three-class fractional cover classification.

## 9 Open Research
All code described here is open source.  All publications sponsored by the EMIT mission related to this code are open access. 

## 10 Acknowledgements

## 11 Contact Details

K. Dana Chadwick 

ORCID: 0000-0002-5633-4865 

Email: dana.chadwick@jpl.nasa.gov 

Role(s) related to this ATBD: writing - original and revision, methodology, software. 

Affiliation – Jet Propulsion Laboratory, California Institute of Technology 

--

Red Willow Coleman

ORCID: 0000-0001-6582-4922

Email: willow.coleman@jpl.nasa.gov 

Role(s) related to this ATBD: writing - original and revision, methodology, software. 

Affiliation – Jet Propulsion Laboratory, California Institute of Technology 


## References

* Green, R. (2022a). <i>EMIT L1B At-Sensor Calibrated Radiance and Geolocation Data 60 m V001</i> [Data set]. NASA Land Processes Distributed Active Archive Center. https://doi.org/10.5067/EMIT/EMITL1BRAD.001

* Green, R. (2022b). <i>EMIT L2A Estimated Surface Reflectance and Uncertainty and Masks 60 m V001</i> [Data set]. NASA Land Processes Distributed Active Archive Center. https://doi.org/10.5067/EMIT/EMITL2ARFL.001

* Hall,  Dorothy K. and Riggs,  George A. and Salomonson,  Vincent V. (1995). <i>Development of methods for mapping global snow cover using moderate resolution imaging spectroradiometer data</i>. Remote Sensing of Environment, 54, 0034-4257. http://dx.doi.org/10.1016/0034-4257(95)00137-P

* Hall, Dorothy K. and Riggs, George A. and Román, Miguel O. (2015). <i> VIIRS Snow Cover Algorithm Theoretical Basis Document (ATBD) </i>. https://viirsland.gsfc.nasa.gov/PDF/VIIRS_snow_cover_ATBD_2015.pdf

* Wessel, P., and W. H. F. Smith (1996), <i>A global, self-consistent, hierarchical, high-resolution shoreline database</i>. Journal of Geophysical Research, 101(B4), 8741–8743. https://doi.org/10.1029/96JB00104

* Zanaga, D., Van De Kerchove, R., De Keersmaecker, W., Souverijns, N., Brockmann, C., Quast, R., Wevers, J., Grosu, A., Paccini, A., Vergnaud, S., Cartus, O., Santoro, M., Fritz, S., Georgieva, I., Lesiv, M., Carter, S., Herold, M., Li, Linlin, Tsendbazar, N.E., Ramoino, F., Arino, O., 2021. <i>ESA WorldCover 10 m 2020 v100</i>. https://doi.org/10.5281/zenodo.5571936

