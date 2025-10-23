# EMIT Fractional Cover Product Delivery - Algorithm Theoretical Basis Document (ATBD)

*Philip Brodrick*<sup>1</sup>, *Francisco Ochoa*<sup>1,2</sup>, *Greg Okin*<sup>2</sup>, *Red Willow Coleman*<sup>1,2</sup>, *K.D. Chadwick*<sup>1</sup>

<sup>1</sup>Jet Propulsion Laboratory, California Institute of Technology
<sup>2</sup>University of California, Los Angeles

Corresponding author: Philip Brodrick (philip.brodrick@jpl.nasa.gov)

**Key Points:**
   1. Please note that this ATBD will be updated on an ongoing basis as the EMIT extended mission progresses. This is intended to be a place where the community can find the most up-to-date information on the current plans for algorithm development and offer contributions.
   2. This is a three component model with a series of quality assessment (QA) flags. The mission will be producing a more detailed V2 of this product which will supersede this version during FY26. 
   3. If you identify issues with this product not current outlined in the Known Issues section, please contribute that information here to help the community.

**Version:** 1.0

**Release Date:** TBD

**DOI:** TBD

## Abstract

## Plain Language Summary

The terrestrial surface of the Earth is comprised of many substances. Using imaging spectroscopy, we work to classify these substances into groups in order to help with the interpretation of various algorithms. In this Version 1 of the EMIT Fractional Cover Product, we are focused on being able to determine the fractions comprised of bare rock/soil (areas where mineral determination algorithms are appropriate for interpretation), photosynthetic vegetation (largely identified by their characteristic red edge and low reflectance in the SWIR, appropriate for interpretation of trait estimation algorithms), and non-photosynthetic vegetation (which can include leaf litter, wood, some lichens and biocrusts). All other surface types, as well as areas that are obscured by clouds, we seek to flag using the QA bands associated with this product.

We are providing this Algorithm Theoretical Basis Document in a github markdown format in order to provide a record of ongoing updates as algorithms improve via the commit record, as well as allowing the community to engage more directly in the process of documentation in keeping with NASA's commitment to open science. 

### Keywords: hyperspectral imaging, imaging spectroscopy, fractional cover, terrestrial

## 1 Version Description

This is Version 1.0 of the EMIT fractional cover and fractional cover quality assessment (QA) product.

## 2 Introduction
Mineral dust aerosols originate as soil particles lifted into the atmosphere by wind erosion. Mineral dust created by human activity makes a large contribution to the uncertainty of direct radiative forcing (RF) by anthropogenic aerosols (USGCRP and IPCC). Mineral dust is a prominent aerosol constituent around the globe. However, we have a poor understanding of its direct radiative effect, partly due to uncertainties in the dust mineral composition. Dust radiative forcing is highly dependent on its mineral-specific absorption properties. The current range of iron oxide abundance in dust source models translates into a large range of values, even changing the sign of the forcing (-0.15 to 0.21 W/m2) predicted by Earth System Models (ESMs) (Li et al., 2020). The National Aeronautics and Space Administration (NASA) selected the Earth Surface Mineral Dust Source Investigation (EMIT) to close this knowledge gap. EMIT is an instrument on the International Space Station (ISS) that can directly measure and map the soil mineral composition of critical dust-forming regions worldwide.

## 3 Context/Background
As part of the EMIT mission, a three-class fractional cover product is being developed to provide estimates of the fractional cover of photosynthetic vegetation (PV or GV), non-photosynthetic vegetation (NPV), and bare soil (S) within each 60 m EMIT pixel. This product will be used to help interpret the surface mineralogy results from the EMIT mission, as well as provide a valuable dataset for the broader science community.

Example spectra from the three classes are shown in Figure 3 from Ochoa et al, 2025. 

![Example spectra from the three classes](figs/Ochoa_spec.jpg "Example spectra from the three classes")

*Figure 3-1. A random selection of spectra from the endmember library used for unmixing. Figure from Ochoa et al., 2025.* 

### 3.1 Historical Perspective
Fractional cover products have been used for a variety of applications, including land cover classification, vegetation monitoring, and soil moisture estimation. The EMIT fractional cover product builds on this previous work by using the unique capabilities of the EMIT instrument to provide high-resolution estimates of fractional cover across the globe as the purview of the mission expands in its extended mission activities. The EMIT fractional cover product was developed as part of the prime mission activities but was not formally delivered to the DAAC. The product is being refined and formally delivered to the DAAC as part of the EMIT extended mission activities.
### 3.2 Additional Information
The fractional cover product described here is a three component model which is masked by a series of quality assessment (QA) flags. The three components are photosynthetic vegetation (PV), non-photosynthetic vegetation (NPV), and bare soil (S). The QA flags include cloud, urban, water, and snow/ice. The fractional cover product is derived using a Monte Carlo Spectral Unmixing approach, which is described in detail in Section 4. The QA flags are derived using a combination of the EMIT L2A surface reflectance product (Green, 2022b), the ESA WorldCover land cover product (Zanaga et al., 2021), and the GSHHG global database of coastlines and rivers (Wessel and Smith, 1996). The QA flags are described in detail in Section 4.5.
## 4 Algorithm Description

### 4.1 Scientific Theory
The fractional cover product is derived using a Monte Carlo Spectral Unmixing approach. This approach uses a spectral library of endmembers to estimate the fractional cover of each component within each pixel. The spectral library used for the EMIT fractional cover product is derived from field and laboratory measurements of soils and vegetation. The spectral library includes endmembers for photosynthetic vegetation, non-photosynthetic vegetation, and bare soil. The spectral library is described in detail in Ochoa et al, 2025. 

### 4.2 Mathematical Theory
To estimate fractional cover, we use a Monte Carlo Spectral Unmixing strategy, based on decades of literature (e.g., Roberts et al. 1998, Asner and Lobell 2000, and Dennison et al., 2019). Several key parameters, including the endmember selection strategy, observation normalization techniques, and the number of bootstrap samples were investigated. Simulation experiments comparing over one million synthetic spectra constructed with endmember holdout sets were utilized to select parameter values (generally following the approach of Okin et al., 2015). Selected parameter values are shown in Table 4.2.1-1. Parameters were chosen based on a combination of mean squared error, prediction variance, prediction bias, and computation time. All values can be tested through parameter selection in the unmix.jl script provided in the SpectralUnmixing EMIT SDS repository (https://github.com/emit-sds/SpectralUnmixing). A sample comparison between
two scenarios is shown in Figure 4.2.-1.

**Table 4.2-1.** _Unmixing parameter value selection_
| Parameter Name | Tested Values | Selected Values |
| --- | --- | --- |
| Endmember Selection | MESMA, Monte Carlo SMA | Monte Carlo SMA |
| Endmember Selection Strategy | Class-even, Random | Class-even |
| Normalization | None, Brightness, 1070, 1500, 1756, 2030 | Brightness |
| Maximum Number of Combinations (MESMA only) | 10, 100, 500, 1000 | NA |
| Number of Endmembers per Class per Bootstrap Draw (SMA only) | 5, 10, 30, 50 | 10 |
| Number of Bootstrap Draws | 1, 5, 10, 20, 50, 100, 200 | 50 |



![Alt text for the figure](figs/Comparisons.png "Comparisons")
*Figure 4.2-1: Comparison between two parameter combinations. On top is a Monte Carlo unmixing analysis with 10 endmembers per class per bootstrap draw, brightness normalization, class-even endmember selection, and 50 bootstrap draws. On the bottom is a MESMA-style unmixing analysis with a maximum of 1000 class-even selected endmembers, brightness normalization, and 50 bootstrap draws.*

### 4.3 Fractional Cover Algorithm Input Variables

The required input files for fractional cover production are in Table 4.3-1.

**Table 4.3-1.** _Input variables_
| Name | Description | Unit | Required |
| --- | --- | --- | --- |
| HDRF (Surface Reflectance) | hemispherical-directional reflectance factor per wavelength | unitless | true |
| Endmember Spectral Library | collection of spectral surface endmembers (HDRF per wavelength) | unitless | true |
| Observation Geometry | solar zenith angle, view zenith angle, relative azimuth angle | degree | false |

### 4.4 Fractional Cover Algorithm Output Variables
The EMIT output data products delivered to the DAAC use their formatting conventions, the system operates internally on data products stored as binary data cubes with detached human-readable ASCII header files. For the fraction cover product, the output variables are: 
1. Fractional cover, provided as an n x c x 3 BIL interleave data cube, with c columns and n lines. Each channel contains the fractional cover as calculated by Monte Carlo SMA (see section 4.2.1).
2. Fractional cover uncertainty, provided as an n x c x 3 BIL interleave data cube, with c columns and n lines. Each channel contains the estimated uncertainty of the fraction cover, as defined in section 6.2.

These products are consistent with the auxiliary data products described in the EMIT L3ASA ATBD, section 4.4.2 (Brodrick et al., 2023).

### 4.5 Fractional Cover QA Product Input Variables
The required input files for fractional cover QA production are in Table 4.5-1.

**Table 4.5-1.** _QA input variables_
| Name | Description | Spatial Resolution | Required | Type |
| --- | --- | --- | --- | --- |
| EMIT L2A Surface Reflectance | hemispherical-directional reflectance factor (HDRF) per wavelength | 60 m | true | raster
| EMIT L2A Mask | estimated surface reflectance uncertainty and masks | 60 m | true | raster
| ESA WorldCover | global landcover classification derived from Sentinel-1 and Sentinel-2 data | 10 m | true | raster
| GSHHG | the Global Self-consistent, Hierarchical, High-resolution Geography Database for coastlines and rivers | <60 m | true | vector

#### 4.5.1 Clouds

Flag all pixels identified as "cirrus" or "cloud" by the EMIT L2A Mask product (Green, 2022a) as **cloud**.

#### 4.5.2 Urban

Flag all pixels with an ESA WorldCover raster value of 50 as **urban**, which corresponds to the "built-up" class. ESA WorldCover documentation defines the built-up class as: "Land covered by buildings, roads and other man-made structures such as railroads. Buildings include both residential and industrial building. Urban green (parks, sports facilities) is not included in this class. Waste dump deposits and extraction sites are considered as bare" (Zanaga et al., 2021). 


#### 4.5.3 Water

Flag all pixels identified as "water" by the EMIT L2A Mask product (Green, 2022a) and all pixels that intersect with the GSHHG global database of coastlines and rivers (Wessel and Smith, 1996) as **water**. 

#### 4.5.4 Snow/Ice

The Normalized Difference Snow Index (NDSI) is an index-based metric used to identify pixels that are most likely to be snow and/or ice (Hall et al., 1995). The following equation is used to calculate NDSI from the EMIT L2A surface reflectance product (Green, 2022b): 

$$
NDSI = \frac{Green - SWIR}{Green + SWIR}
$$

Where 560 nm is the "Green" wavelength band and 1600 nm is the "SWIR (shortwave-infrared)" wavelength band. Following recommendations in the literature that suggest a global NDSI threshold of 0.4 (Hall et al., 2015), flag all pixels with an NDSI value greater than 0.4 as **snow/ice**.

### 4.6 Fractional Cover QA Product Output Variables
The EMIT output data products delivered to the DAAC use their formatting conventions, the system operates internally on data products stored as binary data cubes with detached human-readable ASCII header files.

The QA product is a single band cloud-optimized GeoTIFF (COG), where each flagged QA pixel is assigned one of the following values with colors associated with figures below in parentheses for reference:  
 * 1 = Cloud (orange)
 * 2 = Urban (green)
 * 3 = Water (red)
 * 4 = Snow/Ice

 For pixels that contain multiple QA flags (e.g., a water pixel covered by clouds), the following hierarchy is employed: 
 1. If the pixel contains clouds, QA = 1 
 2. If the pixel contains built-up material, QA = 2 
 3. If the pixel contains water or is a coastal pixel, QA = 3
 4. If the pixel is classified as snow/ice, QA = 4

This hierarchy order minimizes incorrect classification of pixels with NDSI thresholding, which is known to over-identify liquid water as snow/ice. 

## 5 Algorithm Usage Constraints

Note that the QA flags are not mutually exclusive and also the QA flags are not intended as a comprehensive cloud, urban, water, or snow/ice identification product. Rather, it is intended to flag pixels that are likely to be misclassified by the fractional cover algorithm and it is not recommended to use the QA product as a standalone cloud, urban, water, or snow/ice product.

## 6 Performance Assessment

### 6.1 Validation Methods

<i>[validation methodologies should reference Francisco's work here, but I don't have detailed descriptions - would recommend showing validation plots from Ochoa paper?] </i>

### 6.2 Uncertainties

To estimate the uncertainty of the Monte Carlo SMA results, we run 50 Monte Carlo simulations. During each simulation, the endmember selection is seeded differently (representing model error) and the reflectance is perturbed by a (per-wavelength) random deviation proportionate to the channelized reflectance uncertainty provided by L2A. The standard deviation of the soil fractional cover from the different simulations (𝜎 <_>𝑠 𝑗) is then used as the uncertainty.

We neglect the uncertainty associated with QA assessment for the purposes of providing a mask for this product.

### 6.3 Known Issues
The spectral library used for unmixing may not be representative of all surface types globally, leading to potential inaccuracies in fractional cover estimates in some regions. Here we describe known issues with both the fractional cover and fractional cover QA products for users to be aware of. Many of these issues are being addressed in ongoing work to improve the products and will be resolved in a planned v2 release.

#### 6.3.1 Fractional Cover 
The following issues have been identified with the fractional cover product in areas that are not flagged by the current QA product. This section will be updated as issues with the QA product are addressed or the fractional cover product is improved.
* At cloud edges and in cloud shadows regions are often inaccurately estimated as high NPV (Figures 1 & 2). Users may opt to use the cloud buffer mask provided in the EMIT L2A Mask product (Green, 2022a) to remove these pixels from analysis.

![Example of NPV misclassification](figs/NPV_KnownIssue1.png "Example of NPV misclassification")

*Figure 1: Example of NPV misclassification along cloud edges and in cloud shadows on the north side of clouds present in the image.*

![Example of NPV misclassification](figs/NPV_KnownIssue2.png "Example of NPV misclassification")

*Figure 2: Example of NPV misclassification along cloud edges on the north side of clouds present in the image.*


#### 6.3.2 Fractional Cover QA 
* The GSHHG dataset does not include major inland lakes, including the Great Lakes, which should be classified as QA water pixels
* EMIT L2A mask product may fail to identify liquid water for acquisitions with substantial sun glint, leading to misclassified water pixels in the QA product
* NDSI threshold value of 0.4 may over-identify liquid water pixels as snow/ice pixels in the QA product
* Shaded snow/ice pixels may be identified as water in the QA product.
* Bright snow/ice pixels may be misidentified as clouds in the QA product.
* Some clouds over water may result in no QA flag being assigned to the pixel.

## 7 Algorithm Implementation

### 7.1 Algorithm Availability

Fractional cover algorithms used to generate this product are available in the EMIT SDS project on github (https://github.com/emit-sds/emit-sds-frcov). 

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

Role(s) related to this ATBD: writing - original and revision, methodology, quality assessment. 

Affiliation – Jet Propulsion Laboratory, California Institute of Technology 

--

Red Willow Coleman

ORCID: 0000-0001-6582-4922

Email: willow.coleman@jpl.nasa.gov 

Role(s) related to this ATBD: writing - original and revision, methodology, software. 

Affiliation – Jet Propulsion Laboratory, California Institute of Technology 


## References

* Asner, Gregory P., and David B. Lobell. <i>A biogeophysical approach for automated SWIR unmixing of soils and vegetation.</i> Remote sensing of environment 74.1 (2000): 99-112.

* Dennison, P.E., Qi, Y., Meerdink, S.K., Kokaly, R.F., Thompson, D.R., Daughtry, C.S., Quemada, M., Roberts, D.A., Gader, P.D., Wetherley, E.B. and Numata, I., 2019. <i>Comparison of Methods for Modeling Fractional Cover Using Simulated Satellite
Hyperspectral Imager Spectra.</i> Remote Sensing, 11(18), p.2072.

* Green, R. (2022a). <i>EMIT L1B At-Sensor Calibrated Radiance and Geolocation Data 60 m V001</i> [Data set]. NASA Land Processes Distributed Active Archive Center. https://doi.org/10.5067/EMIT/EMITL1BRAD.001

* Green, R. (2022b). <i>EMIT L2A Estimated Surface Reflectance and Uncertainty and Masks 60 m V001</i> [Data set]. NASA Land Processes Distributed Active Archive Center. https://doi.org/10.5067/EMIT/EMITL2ARFL.001

* Hall,  Dorothy K. and Riggs,  George A. and Salomonson,  Vincent V. (1995). <i>Development of methods for mapping global snow cover using moderate resolution imaging spectroradiometer data</i>. Remote Sensing of Environment, 54, 0034-4257. http://dx.doi.org/10.1016/0034-4257(95)00137-P

* Hall, Dorothy K. and Riggs, George A. and Román, Miguel O. (2015). <i> VIIRS Snow Cover Algorithm Theoretical Basis Document (ATBD) </i>. https://viirsland.gsfc.nasa.gov/PDF/VIIRS_snow_cover_ATBD_2015.pdf

* Ochoa, F., Brodrick, P. G., Okin, G. S., Ben-Dor, E., Meyer, T., Thompson, D. R., & Green, R. O. (2025). <i>Soil and vegetation cover estimation for global imaging spectroscopy using spectral mixture analysis.</i> Remote Sensing of Environment, 324, 114746.

* Roberts, D., Gardner, M., Church, R., Ustin, S., Scheer, G., & Green, R. 1998, Remote Sensing of Environment, 65, 267.

* Wessel, P., and W. H. F. Smith (1996), <i>A global, self-consistent, hierarchical, high-resolution shoreline database</i>. Journal of Geophysical Research, 101(B4), 8741–8743. https://doi.org/10.1029/96JB00104

* Zanaga, D., Van De Kerchove, R., De Keersmaecker, W., Souverijns, N., Brockmann, C., Quast, R., Wevers, J., Grosu, A., Paccini, A., Vergnaud, S., Cartus, O., Santoro, M., Fritz, S., Georgieva, I., Lesiv, M., Carter, S., Herold, M., Li, Linlin, Tsendbazar, N.E., Ramoino, F., Arino, O., 2021. <i>ESA WorldCover 10 m 2020 v100</i>. https://doi.org/10.5281/zenodo.5571936

