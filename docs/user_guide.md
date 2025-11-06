# Earth Surface Mineral Dust Source Investigation (EMIT)

## Level 2B Fractional Cover Data Product User Guide

### Draft v01

### 1	Introduction

#### 1.1	Identification
This document describes information about the file structure and datasets provided in the EMIT L2BFRCOV data product. The algorithms and data content of the L2BFRCOV data products are described briefly in this guide, with the purpose of providing the user with sufficient information about the content and structure of the data files to enable the user to access and use the data, in addition to understanding the uncertainties involved in the products.

#### 1.2	Overview
The EMIT Project delivers space-based measurements of surface mineralogy of the Earth’s arid dust source regions. These measurements are used to initialize Earth System Models (ESM) of the dust cycle, which describe the generation, lofting, transport, and deposition of mineral dust. Earth System Models incorporate the dust cycle to estimate the impacts of mineral dust on the optical and radiative properties of the atmosphere, and a variety of environmental and ecological processes. EMIT on the ISS makes measurements over the sunlit Earth’s surface in the range of ±52° latitude. EMIT-based maps of the fractional cover of surface classes is an essential product needed for analysis of the relative abundance of source minerals to address the prime mission science questions, as well as supporting additional science and applications uses.

The EMIT instrument is a Dyson imaging spectrometer that uses contiguous spectroscopic measurements in the visible to short wavelength infrared region of the spectrum to resolve absorption features of dust-forming minerals. From the instrument’s focal plane array, on-board avionics reads out raw detector counts at 1.6 Gbps, then digitizes and stores this data to a high-speed Solid-State Data Recorder (SSDR).  From there, the avionics software reads the raw uncompressed data, packages this data into frames of 32 instrument lines, screens for cloudy pixels within the frames, and performs a lossless 4:1 compression of the frame’s science data before storing the processed, compressed data back onto the SSDR.  The data is later read from the SSDR, wrapped in CCSDS packets and then formatted as ethernet packets for transmission over the International Space Station (ISS) network and downlinked to the EMIT Instrument Operation System (IOS).  Once on the ground, the EMIT IOS delivers the raw ethernet data to the SDS where Level 0 processing removes the Huntsville Operations and Support Center (HOSC) ethernet headers, groups CCSDS packet streams by APID, and sorts them by course and fine time.

The Level 2B Fractional Cover (L2BFRCOV) data product contains geolocated fractional cover estimates of surface classes derived from the Level 2A surfance reflectance data. The fractional cover estimates are derived using the EndMember Combination Monte Carlo, E(MC)<sup>2</sup>, documented in Ochoa et al. (2025), which models each pixel's spectrum as a linear combination of endmember spectra representing different surface classes. The L2BFRCOV product provides the fractional cover estimates along with associated uncertainty metrics and masking quality flags to help users assess the reliability of the data. In addition, the geolocation of all pixel centers is included as well as the calculation of observation geometry and illumination angles on a pixel-by-pixel basis. Each image line of the data product is also UTC time-tagged.

The EMIT Fractional Cover products are provided as Cloud Optimized GeoTIFF (COG) files, with quicklooks as PNG files.

#### 1.3 File Formats

##### 1.3.1    Metadata Structure

EMIT is operating from the ISS, orbiting Earth approx.16 times in a 24-hour day period. EMIT starts and stops data recording based on a surface coverage acquisition mask. The top-level metadata identifier for EMIT data is an orbit, representing a single rotation of the ISS around Earth. Within an orbit, a period of continuous data acquisition is called a scene and is identified by a scene number. An orbit contains multiple scenes, where each scene can cover up to thousands of kilometers down-track, depending on the acquisition mask map. Each scene is broken into granules of 1280 lines down-track. The last granule in a scene is merged into the one before, making the last granule to be between 1280 and 2560 lines down-track. Granules can be downloaded as COG files, and are identified by a date-time string in the COG file name.  

The EMIT Level 2B Fractional Cover collection (EMITL2BFRCOV) contains 7 COG files, and 1 quicklook PNG file (Browse). with file content described in Table 1-1.

Table 1 1: EMITL2BFRCOV collection file list and naming convention

|  File  |  Description  |
|--------|----------------|
| EMIT_L2B_FRCOVPV_&lt;VVV&gt;_&lt;YYYYMMDDTHHMMSS&gt;_&lt;OOOOO&gt;_&lt;SSS&gt;.tif  | Fractional Cover, Photosynthetic Vegetation (aerial fraction)|
| EMIT_L2B_FRCOVNPV_&lt;VVV&gt;_&lt;YYYYMMDDTHHMMSS&gt;_&lt;OOOOO&gt;_&lt;SSS&gt;.tif  | Fractional Cover, Non-Photosynthetic Vegetation (aerial fraction) |
| EMIT_L2B_FRCOVBARE_&lt;VVV&gt;_&lt;YYYYMMDDTHHMMSS&gt;_&lt;OOOOO&gt;_&lt;SSS&gt;.tif  | Fractional Cover, Bare Soil (aerial fraction) |
| EMIT_L2B_FRCOVPVUNC_&lt;VVV&gt;_&lt;YYYYMMDDTHHMMSS&gt;_&lt;OOOOO&gt;_&lt;SSS&gt;.tif  | Fractional Cover Uncertainty, Photosynthetic Vegetation (aerial fraction) |
| EMIT_L2B_FRCOVNPVUNC_&lt;VVV&gt;_&lt;YYYYMMDDTHHMMSS&gt;_&lt;OOOOO&gt;_&lt;SSS&gt;.tif  | Fractional Cover Uncertainty, Non-Photosynthetic Vegetation (aerial fraction) |  
| EMIT_L2B_FRCOVBAREUNC_&lt;VVV&gt;_&lt;YYYYMMDDTHHMMSS&gt;_&lt;OOOOO&gt;_&lt;SSS&gt;.tif  | Fractional Cover Uncertainty, Bare Soil (aerial fraction) |
| EMIT_L2B_FRCOVQC_&lt;VVV&gt;_&lt;YYYYMMDDTHHMMSS&gt;_&lt;OOOOO&gt;_&lt;SSS&gt;.nc  | QC Flag Bands (integer)|
| EMIT_L2B_FRCOV_&lt;VVV&gt;_&lt;YYYYMMDDTHHMMSS&gt;_&lt;OOOOO&gt;_&lt;SSS&gt;.png  | Browse |



&lt;VVV&gt; gives the software version number, e.g., 001

&lt;YYYYMMDDTHHMMSS&gt; is a time stamp, e.g., 20220101T083015

&lt;OOOOO&gt; is the orbit identification number, e.g., 12345

&lt;SSS&gt; is the scene identification number, e.g., 007

##### 1.3.2    L2BFRCOV Data Products
The EMIT L2B Fractional Cover collection contains estimated surface fractional cover and quality assurance bands as Cloud Optimized GeoTIFF (COG) files, which are organized to allow efficient access and retrieval of data subsets.

#### 1.4 Product Availability
The EMIT L2BFRCOV products will be available at the NASA Land Processes Distributed Active Archive Center (LP DAAC, https://lpdaac.usgs.gov/) and through NASA Earthdata (https://earthdata.nasa.gov/).


### 2	Fractional Cover Estimation and Masking
EMIT's Level 2B Fractional Cover (L2BFRCOV) product is generated using the EndMember Combination Monte Carlo, E(MC)<sup>2</sup>, documented in Ochoa et al. (2025). This method models each pixel's spectrum as a linear combination of endmember spectra representing different surface classes. The fractional cover estimates are derived by solving for the coefficients of these endmembers that best fit the observed spectra.
The L2BFRCOV QC Flags include several flags to ensure the reliability of the fractional cover estimates. Pixel flags are provided for cloud contamination, urban cover, water, and snow/ice masks. A more detailed description of the masking procedure and quality control flags can be found in the EMIT Level 2BFRCOV ATBD in this repository.

### 3 References

* Ochoa, F., Brodrick, P. G., Okin, G. S., Ben-Dor, E., Meyer, T., Thompson, D. R., & Green, R. O. (2025). <i>Soil and vegetation cover estimation for global imaging spectroscopy using spectral mixture analysis.</i> Remote Sensing of Environment, 324, 114746.

### 4	Acronyms
| Acronym | Definition |
|---------|------------|
| ATBD    | Algorithm Theoretical Basis Document |
| APID    | Application Process Identifier |
| CCSDS   | Consultative Committee for Space Data Systems |
| COG     | Cloud Optimized GeoTIFF |
| DAAC    | Distributed Active Archive Center |
| ESM     | Earth System Models |
| EMIT    | Earth Surface Mineral Dust Source Investigation |
| E(MC)<sup>2</sup> | EndMember Combination Monte Carlo |
| HOSC    | Huntsville Operations and Support Center |
| IOS     | Instrument Operation System |
| ISS     | International Space Station |
| L2BFRCOV | Level 2B Fractional Cover |
| LP DAAC | Land Processes Distributed Active Archive Center |
| NPV     | Non-Photosynthetic Vegetation |
| PNG     | Portable Network Graphics |
| PV      | Photosynthetic Vegetation |
| QC      | Quality Control |
| SDS     | Science Data System |
| SSDR    | Solid-State Data Recorder |

