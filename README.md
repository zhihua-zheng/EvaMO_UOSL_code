[![DOI](https://zenodo.org/badge/287464756.svg)](https://zenodo.org/badge/latestdoi/287464756)

**This repository contains the MATLAB scripts used in the paper "Evaluating Monin-Obukhov scaling in the unstable oceanic surface layer"**

You can reach me at [zhihua@uw.edu](mailto:zhihua@uw.edu) for questions related to the usage of this code. 

### Description
This set of scripts fufills the analysis presented in our paper in three steps:
- Load and process original data
- Calculate key parameters
- Produce figures

### Prerequisites
This code is tested in MATLAB R2019b.

It also uses some external toolboxes/functions, as listed below:
- [Gibbs-SeaWater (GSW) Oceanographic Toolbox](http://www.teos-10.org/software.htm)
- [Air-sea toolbox](https://github.com/sea-mat/air-sea)
- [cbrewer](https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab)
- [Custom Colormap](https://www.mathworks.com/matlabcentral/fileexchange/69470-custom-colormap)
- [RGB](https://www.mathworks.com/matlabcentral/fileexchange/24497-rgb-triple-of-color-name-version-2)
- [MarkerTransparency](https://www.mathworks.com/matlabcentral/fileexchange/65194-peterrochford-markertransparency)
- [suplabel](https://www.mathworks.com/matlabcentral/fileexchange/7772-suplabel)
- [tight_subplot](https://www.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w)

### Data
The data used in this work are mostly publicly available (see the "Data availability statement" in our paper). For convience, a copy of the original data can be downloaded either from [this Google Drive link](https://drive.google.com/file/d/13UYYOT9AXFufjMw6_wr4-hoNv7M3tT7v/view?usp=sharing), or from the most recent [release](https://github.com/zhihua-zheng/EvaMO_UOSL_code/releases/tag/v1.0) of this repository.

### Work flow
1. Download this repository to your local computer. Alternatively, you can clone it using the command line in Terminal,

```bash
git clone https://github.com/zhihua-zheng/EvaMO_UOSL_code.git
```
Add this repository to your MATLAB search path

```matlab
addpath(genpath('<path to this repository>'))
```

2. Download the data file `EvaMO_UOSL_data.zip` as described in the previous [Data](#-data) section. Unzip it into the `EvaMO_UOSL_code` directory.

3. Run the main analysis script `evaMO_main.m`.

### Attribution
This code is freely available for reuse as described in the MIT License. However, if you use this code in an academic publication, a propriate citation to the source would be appreciated:

* Zhihua Zheng. (2020, August 17). Analysis scripts for "Evaluating Monin-Obukhov scaling in the unstable oceanic surface layer" (Version v1.0). Zenodo. doi:10.5281/zenodo.3988503
