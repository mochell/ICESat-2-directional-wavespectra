# Code for Directional Surface Wave Spectra And Sea Ice Structure from ICEsat-2 Altimetry

Momme Hell,
July 2022,
Brown University

This is a minimal Working example for „“. We provide the original ATL03 and ATL07 data from XX downloaded on:
The in Hell and Horvat described analysis can be reproduced for the example tracks following the guide below

This package uses and modifies code from icesat2_toolkit
https://read-icesat-2.readthedocs.io/
https://github.com/tsutterley/read-ICESat-2
(One might need to follow their authentication steps to download data ...)

## python environment:
The python environment is provided in
```
2021_ICESat2_tracks_py37.yml
```
Instructions of how to install can be found here:
https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html

## setting up environment

After installing and activating the python environment set the base path at the following places:
#### in config/2021_IceSAT2_startup.py
base_path = "your/path/to/base/"

#### in config/config.json
replace "/Users/Shared/Projects/2021_IceSAT2_tracks/release" with "your/path/to/base/"

specifically for the following variables in config.json:
```
"work": "/Users/Shared/Projects/2021_IceSAT2_tracks/release/data/work/",
"scratch": "/Users/Shared/Projects/2021_IceSAT2_tracks/release/data/scratch/",
"plot": "/Users/Shared/Projects/2021_IceSAT2_tracks/release/plots/",
"processed": "/Users/Shared/Projects/2021_IceSAT2_tracks/release/data/processed/",
"local_script": "/Users/Shared/Projects/2021_IceSAT2_tracks/release/modules/",
"config": "/Users/Shared/Projects/2021_IceSAT2_tracks/release/config/",
"base": "/Users/Shared/Projects/2021_IceSAT2_tracks/release/"
```

## Run code:
The provided code is a simplified from a processing structure Maintained via make files. Here we use simple manual execution.
The easiest way to run the code is through command lines. The processing steps are numbered from A01, A02, .. to B01, B03, .. The A section refers to downloading and pre-handling the data and the B section to the part where the algorithm(s) are applied.

### A01b) Find ATL03 and create case ID
This code find the ATL03 tracks to the corresponding ATL07/10 tracks. Here the ATL03 and ATL07 tracks are already downloaded, but to check one can run in the analysis_db/ folder:
```
python3 A01b_ALT07_SHNH_variance_tester.py 20190224012038_08800201_005_01 SH_publish True
```
(uses multicore processing if available)
This also creates the ID_files in work/SH_publish/A01b_ID/. These files are needed to run the rest of the processing. They constrain all needed references and parameters for each case.

The data for the following track IDs are provided, such that the above command can be applied to all the IDs listed below:

```
20190502040726_05180301_005_01
20190219063727_08070201_005_01
20190224012038_08800201_005_01
20190502005851_05160301_005_01
```

Each of these track IDs create two IDs one for the dominantly descending or ascending part of the track. The manuscript mainly uses the following three IDs:
```
SH_20190224_08800210
SH_20190219_08070210
SH_20190502_05160312
```

while additional cases can also be created, since they are rooting from the same ATL03 data.

```
SH_20190219_08070212
SH_20190224_08800212
SH_20190502_05160310
SH_20190502_05180310
SH_20190502_05180312
```

## Re-create the database
The plots will occur at plots/SH/SH_publish/

### A01c) merging the ATL03 data file to one new file and apply corrections:

```
python3 A01c_mergefiles.py SH_20190224_08800210 SH_publish True
```
The output is saved in data/scratch/SH_publish/A01c_ ...

### A02b) deriving WW3 priors:

```
python3 A02b_WW3_hindcast_prior.py SH_20190224_08800210 SH_publish True
```
The output is saved in data/scratch/SH_publish/A01c_ ...


### B01) rebinning the data in 20 meter stencils and apply local reference systems
(uses multicore processing if available)
```
python3 B01_filter_regrid_segments.py SH_20190224_08800210 SH_publish True
```
The output is saved at work/SH_publish/B01_regrid/



### B02) This applies the GFT to all segments and all beams.
(uses multicore processing if available)
```
python3 B02_make_spectra_gFT.py SH_20190224_08800210 SH_publish True
```
output is at work/SH_publish/B02_spectra/



### B03) Plot results of the GFT
```
python3 B03_plot_spectra_ov.py SH_20190224_08800210 SH_publish True
```
plots are found in the respective plot folder


### B04) Angle inversion using MCMC sampling.
(uses multicore processing if available)

requires results from B03 and A02b
```
python3 B04_angle.py SH_20190224_08800210 SH_publish True
```
output is at work/SH_publish/B04_angle/
plots are found in the respective plot folder


### B05) Define angle based on the B04 result
(uses multicore processing if available)
```
python3 B05_define_angle.py SH_20190224_08800210 SH_publish True
```
output is at work/SH_publish/B04_angle/
plots are found in the respective plot folder


### B06) Find cur off frequency, create the corrected 2D-wavenumber spectra, reconstruct the surface heights and store final product in files that contain all data
```
python3 B06_correct_separate_var.py SH_20190224_08800210 SH_publish True
```
output is at work/SH_publish/B04_angle/
plots are found in the respective plot folder

These files constain:
- The original re-binned photon data
- The modeled wave heights in real or spectral space including their uncertainties.
- The variance decomposition in real space
- The spectral estimates with angle corrected coordinates
- The "best guess" 2D spectra along the track.
