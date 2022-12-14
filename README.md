# urea_physchem
This repository contains raw spectral data and analysis script on optical properties of carbon nanoparticles obtained from urea or thiourea and citric acid. The major results of the analysis are collected in a manuscript submitted to [Physchem](https://www.mdpi.com/journal/physchem). The paper link will be added here upon the manuscript acceptance.  

## Files
* `composition.R` General script used for data import, data processing, and plots generation. The script is self-explanatory, and the clue points are described in the comments
* `data/QY_oct_22/abs` Absorbance spectra of the samples. The filenames code the samples as follows: `X_ca_Y_Z_M.txt`, where **X** is `t` for thiourea and `u` for urea, **Y** and **Z** reflect the mass ratio of (thio)urea to citric acid in the sample during preparation, and **M** is `d` for the dialysate (low-molecular fraction) and `p` for purified (colloidal fraction). `ref.txt` file containes the reference (quinine sulfate) spectrum.
* `data/QY_oct_22/fluoro` Emission spectra for the same set of samples. The filenames are consistent with the absorbance spectra.
* `data/cnd_vario_Oct22_X.xlsx` Excitation and emission spectra of selected samples at variable excitation/emission wavelengths obtained using a microplate reader. **X** is the plate number. The samples decription is in the `samples` sheet of the Excel files.
* `urea_3_1_distribution.pdf` and `urea_3_1_distribution_clean.txt` Instrument report on the particles size distribution (NTA data) in a selected sample and the data used for generation of the plot in the manuscript.

## Disclaimer on data integrity
The `composition.R` script is supposed to be slightly refactored shortly, yet keeping the analysis track. Further analysis (if any) will be performed in a separate script.  
As of Dec 12, 2022, the plots generated by the `composition.R` script are just identical to these in the submitted manuscript, however it is not guaranteed that the plots appearance is kept during refactoring of the script.
