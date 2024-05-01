# Expected Net Benefit of Sampling for Survival Data : A Tutorial
This GitHub repository contains the R code and analysis scripts accompanying the paper:
- Calculating the Expected Net Benefit of Sampling for Survival Data: A Tutorial and Case Study. 2024.

## Repository Structure

- `R/`: Contains the function files used by the analysis scripts.
- `data/`: Contains the reconstructed IPD for OS and PFS, general population mortality data and postscript files of the KM plots.
- `analysis/`: Contains the analysis scripts referenced in the paper:
   - [Appendix A: EVPPI, EVSI and ENBS for OS - Ongoing Trial](https://github.com/matverv/enbs-survival-tutorial/blob/main/analysis/S1_main_paper.R)
   - [Appendix B: EVPPI, EVSI and ENBS for OS and PFS - Ongoing Trial](https://github.com/matverv/enbs-survival-tutorial/blob/main/analysis/S2_os_and_pfs.R)
   - [Appendix C: EVPPI, EVSI and ENBS for OS and PFS - New Trial](https://github.com/matverv/enbs-survival-tutorial/blob/main/analysis/S3_new_trial.R)
   - [Appendix D: Reconstruction of Individual Patient Data for OS and PFS](https://github.com/matverv/enbs-survival-tutorial/blob/main/analysis/S4_reconstruct_ipdR.R)

## How to Use

1. Clone or download this repository to your local machine:

```bash
git clone https://github.com/matverv/enbs-survival-tutorial.git
```
  or download as a ZIP file and extract it.

2. Install R and RStudio, if not already installed:
- R: https://cran.r-project.org/
- RStudio: https://www.rstudio.com/products/rstudio/download/

3. Open RStudio and set your working directory to the cloned or downloaded repository folder:

```r
setwd("path/to/enbs-survival-tutorial")
```

4. Open and run the desired analysis script from the analysis/ folder directly in RStudio. The scripts will automatically install any required packages that are not already installed on your system and then load them for use.
 

## Citation

Please cite the following when using this code:

    Calculating the Expected Net Benefit of Sampling for Survival Data: A Tutorial and Case Study. (2024). [Journal Information].

