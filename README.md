# Expected Net Benefit of Sampling for Survival Data : A Tutorial
This GitHub repository contains the R code and analysis scripts accompanying the paper:
- Vervaart et al. Calculating the Expected Net Benefit of Sampling for Survival Data: A Tutorial and Case Study. 2024.

## Repository Structure

- `R/`: Contains the function files used by the analysis scripts.
- `data/`: Contains the reconstructed IPD for OS and PFS, general population mortality data and postscript files of the KM plots.
- `analysis/`: Contains the analysis scripts referenced in the paper:
   - [Supp. 1: EVPPI, EVSI and ENBS for OS - Ongoing Trial](https://github.com/matverv/enbs-survival-tutorial/blob/main/analysis/S1_main_paper.R)
   - [Supp. 2: EVPPI, EVSI and ENBS for OS and PFS - Ongoing Trial](https://github.com/matverv/enbs-survival-tutorial/blob/main/analysis/S2_os_and_pfs.R)
   - [Supp. 3: EVPPI, EVSI and ENBS for OS and PFS - New Trial](https://github.com/matverv/enbs-survival-tutorial/blob/main/analysis/S3_new_trial.R)
   - [Supp. 4: Reconstruction of Individual Patient Data for OS and PFS](https://github.com/matverv/enbs-survival-tutorial/blob/main/analysis/S4_reconstruct_ipdR.R)

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

5. To apply the analysis to your own case study, you will need to undertake the following steps:

    - Replace the case-specific probabilistic analysis (PA) function `run_pa_pembro` and survival data `data/ipd_os.RData` and if relevant `data/pfs_os.RData` with your own PA function and survival datasets, ensuring that your PA results and survival data conform to the formats expected by the analysis scripts. Detailed instructions on the required data formats and analysis settings can be found in the main paper or by examining the case study provided in the repository. If you do not have access to (reconstructed) IPD, you can still implement the analysis by utilizing the median follow-up time in the trial and making assumptions about trial dropout as explained in the main paper.

    - Review and adjust any analysis settings such as the number of trial dropouts, incidence, prevalence, trial costs, and other case-specific details to accurately reflect your own case study.

    **For Excel Users:**

    If your cost-effectiveness model and PA are conducted in Excel and you wish to use these results with our R scripts, follow these additional steps:

      - Store your PA results as seperate CSV files, ensuring that the data structure for costs (m_c), QALYs (m_e), and survival probabilities (e.g. m_os1, m_os2, m_pfs1, m_pfs2) conforms to the expected format. Below are illustrative screenshots showing the format for `m_c` (costs), `m_e` (QALYs), and an example of survival probabilities for one treatment (`m_os1`). Ensure your data follows a similar structure. 
       <p float="left">
        <img src="https://github.com/matverv/enbs-survival-tutorial/assets/58030182/3fe35bd2-45a7-4e20-b766-8187a0fdc98a" alt="Costs Example"/>
        <img src="https://github.com/matverv/enbs-survival-tutorial/assets/58030182/cb0b37ed-6dab-4ed9-8cb2-299e9e52d36e" alt="QALYs Example"/> 
        <img src="https://github.com/matverv/enbs-survival-tutorial/assets/58030182/48d82676-b2cb-4e8c-a3b8-365f13276e19" alt="Survival Probabilities Example"/>
      </p>

      - Import the CSV files into R using the following commands:

      ```r
      m_c    <- read.csv("path/to/your/costs.csv", header = TRUE) # rows index the PA simulations and columns index the treatment strategies
      m_e    <- read.csv("path/to/your/qalys.csv", header = TRUE) # rows index the PA simulations and columns index the treatment strategies
      m_os1  <- read.csv("path/to/your/os1.csv", header = TRUE)   # rows index the model cycles and columns index the PA simulations
      m_os2  <- read.csv("path/to/your/os2.csv", header = TRUE)   # rows index the model cycles and columns index the PA simulations
      m_pfs1 <- read.csv("path/to/your/pfs1.csv", header = TRUE)  # rows index the model cycles and columns index the PA simulations
      m_pfs2 <- read.csv("path/to/your/pfs2.csv", header = TRUE)  # rows index the model cycles and columns index the PA simulations

      l_surv <- list(m_os1 = m_os1, m_os2 = m_os2, m_pfs1 = m_pfs1, m_pfs2 = m_pfs2)
      l_pa   <- list(m_c = m_c, m_e = m_e, l_surv = l_surv)
      ```


## Citation

Please cite the following when using this code:

   Vervaart et al. Calculating the Expected Net Benefit of Sampling for Survival Data: A Tutorial and Case Study. (2024). [Journal Name], [Volume(Issue)], [Page Range]. [DOI]

