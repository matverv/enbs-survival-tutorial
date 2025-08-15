# Expected Net Benefit of Sampling for Survival Data : A Tutorial
This GitHub repository contains the R code and analysis scripts accompanying the paper:
- Vervaart, M. (2024). Calculating the Expected Net Benefit of Sampling for Survival Data: A Tutorial and Case Study. Medical Decision Making, in press.

### Updates & corrections
- August 2025: Correction to discounting and time horizon logic
  The logic for handling decision horizons was replaced with a more direct calculation, and a new `r pv_annuity` function was added to improve the accuracy of the discounting.

***

### Getting started: Running the analysis

Follow these steps to run the case study analysis included in the paper.

1.  **Clone the repository**

    Clone or download this repository to your local machine.
    ```bash
    git clone https://github.com/matverv/enbs-survival-tutorial.git
    ```

2.  **Open the RStudio project**

    Navigate into the cloned folder and open the **`enbs-survival-tutorial.Rproj`** file. This will launch RStudio and automatically set the correct working directory.

3.  **Run an analysis script**

    In RStudio, open any of the scripts from the `analysis/` folder and run them. The scripts will automatically install any required packages.

    - [`S1_main_paper.R`: EVPPI, EVSI and ENBS for OS - Ongoing Trial](https://github.com/matverv/enbs-survival-tutorial/blob/main/analysis/S1_main_paper.R)
    - [`S2_os_and_pfs.R`: EVPPI, EVSI and ENBS for OS and PFS - Ongoing Trial](https://github.com/matverv/enbs-survival-tutorial/blob/main/analysis/S2_os_and_pfs.R)
    - [`S3_new_trial.R`: EVPPI, EVSI and ENBS for OS and PFS - New Trial](https://github.com/matverv/enbs-survival-tutorial/blob/main/analysis/S3_new_trial.R)
    - [`S4_reconstruct_ipd.R`: Reconstruction of Individual Patient Data](https://github.com/matverv/enbs-survival-tutorial/blob/main/analysis/S4_reconstruct_ipdR.R)

***


## Using your own data

To apply the analysis to your own case study, you will need to undertake the following steps:

-   Replace the case-specific probabilistic analysis (PA) function `run_pa_pembro` and survival data `data/ipd_os.RData` and if relevant `data/pfs_os.RData` with your own PA function and survival datasets, ensuring that your PA results and survival data conform to the formats expected by the analysis scripts. Detailed instructions on the required data formats and analysis settings can be found in the main paper or by examining the case study provided in the repository. If you do not have access to (reconstructed) IPD, you can still implement the analysis by utilizing the reported follow-up time in the trial and making assumptions about trial dropout as explained in the main paper.

-   Review and adjust any analysis settings such as the number of trial dropouts, incidence, prevalence, trial costs, and other case-specific details to accurately reflect your own case study.

### For Excel users

If your cost-effectiveness model and PA are conducted in Excel and you wish to use these results with our R scripts, follow these additional steps:

First, ensure your CSV files are structured correctly. The rows and columns must match the required format for costs, QALYs, and survival probabilities.

| Costs (`m_c.csv`) | QALYs (`m_e.csv`) | Survival Probabilities (`m_os1.csv`) |
| :---: | :---: | :---: |
| ![Costs Example](https://github.com/matverv/enbs-survival-tutorial/assets/58030182/3fe35bd2-45a7-4e20-b766-8187a0fdc98a) | ![QALYs Example](https://github.com/matverv/enbs-survival-tutorial/assets/58030182/cb0b37ed-6dab-4ed9-8cb2-299e9e52d36e) | ![Survival Probabilities Example](https://github.com/matverv/enbs-survival-tutorial/assets/58030182/48d82676-b2cb-4e8c-a3b8-365f13276e19) |

Next, use the following R code in your script to load the data. For best results, place your CSVs in the `data/` folder.

```r
# Load PA results from CSV files
m_c    <- read.csv("data/your_costs.csv")
m_e    <- read.csv("data/your_qalys.csv")
m_os1  <- read.csv("data/your_os1.csv")
m_os2  <- read.csv("data/your_os2.csv")
m_pfs1 <- read.csv("data/your_pfs1.csv")
m_pfs2 <- read.csv("data/your_pfs2.csv")

# Structure the data into the required list format
l_surv <- list(m_os1 = m_os1, m_os2 = m_os2, m_pfs1 = m_pfs1, m_pfs2 = m_pfs2)
l_pa   <- list(m_c = m_c, m_e = m_e, l_surv = l_surv)
```

***


### Repository structure

- `R/`: Contains the function files used by the analysis scripts.
- `data/`: Contains the case study data.
- `analysis/`: Contains the main analysis scripts.

***

## Citation

Please cite the following when using this code:

  Vervaart, M. (2024). Calculating the Expected Net Benefit of Sampling for Survival Data: A Tutorial and Case Study. Medical Decision Making, in press.
  