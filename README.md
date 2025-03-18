# Environmental shark extinction drivers
This is the repo for the current project from the Pimiento lab under the lead of G. Mathes, where we try to estimate the true (unbiased) effect of various environmental parameters on shark extinction risk using a newly compiled database of fossil shark occurrences. 

## Related manuscript
**Climate change drives neoselachian extinctions over geological timescales**  

*Gregor H. Mathes<sup>1</sup>, Daniele Silvestro, Kristína Kocáková, Jaime Villafaña, Catalina Pimiento*
  
<sup>1</sup>Corresponding author

# Causal Shark Extinctions Analysis

## Software Dependencies and Operating Systems
The analyses in this study were conducted in **R version 4.2.2**. The primary statistical modeling was performed using the **brms R package (version 2.19.0)** for Bayesian regression modeling. Other R packages used include:
- **dagitty (version 0.3.1)** for causal diagram analysis.
- **ggm (version 2.5)** for partial correlation testing.
- **PyRate** for estimating preservation potential and extinction risk.

**Operating system:**  
The code has been tested on **Windows 10 and 11**.

## Versions the Software Has Been Tested On
- **R version:** 4.2.2  
- **brms:** 2.19.0  
- **dagitty:** 0.3.1  
- **ggm:** 2.5  
- **PyRate:** Latest stable release (2024)  

## Required Non-Standard Hardware
- No non-standard hardware is required. The analysis can be run on a standard personal computer with **at least 8GB of RAM**.

## Installation Guide
1. Install **R (4.2.2 or higher)** from [CRAN](https://cran.r-project.org/).
2. Install **Stan** backend for Bayesian modeling if not already installed.
3. Install the necessary R packages by running the following commands:
   ```r
   install.packages("brms")
   install.packages("dagitty")
   install.packages("ggm")
   ```

Install PyRate following instructions from [PyRate's official repository](https://github.com/dsilvestro/PyRate).



## To do  

- [X] Update data readme
- [X] Calculate geographic range per bin
- [X] Repeat analysis with 1 myr bins
- [X] Use unambiguous taxonomic names
- [X] Add all models to productivity beta estimation
- [X] Repeat analysis with genus instead of species
- [X] Repeat analysis with Cenozoic subset at higher resolution
- [X] Repeat analysis with Superorders
- [X] Omit the last temporal stage 
- [X] Calculate log-odds for temperature per bin by using a mixed effect model
- [X] Clean-up and streamline according to temperature analysis
