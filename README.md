# General Insurance Risk Analysis Project

## Overview

This project involved the analysis of insurance risk models, focusing on loss severity data and the financial stability of an insurance company. The project was divided into two main parts: **Part 1** focused on fitting probability distributions to loss severity data, while **Part 2** analyzed the insurer's surplus process and the impact of reinsurance strategies on ruin probability.

## Part 1: Loss Severity Analysis

### Objective
The goal was to analyze a dataset of 2,000 property insurance claims to determine the most appropriate probability distribution for modeling loss severity. The analysis included both complete and censored data scenarios.

### Methodology
1. **Model Fitting:** Maximum Likelihood Estimation (MLE) was used to fit four probability distributions (**Log-normal**, **Exponential**, **Gamma**, and **Pareto**) to the loss data.
2. **Goodness-of-Fit Tests:** Graphical approaches and hypothesis tests (K-S and A-D tests) were employed to evaluate the quality of the fitted models.
3. **Censored Data Analysis:** MLE estimates were calculated for the **Exponential** and **Pareto** distributions under a policy limit of $10,000. QQ plots were produced to visually assess the fit.

### Key Findings
- The best-fitting distribution for the loss severity data was identified based on statistical criteria and graphical analysis.
- For censored data, the **Exponential** and **Pareto** distributions were evaluated, and their likelihood functions were derived.

## Part 2: Insurer's Surplus Process and Reinsurance Analysis

### Objective
The task focused on assessing the financial stability of an insurance company by analyzing its surplus process and ruin probability. The impact of two reinsurance strategies—**Proportional Reinsurance** and **Excess of Loss (EoL) Reinsurance**—was evaluated.

### Methodology
1. **Ruin Probability Estimation:** The ruin probability over 5 years was estimated under two scenarios:
   - **Scenario I:** Exponential distribution for claim sizes.
   - **Scenario II:** Gamma distribution for claim sizes.
2. **Reinsurance Analysis:** The impact of reinsurance on the ruin probability was analyzed by:
   - Calculating the ruin probability with reinsurance.
   - Determining the range of retention levels (α for Proportional Reinsurance and d for EoL Reinsurance) to ensure financial stability.
   - Maximizing the adjustment coefficient to minimize the ruin probability.
3. **Sensitivity Analysis:** The sensitivity of the ruin probability to key parameters (initial surplus, premium loading factors) was analyzed to recommend adjustments to the reinsurance strategy.

### Key Findings
- The ruin probability was estimated under both scenarios, and the effectiveness of reinsurance in reducing ruin probability was demonstrated.
- The optimal retention levels for Proportional and EoL Reinsurance were determined to maximize financial stability.
- Sensitivity analysis provided insights into how changes in key parameters impact the ruin probability, guiding recommendations for premium rates and reinsurance strategies.

## Conclusion

This project successfully applied statistical and actuarial techniques to analyze loss severity data and assess the financial stability of an insurance company. The findings provided valuable insights into the selection of appropriate probability distributions for loss modeling and the impact of reinsurance strategies on ruin probability. The results were supported by rigorous statistical analysis and sensitivity testing, offering actionable recommendations for improving the insurer's financial stability.

