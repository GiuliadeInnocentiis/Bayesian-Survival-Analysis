# Bayesian-Survival-Analysis
Bayesian models were implemented for the survival analysis of patients with breast cancer. To evaluate the effect of covariates on survival time, Exponential, Weibull, and Log-Normal models with fixed effects, as well as Weibull and Log-Normal models with random effects, were estimated.

## Overview

Survival analysis represents a set of fundamental statistical techniques for studying the time that elapses between a starting point such as the beginning of a clinical study—and the occurrence of an event of interest, such as a patient's death or disease recurrence. This type of analysis is widely used in biomedical, epidemiological, and biological research.

The objective of survival analysis is to describe the distribution of time-to-event in the population under study and to quantify the impact of explanatory factors on survival.

Traditionally, this field has relied on frequentist approaches, such as the Cox model or classical parametric distributions, which provide powerful tools but have limitations in handling structural complexity in the data. In recent decades, the Bayesian approach has gained increasing attention. By explicitly using prior distributions and updating prior knowledge through the posterior distribution, Bayesian modeling offers a more transparent treatment of uncertainty and greater flexibility in model specification.

In this analysis, three Bayesian parametric models and one semi-parametric model are presented. The differences among them are discussed, the effects of covariates on time-to-event are interpreted, and their predictive performance is compared.
