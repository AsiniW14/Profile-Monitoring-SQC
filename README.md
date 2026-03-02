Profile Monitoring in Statistical Quality Control

This repository contains an R implementation of Profile Monitoring, a specialized branch of Statistical Process Control (SPC). While traditional SPC monitors single quality characteristics, Profile Monitoring tracks the stability of the functional relationship between a dependent variable ($y$) and one or more independent variables ($x$) over time.

Project Overview

In many modern manufacturing processes, quality is defined by a profile (a curve or function) rather than a single point. This project focuses on Simple Linear Profiles, where the relationship is modeled as:

$$y_{ij} = B_0 + B_1 x_{ij} + \epsilon_{ij}$$

The code simulates a two-phase monitoring approach to ensure process stability and detect shifts in these parameters.

Features

Phase I Analysis: Establishes a stable "In-Control" (IC) baseline by estimating parameters from 30 initial profiles.
Phase II Monitoring: Performs real-time monitoring of new incoming profiles to detect process disturbances, specifically shifts in the slope ($\beta_1$).
Hotelling’s $T^2$ Control Chart: Utilizes a multivariate control chart to monitor the intercept and slope simultaneously.
Visualization: Generates clear, professional plots using ggplot2 to visualize Phase I stability and Phase II out-of-control signals.

Code Structure

The R script SQC_147972_RCode.R is organized as follows:
Global Parameters: Sets the seed for reproducibility and defines the true IC parameters ($\beta_0 = 5, \beta_1 = 2$).
Phase I Generation: Simulates 30 profiles, fits linear models to each, and calculates the Upper Control Limit (UCL) based on the $F$-distribution.
Phase II Simulation: Introduces a shift in the slope (from 2.0 to 2.2) starting at profile 11 to test the sensitivity of the monitoring system.
Results & Visualization: Plots the $T^2$ statistics against the UCL to identify exactly when the process goes "Out-of-Control".

Getting Started

Prerequisites
You will need R installed on your system, along with the ggplot2 library. - install.packages("ggplot2")

Running the Analysis
Simply source the script in your R environment: source("SQC_147972_RCode.R")

Results
The simulation successfully detects the shift in the process. In the Phase II chart, you will observe the $T^2$ statistic exceeding the UCL shortly after the shift is introduced at Profile 11, demonstrating the effectiveness of the Hotelling $T^2$ approach for linear profiles.

Developed as part of the Statistical Quality Control Seminar at Europa Universität Viadrina.
