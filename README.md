---
title: "README.md"
author: "Stephen Tivenan"
date: "2026-05-08"
output: html_document
---


Global Envelope Test: Interactive Simulation

This repository contains a Shiny dashboard designed to visualize the **Maximum Absolute Deviant Global Envelope Test** using B-splines and Gaussian Process (GP) noise.

How to Run
To run this application locally, ensure you have R installed and run:

```R
install.packages(c("shiny", "MASS", "splines2", "plotly", "dplyr"))
shiny::runGitHub("YOUR_USERNAME/Global-Envelope-Simulation")
