# DataXflowGen
DataXflowGen is a GenAI-assisted workflow for constructing, fitting and iteratively refining signed gene regulatory networks (sgGRNs) based on ordinary differential equations (ODEs). The pipeline combines large language models (LLMs) for hypothesis generation with deterministic data processing and model fitting.

## Overview

The core idea of DataXflowGen is to use GenAI where it is strong (literature‑ and knowledge‑based suggestions of genes and interactions) and to keep all quantitative validation in a mechanistic ODE framework. In practice, the workflow:

- starts from a set of genes of interest (e.g. derived from DEG analysis or ML feature selection),
- uses GenAI to propose additional candidate regulators and interactions,
- converts these suggestions into a signed interaction network,
- translates the network into an ODE model that is fitted to time‑resolved data, and
- iteratively refines the network using a back‑loop driven by model misfit.

This repository provides the scripts that implement this workflow.

## Pipeline

The pipeline consists of the following steps:

1. **GenAI‑assisted gene and interaction proposal**  
   Prompts to a GenAI assistant expand an initial list of genes of interest and generate candidate interactions. The output is a cleaned gene list and a table of gene–gene pairs.

2. **Pairwise interaction querying and encoding**  
   Each gene pair is queried individually to decide whether one gene activates, inhibits or does not regulate the other (directly or indirectly). The answers are encoded as signed integers and stored in a standardized CSV file.

3. **Network construction and pruning**  
   The interaction table is converted into a weighted adjacency matrix. Optional pruning removes redundant indirect edges if an equivalent path of direct edges exists, resulting in a more compact, interpretable network.

4. **Model fitting with Data2Dynamics (D2D)**  
   The adjacency matrix is translated into an ODE model within the Data2Dynamics framework. Parameters are fitted to experimental data, and regularization is used to simplify the model while preserving its ability to explain the data.

5. **Visualization**  
   Plots compare experimental measurement data and model simulations for each gene and summarize per‑gene chi‑square contributions to highlight where the model fits poorly.

6. **Back‑loop refinement with GenAI**  
   Genes with the highest misfit are identified, and GenAI suggests upstream regulators for these targets. The network and ODE model are updated, refit, and reevaluated until the network is stable and data-consistent.

## Repository Structure

- `Initial_model_generation_pipeline/` – GenAI-based expansion of the gene list, querying pairwise interactions, creation of an adjacency matrix, and export of the topology for the first network.
- `back_loop/` – Chi-square aggregation, regulation assignment, GenAI-guided upstream regulator search, network extension, and loop-specific adjacency updates.
- `evaluation/` – Model evaluation against random models and response consistency checks for GenAI outputs.
- `example/` – Example input data, toy settings, and a minimal working case for reproducing the pipeline.
- `extras/` – Additional helper scripts, utilities, or optional analyses not part of the core workflow.
- `filtering/` – Scripts for cleaning outputs, removing redundant or indirect interactions, and preparing filtered network files.
- `normalisation/` – Data-specific normalization of scRNA-seq time series, including log transform, scaling to [0,1], and effective standard error calculation.
- `output_preparation/` – Export of cleaned matrices, PEtab/D2D-ready files, yED topology files, and intermediate CSV outputs.
- `visualisation/` – Plotting of chi-square diagnostics and data-vs-simulation figures for publication-quality reporting.


## Getting Started

1. **Clone the repository**

```bash
git clone https://github.com/Sunbio1/DataXflowGen.git
cd DataXflowGen
```

2. **Set up the environment**

- Install Python (version as used in the manuscript).
- Install required Python packages (e.g. `openai` and standard scientific stack).  
  A `requirements.txt` or environment description can be added here.

3. **Run the pipeline**

- Start with example/ to reproduce the workflow with the provided sample data.
- Follow the scripts in Initial model generation/ to create the first signed GRN hypothesis.
- Use back_loop/ to iteratively extend and refine the network based on chi-square values.
- Inspect the outputs in output_preparation/, and visualisation/.

4. **Adapt to your own data**

- Replace the example input files in example/ with your own time-series expression data and metadata.
- Adjust the configuration, prompts, and filtering steps to match your biological context.
- Re-run the workflow to generate and validate a network tailored to your data.

## Interpretability

DataXflowGen is designed to keep the resulting ODE-based GRN interpretable: each node represents a gene, each edge represents activation or inhibition, and the fitted parameters describe mechanistic regulatory influences. GenAI is used for hypothesis generation and network extension, while model fitting remain data-driven and transparent.
