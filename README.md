# DataXflowGen
DataXflowGen is a GenAI-assisted workflow for constructing, fitting and iteratively refining signed gene regulatory networks (sdGRNs) based on ordinary differential equations (ODEs). The pipeline combines large language models (LLMs) for hypothesis generation with deterministic data processing and model fitting.

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

5. **Visualization and diagnostics**  
   Diagnostic plots compare experimental trajectories and model simulations for each gene and summarize per‑gene chi‑square contributions to highlight where the model fits poorly.

6. **Back‑loop refinement with GenAI**  
   Genes with the highest misfit are identified, and GenAI is used again to suggest upstream regulators for these targets. The network and ODE model are extended accordingly, fitted again, and re‑evaluated. This loop can be repeated until a stable, data‑consistent network is obtained.

## Repository Structure

A possible folder structure (adapt to your repo if needed):

- `example/` – Example configuration and data (toy project or case study).
- `genai/` – Scripts for prompting GenAI, parsing responses, and generating gene lists and interaction tables.
- `network/` – Scripts for building and updating adjacency matrices, pruning indirect interactions, and exporting topologies (e.g. for yED).
- `d2d/` – Files and scripts related to Data2Dynamics (model definitions, setup, fitting, L1 regularization, intervention analysis).
- `visualization/` – Scripts for plotting chi‑square diagnostics and data vs. simulation.
- `backloop/` – Scripts for chi‑square aggregation, regulation assignment (activation/inhibition), and GenAI‑assisted network extension.
- `data/` – Example input data (expression time series, metadata) and intermediate CSV outputs.

If your structure differs, you can replace this section with a short bullet list of the actual folders you use.

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

3. **Run the example**

- Follow the instructions in `example/` (or add a short step‑by‑step description once your example is fixed).
- Inspect the intermediate CSV files (gene lists, interaction tables, adjacency matrices) and the generated plots.

4. **Adapt to your own data**

- Replace the example gene list and data files with your own.
- Adjust prompts and configuration files to your biological context and data type.
- Re‑run the pipeline to construct and refine a GRN tailored to your use case.

## Interpretability

DataXflowGen uses ODE‑based GRNs to keep the model structure and parameters biologically interpretable: each node corresponds to a gene, each edge to a regulatory influence (activation or inhibition), and each parameter to a meaningful quantity such as interaction strength. GenAI is used to generate and extend hypotheses, while the data‑driven validation and refinement are performed in a transparent, mechanistic modeling framework.
