## Welcome to the GS4PB wiki!

### The GS4PB (Previously SOYGEN2) App implements a genomic selection pipeline to integrate genotyping, phenotyping, and envirotyping data. The pipeline supports GS for single traits and multi-traits across single and multiple environments [(GS4PB Manuscript)](https://acsess.onlinelibrary.wiley.com/doi/10.1002/tpg2.70150). 
 

![GS4PB_Intro](https://github.com/user-attachments/assets/6fb38883-0348-4fbc-a70f-e2130a6b148c)

## Table of Contents

| # | Version | Description |
|---|---------|-------------|
| [I](#i-using-the-docker-container) | **Docker Container** | Run the Shiny app via Docker — recommended for most users |
| [II](#ii-installing-and-using-the-gs4pb-r-package) | **GS4PB R Package** | Install the package in R; run the Shiny app or use functions programmatically in scripts / HPC pipelines |
| [III](#iii-running-the-app-on-an-hpc-cluster-singularity--apptainer) | **HPC (Singularity / Apptainer)** | Run the containerised Shiny app on an HPC cluster |
| [IV](#iv-cyverse-implementation) | **CyVerse** | CyVerse cloud implementation *(details coming soon)* |
| — | [**Trial Data**](#trial-data) | Download example SE and ME datasets |
| [V](#v-detailed-documentation) | **Detailed Documentation** | Full wiki and function reference |

---

## Recommended Method to Run the Application

### I. Using the Docker Container

1. **Install Docker Engine**  
   - Download and install the Docker engine for your operating system from [here](https://docs.docker.com/engine/install/).  
   - Ensure the Docker engine is running when executing Docker commands.  
   - For a quick introduction to running Docker containers, check out this short video: [Docker Desktop Introduction](https://docs.docker.com/get-started/introduction/get-docker-desktop/).

2. **Pull the Docker Image**  
   ```bash
   docker pull umnlorenzgroup/gs4pb:updated
   ```

3. **Run the Docker Container**  
   a) **On Git Bash (Windows)**  
   - **With Local Directory Binding**  
     Replace `source` with your target directory:
     ```bash
     winpty docker run --mount type=bind,source=$HOME/Downloads/,target=/root/Results -it -p 3838:3838 umnlorenzgroup/gs4pb:updated
     ```

   - **Without Local Directory Binding**  
     ```bash
     winpty docker run -it -p 3838:3838 umnlorenzgroup/gs4pb:updated
     ```
     After execution, copy the results to your local directory using:
     ```bash
     docker cp <container_id>:/root/Results/ $HOME/Downloads/
     ```

   b) **On Other Systems**  
   Replace `source` with your target directory:
   ```bash
   docker run --mount type=bind,source=$HOME/Downloads/,target=/root/Results -it -p 3838:3838 umnlorenzgroup/gs4pb:updated
   ```

4. **Access the Application**  
   Open your browser and navigate to [http://localhost:3838/](http://localhost:3838/).

---

### II. Installing and Using the GS4PB R Package

GS4PB can be installed as a standard R package and used either through its Shiny interface or programmatically in R scripts and HPC pipelines.

#### Installation

1. Clone the repository using Git:
   ```bash
   git clone https://github.com/UMN-Lorenz-Group/GS4PB.git
   ```

2. Open RStudio and navigate to the root of the package folder. Ensure `conda` is available (required for AlphaPlantImpute2 imputation):
   ```R
   Sys.which("conda")
   # If empty like "":
   reticulate::install_miniconda()
   ```

3. Build and install the package using devtools:
   ```R
   devtools::install_deps(dependencies = TRUE)
   devtools::document()
   # ℹ Updating GS4PB documentation
   # ℹ Loading GS4PB
   # Warning message:
   #   In nsenv[[f_name]](dirname(ns_path), package) :
   #   Conda environment not found. Run GS4PB:::setup_python_env() to create it.
   ### Ignore this warning for now and build the package
   devtools::build()
   devtools::install()
   ```

4. Set up the Python environment (required only for AlphaPlantImpute2 imputation):
   ```R
   library(GS4PB)
   GS4PB:::setup_python_env()   # creates the "GS4PB_CondaEnv" conda environment
   ```

#### Running the Shiny App

Once installed, launch the interactive Shiny interface:
```R
options(java.parameters = "-Xmx30g")   # set BEFORE loading any Java-dependent package
library(GS4PB)
GS4PB::run_app()
```
Follow the app's guided tabs for sequential pipeline steps.

#### Programmatic / Script Use

GS4PB functions can also be called directly in R scripts — useful for HPC cluster jobs or automated pipelines. Set the Java heap size and initialise rTASSEL before any genotype-loading call:
```R
options(java.parameters = "-Xmx30g")   # adjust to available RAM
library(GS4PB)
rTASSEL::startLogger(fullPath = NULL, fileName = NULL)
rTASSEL::initializeTASSEL()
```

**Workflow Vignettes**

| Workflow | Vignette | Preview | Description |
|----------|----------|---------|-------------|
| **SE: Single-Environment** | `vignette("GS4PB_Pipeline", package = "GS4PB")` | [View on GitHub](https://github.com/UMN-Lorenz-Group/GS4PB/blob/main/vignettes/GS4PB_Pipeline.Rmd) | End-to-end pipeline for single-environment single-trait (SE:ST) and multi-trait (SE:MT) predictions: VCF loading, QC filtering, imputation, optional training-set optimisation, cross-validation, and GEBV ranking. |
| **ME: Multi-Environment** | `vignette("GS4PB_ME_Workflow", package = "GS4PB")` | [View on GitHub](https://github.com/UMN-Lorenz-Group/GS4PB/blob/main/vignettes/GS4PB_ME_Workflow.Rmd) | Advanced multi-environment workflow for HPC: ME cross-validation (CV0/CV1/CV2/CV00), leave-one-test-out CV (LOTCV), LOFO prediction with partial masking, multi-trait extension. Includes parallel `foreach` patterns and SLURM setup. |

A copy of the raw executable pipeline script (no prose) is also installed with the package:
```R
file.edit(system.file("pipeline", "GS4PB_Pipeline.R", package = "GS4PB"))
```

---

### Trial Data 
1. Trial data for single environmental trial can be found here:

   https://data.cyverse.org/dav-anon/iplant/home/umnlorenzgroup/SingleEnv_Trial_Data.zip
   
   This folder contains three files: genotypic data file in VCF format, a '.csv' file
   containing target line IDs and a phenotypic data file in '.csv' format. 

3. Trial data for multi-environmental trials can be found here:

   https://data.cyverse.org/dav-anon/iplant/home/umnlorenzgroup/MultiEnv_Trial_Data.zip 

   The multi-environmental trial data folder contains three files: a genotypic data file in VCF format, a phenotypic data file in '.csv' format and
   location coordinates file for retrieving enviromics (weather) data 

### III. Running the App on an HPC Cluster (Singularity / Apptainer)

For users with access to an HPC cluster, GS4PB can be run as a containerised
Shiny app using [Singularity](https://sylabs.io/singularity/) or
[Apptainer](https://apptainer.org/) — no Docker or local R installation required.

Full setup instructions and launcher scripts are maintained in a dedicated
helper repository:

**[GS4PBAppHPCHelp](https://github.com/UMN-Lorenz-Group/GS4PBAppHPCHelp)**

#### Quick Start

1. Clone the helper repository and place the `gs4pb_updated.sif` container
   image in the same directory:
   ```bash
   git clone https://github.com/UMN-Lorenz-Group/GS4PBAppHPCHelp.git
   cd GS4PBAppHPCHelp
   # copy or symlink gs4pb_updated.sif here
   ```

2. Make the launcher script executable and run it:
   ```bash
   chmod +x run_app.sh
   ./run_app.sh
   ```

3. Access the application at [http://localhost:3838/](http://localhost:3838/).

Configuration options (port, bind address, results directory) are set in
`config.env`. The script automatically creates an isolated `runtime/`
directory for container home, tmp, and results.

> **Requirement:** Apptainer or Singularity must be available as a module or
> system install on your cluster. Check with `which apptainer` or
> `which singularity`.

---

### IV. CyVerse Implementation
(Details to be added.)

---

### V. Detailed Documentation

For more information, refer to the wiki page:  
[GS4PB Wiki Documentation](https://github.com/UMN-Lorenz-Group/GS4PB/wiki)

---
## Introduction

# GS4PB (Previously SoyGen2 : Science Optimized Yield Gains Across Environments - V2) 
## An R shiny application to implement a genomic selection pipeline from quality control of genotypic data to making genomic predictions 
![PipelinePic_for_Wiki_V2](https://github.com/UMN-Lorenz-Group/SoyGen2App/assets/12753252/5e76c000-bf4e-4849-bbad-29df6a6fb22e)
