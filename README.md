## Welcome to the GS4PB wiki!

### The GS4PB (Previously SOYGEN2) App implements a genomic selection pipeline to integrate genotyping, phenotyping, and envirotyping data. The pipeline supports GS for single traits and multi-traits across single and multiple environments.
 

![GS4PB_Intro](https://github.com/user-attachments/assets/6fb38883-0348-4fbc-a70f-e2130a6b148c)

## Recommended Method to Run the Application

### I. Using the Docker Container

1. **Install Docker Engine**  
   - Download and install the Docker engine for your operating system from [here](https://docs.docker.com/engine/install/).  
   - Ensure the Docker engine is running when executing Docker commands.  
   - For a quick introduction to running Docker containers, check out this short video: [Docker Desktop Introduction](https://docs.docker.com/get-started/introduction/get-docker-desktop/).

2. **Pull the Docker Image**  
   ```bash
   docker pull ivanvishnu/gs4pb:updated
   ```

3. **Run the Docker Container**  
   a) **On Git Bash (Windows)**  
   - **With Local Directory Binding**  
     Replace `source` with your target directory:
     ```bash
     winpty docker run --mount type=bind,source=$HOME/Downloads/,target=/root/Results -it -p 3838:3838 ivanvishnu/gs4pb:updated
     ```

   - **Without Local Directory Binding**  
     ```bash
     winpty docker run -it -p 3838:3838 ivanvishnu/gs4pb:updated
     ```
     After execution, copy the results to your local directory using:
     ```bash
     docker cp <container_id>:/root/Results/ $HOME/Downloads/
     ```

   b) **On Other Systems**  
   Replace `source` with your target directory:
   ```bash
   docker run --mount type=bind,source=$HOME/Downloads/,target=/root/Results -it -p 3838:3838 ivanvishnu/gs4pb:updated
   ```

4. **Access the Application**  
   Open your browser and navigate to [http://localhost:3838/](http://localhost:3838/).

---

### II. Installing and Running the Application in RStudio

1. Clone the repository using Git:  
   ```bash
   git clone https://github.com/UMN-Lorenz-Group/GS4PB.git
   ```

2. Open the `App/app.R` file in RStudio and run the application.

3. Follow the app’s guided instructions and move through the tabs for sequential implementation of pipeline steps.

---

### III. CyVerse Implementation
(Details to be added.)

---

### IV. Test the Application on GESIS Notebook Binder

- You can test the application by clicking the **'Launch Binder'** icon.  
  **Note:** Since this is a public demo notebook binder, it may not be frequently updated and could take a while to start.

---
## Introduction

# GS4PB (Previously SoyGen2 : Science Optimized Yield Gains Across Environments - V2) 
## An R shiny application to implement a genomic selection pipeline from quality control of genotypic data to making genomic predictions 
![PipelinePic_for_Wiki_V2](https://github.com/UMN-Lorenz-Group/SoyGen2App/assets/12753252/5e76c000-bf4e-4849-bbad-29df6a6fb22e)

### V. Detailed Documentation

For more information, refer to the wiki page:  
[GS4PB Wiki Documentation](https://github.com/UMN-Lorenz-Group/GS4PB/wiki)

