# GS4PB (Previously SoyGen2 : Science Optimized Yield Gains Across Environments - V2) 
## An R shiny application to implement a genomic selection pipeline from quality control of genotypic data to making genomic predictions 
![PipelinePic_for_Wiki_V2](https://github.com/UMN-Lorenz-Group/SoyGen2App/assets/12753252/5e76c000-bf4e-4849-bbad-29df6a6fb22e)
 
### The recommended method to run the application is via the docker container 
#### 1) Install docker engine in your system and make sure that is running. You can download and install the docker engine for your OS from here: https://docs.docker.com/engine/install/. Once you install it, make sure the docker engine is running when you run the docker commands. An intro to running docker containers can be found in this short video: https://docs.docker.com/get-started/introduction/get-docker-desktop/.
#### 2) docker pull ivanvishnu/gs4pb:updated
#### 3) Run docker 
#### &nbsp; &nbsp; a) On gitbash: winpty docker run -d -p 3838:3838 ivanvishnu/gs4pb:updated 
#### &nbsp; &nbsp;  b) On other systems: docker run -d -p 3838:3838 ivanvishnu/gs4pb:updated 
#### 4) Access through the local link: http://localhost:3838/


### 

### II) To install and run the application in Rstudio, 
#### 1) git clone https://github.com/UMN-Lorenz-Group/GS4PB.git from git bash or terminal 
#### 2) Open App/app.R in Rstudio and run app 
#### 3) Follow the easy-to-follow instructions in the app and move through the tabs for sequential implementation of steps in the pipeline 

### III) Cyverse Implementation




### IV) Test it on gesis notebook binder by clicking the 'launch binder' icon 
#### (Since this is a public notebook binder meant for demos, it is not updated frequently and may take a long time to start)
[![Binder](https://mybinder.org/badge_logo.svg)](https://notebooks.gesis.org/binder/v2/gh/UMN-Lorenz-Group/SoyGen2App/main?urlpath=rstudio)

### V) For detailed documentation, check out the wiki page here:  
https://github.com/UMN-Lorenz-Group/GS4PB/wiki
