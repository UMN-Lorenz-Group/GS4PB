# SoyGen2 (Science Optimized Yield Gains Across Environments - V2) 
## An R shiny application to implement a genomic selection pipeline from quality control of genotypic data to making genomic predictions 
![PipelinePic_for_Wiki_V2](https://github.com/UMN-Lorenz-Group/SoyGen2App/assets/12753252/5e76c000-bf4e-4849-bbad-29df6a6fb22e)
 
### The recommended method to run the application is via the docker container 
#### 1) Install docker engine in your system and make sure that is running 
#### 2) docker pull ivanvishnu/soygen2:updated
#### 3) 
#### a) On gitbash: winpty docker run -d -p 3838:3838 ivanvishnu/soygen2:updated 
####   b) On other systems: docker run -d -p 3838:3838 ivanvishnu/soygen2:updated 
#### 4) Access through the local link: http://localhost:3838/


### 

### II) To install and run the application in Rstudio, 
#### 1) git clone https://github.com/UMN-Lorenz-Group/SoyGen2App.git from git bash or terminal 
#### 2) Open App/app.R in Rstudio and run app 
#### 3) Follow the easy-to-follow instructions in the app and move through the tabs for sequential implementation of steps in the pipeline




### III) Test it on gesis notebook binder by clicking the 'launch binder' icon 
#### (Since this is a public notebook binder meant for demos, it is not updated frequently and may take a long time to start)
[![Binder](https://mybinder.org/badge_logo.svg)](https://notebooks.gesis.org/binder/v2/gh/UMN-Lorenz-Group/SoyGen2App/main?urlpath=rstudio)
