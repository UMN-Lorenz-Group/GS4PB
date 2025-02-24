
if(!require("renv",quiely=TRUE)){
  install.packages("renv")
 }

library(renv)

#### Create a project library & link all the required libraries to the project library
renv::init()

#### Set python version to be used 
renv::use_python()
## Selected 5. C:/Users/ivanv/.pyenv/pyenv-win/shims/python3.bat

reticulate::use_virtualenv("C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2App/App/renv/python/virtualenvs/renv-python-3.12")
virtualenv_install("C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2App/App/renv/python/virtualenvs/renv-python-3.12","numba") 
virtualenv_install("C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2App/App/renv/python/virtualenvs/renv-python-3.12","C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/AlphaPlantImpute2/alphaplantimpute2-1.5.3-py3-none-any.whl") 


#### Inconsistency with one package

remove.packages("Matrix")
renv::install("Matrix")
renv::snapshot(packages = "Matrix")



renv::dependencies() 


# Create a .renvignore file in the project directory
echo "Data/" >> .renvignore
echo "Output/" >> .renvignore


renv::rehash()
renv::restore()

