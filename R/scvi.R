# # do we have reticuate library?
# library(reticulate)
# #what is our python env?
# py_config()
# 
# #what is our conda envs?
# conda_list()
# 
# 
# ## in the terminal 
# 
# # brew install miniconda
# 
# ## Then restart your terminal and initialize Conda:
# 
# # conda init zsh  # Or "conda init bash" depending on your shell
# 
# ## Then restart your terminal and you are in the new Conda env
# 
# # conda --version
# 
# 
# # conda create -n scvi-env python=3.10
# 
# # conda activate scvi-env
# # pip install scvi-tools
# 
# 
# ## back to R studio 
# file.edit("~/.Renviron")
# 
# ## and enter: 
# # RETICULATE_PYTHON=/opt/homebrew/Caskroom/miniconda/base/envs/scvi-env/bin/python
# ## save and restart R
# 
# 
# 
# ### now you should be able to 
# 
# 
# Sys.setenv(RETICULATE_PYTHON = "/opt/homebrew/Caskroom/miniconda/base/envs/scvi-env/bin/python")
# 
# 
# library(reticulate)
# 
# Sys.getenv("RETICULATE_PYTHON")  # Should return Miniconda Python path
# py_config()
# 
# # library(reticulate)
# use_condaenv("scvi-env", required = TRUE)
# use_python("/opt/homebrew/Caskroom/miniconda/base/envs/scvi-env/bin/python", required = TRUE)
# py_config()
# py_module_available("scvi")
# scvi <- import("scvi", convert = FALSE)  # Try loading manually
# 
# 
# system("which python")
# system("python -c 'import scvi; print(scvi.__version__)'")
# 
# 
# 
