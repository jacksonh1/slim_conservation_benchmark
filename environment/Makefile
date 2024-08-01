#################################################################################
# GLOBALS                                                                       #
#################################################################################
.ONESHELL:

CONDA = mamba
SHELL = /bin/bash
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate
PROJECT_NAME = slim_conservation
# PYTHON_VERSION = 3.10

#################################################################################
# COMMANDS                                                                      #
#################################################################################

## create environment
.PHONY: build_environment
build_environment:
	$(CONDA) create --name $(PROJECT_NAME)
	$(CONDA_ACTIVATE) $(PROJECT_NAME)
	$(CONDA) env update --name $(PROJECT_NAME) --file ./environment.yml
	@echo ">>> conda env created. Activate with:\nconda activate $(PROJECT_NAME)"

## create environment for mac x86
.PHONY: build_environment_mac_x86
build_environment_mac_x86:
	CONDA_SUBDIR=osx-64 conda create -n $(PROJECT_NAME)
	$(CONDA_ACTIVATE) $(PROJECT_NAME)
	conda config --env --set subdir osx-64
	$(CONDA) env update --name $(PROJECT_NAME) --file ./environment.yml
	@echo ">>> conda env created. Activate with:\nconda activate $(PROJECT_NAME)"

## check architecture and run appropriate target
.PHONY: environment
environment:
	if [ `uname -m` = "arm64" ]; then make build_environment_mac_x86; else make build_environment; fi




# ## create environment
# .PHONY: build_environment
# build_environment:
# 	$(CONDA) create --name $(PROJECT_NAME)
# 	$(CONDA_ACTIVATE) $(PROJECT_NAME)
# 	$(CONDA) install python=$(PYTHON_VERSION)
# 	$(CONDA) env update --name $(PROJECT_NAME) --file ./slim_conservation_orthogroup_generation/environment.yml
# 	$(CONDA) env update --name $(PROJECT_NAME) --file ./slim_conservation_scoring/environment.yml
# 	@echo ">>> conda env created. Activate with:\nconda activate $(PROJECT_NAME)"

# ## create environment for mac x86
# .PHONY: build_environment_mac_x86
# build_environment_mac_x86:
# 	CONDA_SUBDIR=osx-64 conda create -n $(PROJECT_NAME)
# 	$(CONDA_ACTIVATE) $(PROJECT_NAME)
# 	conda config --env --set subdir osx-64
# 	$(CONDA) install python=$(PYTHON_VERSION)
# 	$(CONDA) env update --name $(PROJECT_NAME) --file ./slim_conservation_orthogroup_generation/environment.yml
# 	$(CONDA) env update --name $(PROJECT_NAME) --file ./slim_conservation_scoring/environment.yml
# 	@echo ">>> conda env created. Activate with:\nconda activate $(PROJECT_NAME)"





