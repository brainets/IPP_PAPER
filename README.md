# Information Processing Primitives

This is a reference implementation of the following paper:

Voges N., et. al. (2023) Decomposing neural circuit function into information processing primitives. JOURNAL

*preprint:*[https://www.biorxiv.org/content/10.1101/2022.08.04.502783v1.abstract](https://www.biorxiv.org/content/10.1101/2022.08.04.502783v1.abstract)

## Platform information and requirements

**Platform:** Ubuntu 22.04.2 LTS

**Python:** 3.9.16

**Nest:** 3.4

**Matplotlib:** 3.7.1

**Seaborn:** 0.12.2

**SciPy:** 1.10.1

**numpy:** 1.24.3

**xarray:** 2023.5.0

**frites:** 4.65.0

**jupyter-lab:** 3.3.0

**conda:**  23.3.1

## Requirements installation

The easiest way to install the requirements needed to run the codes in this repository is to install them using anaconda (https://www.anaconda.com/download) and the enviroment file provided in this repository. Otherwise, the packages can be installed manually using pip or conda. 

**the commands bellow assumes that you are using the terminal in a linux system**

First, clone this repository locally:

```
git clone https://github.com/brainets/IPP_PAPER
```

then, change to the clone directory:

```
cd IPP_PAPER
```

and install the packages via:

```
conda env create -f ipp_env.yml
```

Once the installation is finished, activate the enviroment

```
conda activate ipp
```

Now you should be able to run the codes in this repository.

### FRITES installation

The easiest way to install FRITES (remember to activate the ipp conda env before) is via pip:

```
pip install frites
```

or 

```
pip install git+https://github.com/brainets/frites
```

or 


```
git clone https://github.com/brainets/frites
cd frites
python setup.py develop
```

For more information on FRITES check  [https://github.com/brainets/frites](https://brainets.github.io/frites/).

## Repository content

This repository contains the code implementation for the ring model and the large-scale cortical model in the **src** directory:

#### Ring model

- integration.py: functions to perform numerical integration of the ring model.
- model.py: function to create the ring models and simulate them.
- utils.py: utility functions for the ring model (create simulation batches, stimulus kernels, etc).

#### Large scale model

- meanfieldcircuit.py: Implementation of the model using the Nest simulator [https://www.nest-simulator.org/](https://www.nest-simulator.org/).
- setParams.py: Function to get parameters of the model.

#### Notebooks:

The notebooks can be found in the **Notebooks** subdirectory.

- "1. One Ring Model.ipynb: Notebook to simulate the one ring model."
- "2. Two Rings Model.ipynb: Notebook to simulate the two rings (2FF) model."
- "3. Three Rings Model.ipynb: Notebook to simulate the three rings (3FF) model."
- "4. DataSet.ipynb: Load the dataset from Markov et. al., (2013)."
- "5. Mean Field circuit IPP.ipynb: Load data for the large-scale circuit and produces figure 6 from the paper."

**ATTENTION:** In order to run notebook 5, one should first simulate the large-scale model, this can be done using the **main.py** script in the main directory:

```
python main.py 2
```

but this code can be memory consuming, and you might want to get a ''bigger'' computer before doing that!

