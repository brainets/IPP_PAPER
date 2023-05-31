# Information Processing Primitives

This is a reference implementation of the following paper:

Voges N., et. al. (2023) Decomposing neural circuit function into information processing primitives. JOURNAL

**preprint:** [https://www.biorxiv.org/content/10.1101/2022.08.04.502783v1.abstract](https://www.biorxiv.org/content/10.1101/2022.08.04.502783v1.abstract)

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

## Standalone C/C++ implementation

#### Copyright 2019 Johannes Hausmann, Nicole Voges (drjoe@free.fr, nicole.voges@gmx.com)

The implementation of the analysis focused on the ring model was originally provided as
a standalone C/C++ that we also make available at the **src_C** repository.

Inside this directory you will enconteur the following source codes:

- nmring.h/.cpp: 	main model, integration, and stimulus-related code.
- util.h/.cpp:       	utility macros and file system, averaging, PNG and numpy array writing methods.
- threadpool.h/.cpp:   	thread pool to run the model in a MT environment.
- oneRing_model.cpp:	main driver program for a single ring.
- fb2rings_model.cpp:	main driver program for two rings with feed-forward and -backward couplings.
- ff3rings_model.cpp	main driver program for three rings with feed-forward coupling. 

To use this version of the implementation first chang to the src_C directory:

``` 
cd src_C
``` 

and compile the code you want to execute. For compiling the code for all three models do:

```
make2
```

or a specific one (for instance the one ring configuration) by typing:

```
make oneRing_model
```

Clean compiled code by typing :

```
make clean 
```

To test functionality, comment in the lines at the beginning of each main function marked with
"// minitest of functionality" (and comment out the lines marked with "// the real thing:") and compile.
