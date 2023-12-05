![GOOSE_LOGO_FULL](https://github.com/ryanemenecker/goose/blob/main/images/goose_logo_3.png) 
# GOOSE : Generate disOrdered prOtiens Specifying propErties

### Last updated December 2023

### What is GOOSE?
GOOSE is a python package developed to make the generation of IDRs or IDR variants easy. Basically, we want to make it easy for you to design IDRs or IDR variants for your research.

## What can GOOSE do?

There are four main functionalities currently in GOOSE. 

**1.** Generate synthetic IDRs where you can specify the length and:  
 • properties including *average hydropathy*, *fraction of charged residues (FCR)*, *net charge per residue (NCPR)*, and *kappa* (which is a property that defines opposite charge distribution in a sequence),  
 • fractions of amino acids (you can specify multiple fractions simultaneously), or  
 • ensemble dimensions, either end-to-end distance (Re) or radius of gyration (Rg)  

**2.** Generate IDR variants. There are over 16 different kinds of sequence variants in GOOSE, and they are intended to change your IDR of interest in ways that let you test different hypotheses.  

**3.** Make sequence libraries spanning sequence properties or fractions of amino acids.  

**4.** Analyze sequences.  

## How can I use GOOSE?

You can use GOOSE from Python or from a Google Colab notebook. The Colab notebook can be found at https://colab.research.google.com/drive/1U9B-TfoNEZbbjhPUG5lrMPS0JL0nDB3o?usp=sharing

## Installation - GOOSE takes flight!

Right now you can only install GOOSE through Github. It will be on PyPi to allow for pip installation soon!  

GOOSE has a few requirements **prior** to installation. Just follow the steps below to use GOOSE!  

1. Install cython and numpy.  

	$ pip install cython
	$ pip install numpy

2. Install GOOSE.  

To install directly from the git repository simply run:

	$ pip install git+https://github.com/idptools/goose.git

Or to clone the GitHub repository and install locally run - 

	$ git clone https://github.com/idptools/goose.git
	$ cd goose
	$ pip install .

**Important note**

GOOSE requires the package ``sparrow``. Sparrow should be downloaded automatically just by pip installing GOOSE, but if you have issues, try installing it by running:

	$ pip install git+https://github.com/holehouse-lab/sparrow.git

This will install SPARROW. **Important note**: if your attempted install of SPARROW fails, it may be because you do not have numpy or cython installed. I made them both required for installation of GOOSE, so if you install GOOSE first, you should be ok. See step 1. of Installation for instructions on installing cython and numpy. 


## Documentation

Documentation for GOOSE can be found at https://goose.readthedocs.io/en/latest/index.html.  


## How to cite GOOSE

You can not currently cite GOOSE as we have yet to publish it (hopefully soon!). We would appreciate if you would mention GOOSE in your methods section with a link to the Github page so readers of your paper can understand how you generated the sequences you used.  

## Changes

The section below logs changes to GOOSE.  

#### V0.1.3 - Bug fixes and more tests. (December 2023)

Wrote up tests for nearly all possible combinations of creating sequences by specifying different numbers of sequence properties. Found some bugs along the way and fixed them (as far as I can tell). Generally made sequence generation more robust across parameter space. Fixed some issues with exceptions. Updated minimum kappa parameter to a value of 0.03 because GOOSE often has a hard time getting to 0.0. Also increased max hydropathy value that GOOSE can generate to 6.6 on the adjusted Kyte-Doolittle scale. 

#### V0.1.2 - New features, bug fixes, general improvement. (October 2023)

Made sequence generation by specifying radius of gyration or end-to-end distance user-facing. It was hiding in the backend code for some time now but I finally had time to make the user-facing functionality and write up the docs :). The user facing functions are called ``seq_re()`` and ``seq_rg()``. Added sequence variants that alter Re and Rg called ``re_var()`` and ``rg_var()``, respectively. Added helicity predictions to ``analyze`` functionality.

#### V0.1.1 - Bug fix and new feature (October 2023)

Fixed bug where sequences with high FCR specified could not have a high enough hydropathy value using the ``create.sequence()`` function. Added ability to exclude residues from sequence generated using the ``create.sequence()`` function (cannot exclude charged residues when FCR / NCPR are specified). Improved range of kappa values that could be specified. Improved performance of some backend functionality. Fixed some additional edge case bugs. Changed name of ``shuffle_var()`` to ``region_shuffle_var()`` due to previous addition of other types of shuffle variants. Added more tests.

#### V0.1.0 - Initial release (September 2023)

Initial release. Begin tracking changes.

## Acknowledgements

**Funding**
This project was supported by grants from agencies including the National Science Foundation (NSF), the National Institute of Health (NIH).

**Development**
Development of GOOSE was led primarily by Ryan Emenecker. Numerous ideas were contributed from others along the way including from Shahar Sukenik, Alex Holehouse, and many more. 

**Cookiecutter**
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.6.

### Copyright

Copyright (c) 2023, Ryan Emenecker - Holehouse Lab