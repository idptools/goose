![GOOSE_LOGO_FULL](https://github.com/ryanemenecker/goose/blob/main/images/goose_logo_3.png) 
# GOOSE : Generate disOrdered prOtiens Specifying propErties

### Last updated October 2023

### What is GOOSE?
GOOSE is a python package developed to make the generation of IDRs or IDR variants easy. Basically, we want to make it easy for you to design IDRs or IDR variants for your research.

## What can GOOSE do?

There are four main functionalities currently in GOOSE. 

**1.** Generatae IDRs where you can specify the properties *length*, *average hydropathy*, *fraction of charged residues (FCR)*, *net charge per residue (NCPR)*, and *kappa*, which is a property that defines opposite charge distribution in a sequence.

**2.** Generate IDR variants.

**3.** Make sequence libraries. This includes libraries of sequences spanning different sequence properties or fractions of amino acids. 

**4.** Analyze sequences. 

## How can I use GOOSE?

You can use GOOSE from Python or from a Google Colab notebook. The Colab notebook can be found at https://colab.research.google.com/drive/1U9B-TfoNEZbbjhPUG5lrMPS0JL0nDB3o?usp=sharing

## Installation - GOOSE takes flight!

Right now you can only install GOOSE through Github. It will be on PyPi to allow for pip installation soon!

To install directly from the git repository simply run:

	$ pip install git+https://github.com/idptools/goose.git

Or to clone the GitHub repository and install locally run - 

	$ git clone https://github.com/idptools/goose.git
	$ cd goose
	$ pip install .

**Important note**

GOOSE also requires the package ``sparrow``. Sparrow should be downloaded automatically just by pip installing GOOSE, but if you have issues, try installing it by running:

	$ pip install git+https://github.com/holehouse-lab/sparrow.git

This will install SPARROW. **Important note**: if your attempted install of SPARROW fails, it may be because you do not have numpy or cython installed. I made them both required for installation of GOOSE, so if you install GOOSE first, you should be ok. If you still have problems, try running from terminal:

	$ pip install cython
	$ pip install numpy


## Documentation

Documentation for GOOSE can be found at https://goose.readthedocs.io/en/latest/index.html.


## How to cite GOOSE

You can not currently cite GOOSE as we have yet to publish it (hopefully soon!). We would appreciate if you would mention GOOSE in your methods section with a link to the Github page so readers of your paper can understand how you generated the sequences you used.

## Changes

The section below logs changes to GOOSE. 

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