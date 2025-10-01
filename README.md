![GOOSE_LOGO_FULL](https://github.com/idptools/goose/blob/main/images/goose_logo_3.png) 
# GOOSE : Generate disOrdered prOtiens Specifying propErties

## Last updated October 2025, latest version V0.2.5

### What's new (the highlights)?
#### October 2025 update:
* Added batch functionality in SequenceOptimizer to support batch calculation or batch prediction functionality.
* Cythonized the functionality for generating mutations in SequenceOptimizer (thanks Jeff!) 
* Improved speed of SequenceOptimizer and reduced memory usage
* Changed minimum Python version to 3.9.

#### September 2025 update:
* Complete overhaul of the SequenceOptimizer architecture to add support for optimization of properties that have highly variably seensitivies and scales of values. 
* Complete overhaul to the SequenceOptimizer properties. Added functionality to set targets to be minimum, maximum, or exact values. 
* Introduction of numerous new SequenceOptimizer properties. 
* Added linear profiles of properties in SequenceOptimizer so you can optimize across sliding windows of values across a sequence.
* Added the ability to optimize towards arbitrary vectors for linear profiles. 
* Added the ability to optimize towards arbitrary matrices for properties involving matrix calculations. 
* Update of demo notebooks in /demos to reflect changes in code. 


#### Available Demos
* **Basic optimization**: see /demos/sequence_optimization.ipynb for basic usage. 
* **Custom properties**: see /demos/custom_optimizer_peroperties.ipynb for creating and implementing custom user-defined properties  
* **Design by interaction**: see /demos/generate_sequences_by_interaction.ipynb for designing sequences to interact with a target sequence using epsilon-based properties.
* **Design by linear profiles**: see /demos/linear_profiles.ipynb for designing sequences to match linear profiles of properties like NCPR.
* **Design by interaction matrices**: see /demos/epsilon_matrix_variants.ipynb for designing sequences to match or modify interaction matrices.

### What is GOOSE?
GOOSE is a python package developed to make the generation of IDRs or IDR variants easy. 

## What can GOOSE do?

The main functionalities of GOOSE are:

**1.** Generate synthetic IDRs where you can specify length and...
 • Simultaneously specify *average hydrophobicity*, *fraction of charged residues (FCR)*, *net charge per residue (NCPR)*, and *kappa* (quantifies opposite charge distribution)  
 • Fractions of amino acids (multiple fractions simultaneously)  
 • End-to-end distance (Re) or radius of gyration (Rg)  

 **2.** Generate sequences by *sequence optimization*. There are many predefined properties you can optimzie towards. In addition, you can **define your own functions to design sequences**!

**3.** Generate IDR variants. There are now 16 different kinds of sequence variants in GOOSE, and they are intended to change your IDR of interest in ways that let you test various hypotheses.  

**4.** Analyze sequences. Although this is possible in GOOSE, I recommend using [SPARROW](https://github.com/idptools/sparrow) for sequence analysis. 

## How to use GOOSE

You can use GOOSE from Python. There is also most of the functionality in GOOSE available in a [Google Colab notebook](https://colab.research.google.com/drive/1U9B-TfoNEZbbjhPUG5lrMPS0JL0nDB3o?usp=sharing).

## Python Requirements

GOOSE requires a minimum Python version of 3.9. The primary developers all use MacBooks with Apple Silicon, which only support Python 3.9 or later. Because of this, testing Python 3.8 requires us to use a different machine. This is fine but does mean we can't find bugs specific to Python 3.8 during regular usage. Therefore, we chose to require 3.9 or later as of October 2025.

## Installation - GOOSE takes flight!

Right now you can only install GOOSE through Github. We plan to put it on PyPi at some point to allow for install via pip!  

After moving GOOSE to use pyproject.toml, you should be able to install GOOSE in a single step. 

To install directly from the git repository simply run:

```bash
pip install git+https://github.com/idptools/goose.git
```

Or to clone the GitHub repository and install locally run - 
```bash
git clone https://github.com/idptools/goose.git
cd goose
pip install .
```

**Important note**

GOOSE requires the package ``sparrow``. Sparrow should be downloaded automatically just by pip installing GOOSE, but if you have issues, try installing it by running:

```bash
pip install git+https://github.com/holehouse-lab/sparrow.git
```

This will install SPARROW. **Important note**: if your attempted install of SPARROW fails, it may be because you do not have **numpy or cython** installed. I made them both required for installation of GOOSE, so if you install GOOSE first, you should be ok. See step 1. of Installation for instructions on installing cython and numpy. If you keep having issues, please contact me or raise an issue and I'll get to it as soon as I can.


## Documentation

Documentation for GOOSE can be found [here](https://goose.readthedocs.io/en/latest/index.html).


## How to cite GOOSE

For the time being, you can cite our [preprint](https://www.biorxiv.org/content/10.1101/2023.10.29.564547v2). It would also be helpful to link the GitHub page in your methods section so readers of your paper can understand how you generated the sequences you used.  

## Changes

The section below logs changes to GOOSE.  

#### V0.2.4 and V0.2.5 - More SequenceOptimzier updates! (October 2025)
* Bug fixes in SequenceOptimizer
* Addition of teests for the SequenceOptimizer
* Implementation of Tox for automated testing across different Python versions
* Minimum Python version changed to 3.9.

#### V0.2.3 - More SequenceOptimzier updates! (October 2025)
* Added batch functionality in SequenceOptimizer to support batch calculation or batch prediction functionality.
* Cythonized the functionality for generating mutations in SequenceOptimizer (thanks Jeff!) 
* Improved speed of SequenceOptimizer and reduced memory usage

#### V0.2.2 - SequenceOptimizer update! (September 2025)
* Complete overhaul of the SequenceOptimizer architecture to add support for optimization of properties that have highly variably seensitivies and scales of values.
* Complete overhaul to the SequenceOptimizer properties. Added functionality to set targets to be minimum, maximum, or exact values. 
* Introduction of numerous new SequenceOptimizer properties. 
* Added linear profiles of properties in SequenceOptimizer so you can optimize across sliding windows of values across a sequence.
* Added the ability to optimize towards arbitrary vectors for linear profiles. 
* Added the ability to optimize towards arbitrary matrices for properties involving matrix calculations. 
* Update of demo notebooks in /demos to reflect changes in code. 

#### V0.2.1 - MAJOR UPDATE... again! (July 2025)
* Complete rewrite of variant generation functionality for both the backend code and the front end code. 
* Deprecation of library generation. If this ends up a feature people want back, I'm happy to reimplement it. 
* Addition of 11 model organism proteome IDR amino acid frequencies for use when specifying custom probabilities in sequence generation. 
* Addition of the ability to generate sequences by class of amino acids (in addition to properties and amino acid fractions).
* Overhaul of code used for generation of minimal_variant (again). 

#### V0.2.0 - MAJOR UPDATE! (May 2025)
* Overhaul of sequence generation. Sequence generation now by default does not use weighted lists. Rather, vectorized operations produce many sequences at once allowing creating of much more diverse sequences while maintaining speed of sequence generation. 
* Introduction of customizable probabilities of dictionaries when generating sequences and specifying properties. 
* Introduction of SequenceOptimizer. This allows for users to generate sequences with many combinations of parameters and create their own custom parameters. 
* Release of sequence generation by epsilon (a way to quantify sequence interaction). See https://github.com/idptools/finches/.
* Overhaul of code used for generation of minimal_variant. 
* Update to installation via pyproject.toml to enable compatbility with newer versions of Python. 
* Deprecation of generation of sequences predicted to form alpha helices, beta sheets, or beta strands

#### Version unchaged  - continued work on sequence_optimization functionality and a small bug fix (November 2024)

* Continued work on adding the ability to optimized sequences based on any user-defined parameter including ones not hardcoded into GOOSE. 
* In line documentation was completed and additional testing of the code was completed. Multiple bugs were fixed. 
* A bug where the error allowed when using the create.seq_rg() function was equal to the error for create.seq_re, which was incorrect. 
* Added functionality for sequence optimization by total epsilon and epsilon attractive / repulsive vectors

#### V0.1.4 - More bug fixes and more tests. (December 2023)

Added in tests to many thousands (maybe millions) of sequences that span sequence parameter space and sequence fraction space. Using these tests, I found some very edge case bugs in GOOSE. These bugs have been fixed. Additionally, some of the backend code has been updated to be more efficient, which will increase the speed at which GOOSE can generate sequences. Also added more information on GOOSE's ability to make sequences at the edge of kappa values. Finally, synthetic sequence creation when kappa is specified has been made to be a bit faster. 

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

Development of GOOSE was led primarily by Ryan Emenecker. Numerous ideas were contributed from others along the way including from Shahar Sukenik, Alex Holehouse, and more. Development of the SequenceOptimizer functionality released in V0.2.0 was led by Jeffrey Lotthammer. GOOSE is not possible without the amazing tools developed by the Holehouse Lab. Shoutout to Garrett Ginell, Dan Griffith, Borna Novak, and Nick Razo for their various contributions.

**Cookiecutter**
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.6.

### Copyright

Copyright (c) 2023, Ryan Emenecker - Holehouse Lab