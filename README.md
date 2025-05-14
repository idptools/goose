![GOOSE_LOGO_FULL](https://github.com/idptools/goose/blob/main/images/goose_logo_3.png) 
# GOOSE : Generate disOrdered prOtiens Specifying propErties

## Last updated May 2025

### What's new (the highlights)?
* Thanks to moving to a system whereby sequences are generated using numpy vectorized operations, we no longer need to use weighted probabilities for sequence generation to keep GOOSE fast. This should dramatically improve the sequence space that GOOSE can explore when generating a protein seuqence. 
* A completely new approach for sequence / sequence variant generation has been launched! Check out the new SequenceOptimizer documentation on ReadTheDocs. 
* SPEED! GOOSE is faster. This is numpy and metapredict V3. 
* Additional functionality for sequence creation (ex. custom probabilities for amino acids when specifying sequences by properties).
* Updated minimal_variant functionality. It is  slower but does a much better job of minimizing the number of amino acids changed when trying to create your desired sequence from a starting sequence. 
* Addition of demo Jupyter Notebooks. See /goose/demos for examples on using SequenceOptimizer, generating sequences based on interactions, and making custom optimizer properties


### What is GOOSE?
GOOSE is a python package developed to make the generation of IDRs or IDR variants easy. 

## What can GOOSE do?

The main functionalities of GOOSE are:

**1.** Generate synthetic IDRs where you can specify length and...
 • Simultaneously specify *average hydrophobicity*, *fraction of charged residues (FCR)*, *net charge per residue (NCPR)*, and *kappa* (quantifies opposite charge distribution)  
 • Fractions of amino acids (multiple fractions simultaneously)  
 • End-to-end distance (Re) or radius of gyration (Rg)  
 • Interactions between an IDR and itself or other IDRs

 **2.** Generate sequences by *sequence optimization*. This is a new approach for sequence  or variant generation in GOOSE. In addition, you can **define your own functions to design sequences**!

**3.** Generate IDR variants. There are over a dozen different kinds of sequence variants in GOOSE, and they are intended to change your IDR of interest in ways that let you test various hypotheses.  

**4.** Make sequence libraries spanning sequence properties or fractions of amino acids.  

**5.** Analyze sequences. Although this is possible in GOOSE, I recommend using[SPARROW](https://github.com/idptools/sparrow) for sequence analysis. 

## How to use GOOSE

You can use GOOSE from Python or from a Google Colab notebook. The Colab notebook can be found at https://colab.research.google.com/drive/1U9B-TfoNEZbbjhPUG5lrMPS0JL0nDB3o?usp=sharing

**NOTE: May 2025** - It seems that the colab notebook is not working properly at this time. I will fix it as soon as I can. It doesn't appear to be an issue with GOOSE but rather an issue with a dependency (which makes it a little harder to debug). Apologies for any inconveniences!

## Installation - GOOSE takes flight!

Right now you can only install GOOSE through Github. We plan to put it on PyPi at some point to allow for install via pip!  

GOOSE has a few requirements **prior** to installation. Just follow the steps below to use GOOSE!  

**First install cython and numpy.**  
```bash
pip install cython
pip install numpy
```

**Now you can install GOOSE.**  

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

Documentation for GOOSE can be found [here](https://goose.readthedocs.io/en/latest/index.html)  


## How to cite GOOSE

For the time being, you can cite our [preprint](https://www.biorxiv.org/content/10.1101/2023.10.29.564547v2). It would also be helpful to link the GitHub page in your methods section so readers of your paper can understand how you generated the sequences you used.  

## Changes

The section below logs changes to GOOSE.  

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