What is GOOSE?
===============
**GOOSE : Generate disOrdered prOtiens Specifying propErties** is a python package developed to generate intrinsically disordered proteins or protein regions (collectively IDRs) and IDR variants. My goal is to make this *as useful as possible to the protein disorder community*. If you have any feature requests, please email me at remenecker@wustl.edu with the subject line **GOOSE feature request** and I'll be happy to see if I can make it happen. 

What can GOOSE do?
--------------------
The main functionalities of GOOSE are:

- Generate synthetic IDRs where you can specify length and:

  #. Simultaneously specify average hydrophobicity, fraction of charged residues (FCR), net charge per residue (NCPR), and kappa (quantifies opposite charge distribution)  
  #. Fractions of amino acids (multiple fractions simultaneously)  
  #. Classes of amino acids (multiple classes simultaneously)  
  #. End-to-end distance (Re) or radius of gyration (Rg)  

- Generate sequences by sequence optimization. This is a new approach for sequence  or variant generation in GOOSE. In addition, you can define your own functions to design sequences!
- Generate IDR variants. There are over a dozen different kinds of sequence variants in GOOSE, and they are intended to change your IDR of interest in ways that let you test various hypotheses.  

What is new in V0.2.1?
-----------------------
The higlights of the new features in GOOSE V0.2.1 are:

* Complete overhaul of the variant generation functionality to make it easier to use and less confusing.
* Major improvements to sequence generation to increase speed. 
* Introduction of custom probabilities for sequence generation. This allows you to specify the probability of each amino acid in your sequence, which can be useful for generating sequences with specific properties or characteristics.
* Addition of probabilities of amino acids in IDRs by organism based on their frequency across their proteomes. 
* Better exploration of sequence space by reducing constraints during sequence generation.
* Much faster kappa optimization and improvements to generate sequences with extreme kappa values.
* Improved error handling and more informative error messages.
* Improved documentation and examples to help users get started with GOOSE more easily.
* Many many bug fixes. 


For full details, see the change log `on GitHub <https://github.com/idptools/goose>`_ 

How can I use GOOSE?
--------------------
You can use GOOSE from Python or from Goole Colab.

We also have a Google Colab notebook. The Colab notebook can be found `in this link <https://colab.research.google.com/drive/1U9B-TfoNEZbbjhPUG5lrMPS0JL0nDB3o?usp=sharing>`_.


Installation - GOOSE takes flight!
===================================
A few notes on how to best use GOOSE.
* We strongly recommend using GOOSE in a virtual environment. This is not required, but it will help you avoid any issues with package dependencies. See conda or venv for more information on how to set up a virtual environment.
* GOOSE was tested largely in Python 3.11, so you will likely run into the fewest issues if you use that version. 
* As of GOOSE v0.2.0, GOOSE should work with Python 3.12 and above. However, if you run into issues, let us know and we will do our best to fix them.
  

The easiest way to install GOOSE is using pip.

.. code-block:: bash

    $ pip install goose

  
To install directly from the git repository simply run:

.. code-block:: bash

    $ pip install git+https://github.com/idptools/goose.git

Or to clone the GitHub repository and install locally run - 

.. code-block:: bash

    $ git clone https://github.com/idptools/goose.git
    $ cd goose
    $ pip install .

**Important note**

GOOSE requires the package ``sparrow``. Sparrow should be downloaded automatically just by pip installing GOOSE, but if you have issues, try installing it by running:

.. code-block:: bash

    $ pip install git+https://github.com/holehouse-lab/sparrow.git

This will install SPARROW. **Important note**: if your attempted install of SPARROW fails, it may be because you do not have numpy or cython installed. I made them both required for installation of GOOSE, so if you install GOOSE first, you should be ok. See step 1. of Installation for instructions on installing cython and numpy. 


Important Limitations
======================
GOOSE has some important limitations that users should be aware of. First, GOOSE makes sequences **predicted** to be disordered based on the disorder predictor metapredict. Although modern disorder predictors have proven to be *quite good*, one should aways keep in mind that predicted disorder is **not** gaurenteed disorder. 

Allowed error in sequence properties
-------------------------------------
GOOSE by default allows a *small* amount of error between some user input properties and the properties of returned sequences. For hydropathy, the allowed error is 0.07, which is honestly negligible. For kappa, allowed error is 0.03. This is a balance between accuracy and speed. You can change these values by specifying kappa_tolerance or hydropathy_tolerance in the relevant functions. In addition, if you install GOOSE locally, you can go into goose/backend/parameters and modify these values globally. Finally, if you choose an NCPR / FCR combination that is mathematically impossible, GOOSE will get as close as it can.

Speed, specified properties, and stochasticity
-----------------------------------------------
The protein disorder field moves fast, and we are not here to slow your research down. It was important for us to make GOOSE as fast as possible. However, because GOOSE incorporates stochasticity into sequence generation, GOOSE still has to do some work when designing your disordered sequence. The stochasticity in sequence generation makes it harder for GOOSE to generate sequences but helps minimize the chance that GOOSE makes the same sequence more than once. This is important because it allows you to create many sequences or sequence variants with the exact same overall properties but different primary sequences. As far as speed goes, *the more properties you specify, or the more constraints you put on sequence design, the more time it will take GOOSE to generate your sequence*. 

Failed sequence generation
---------------------------
Sometimes GOOSE can't make your sequence. However, you can usually just run the code a few more times and GOOSE will eventually land on a solution that matches your specified properties (thanks to the inherent stochasticity in sequence generation). The reason we designed GOOSE this way is to avoid situations where you try to make a sequence that is difficult for GOOSE to generate and GOOSE spends 10+ minutes working it out. If you still can't get a sequence you want, try *slightly* adjusting your properties or reducing the disorder cutoff value. 

Limits on specifying sequence properties
-----------------------------------------
GOOSE will only return sequences with disorder values above the cutoff disorder threshold. Some sequence compositions (for example, very high mean hydrophobicity) are simply not predicted to be disordered. GOOSE will not by default return these sequences to you. Apart from sequences not predicted to be disordered, it is also important to note that some combinations of sequence properties are not mathematically possible. GOOSE uses a rescaled Kyte Doolittle hydropathy scale for calculating mean hydrophobicity. This scale goes from 0 to 9 where higher values are more hydrophobic. The charged residues have low hydrophobicity values (R = 0, K = 0.6, D = 1, E = 1). Therefore, if you have a sequence with too many charged residues, you limit how high the mean hydrophobicity can go. If you specify a high FCR and a high hydrophobicity, that sequence may be mathematically impossible to make. GOOSE will return an error if you do this. 

Best practices when using GOOSE
--------------------------------
It is best practice to double check that the sequences you make using GOOSE are what you intended. You can do this using the *analyze* module included with GOOSE. Although we have done extensive testing on GOOSE functionality, due to the massive sequence space that is possible when generating an IDR, you may encounter bugs. We would appreciate if you would report these bugs, and we will do our best to fix them as quickly as possible.
