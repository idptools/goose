What is GOOSE?
===============
**GOOSE : Generate disOrdered prOtiens Specifying propErties** is a python package developed to generate intrinsically disordered proteins or protein regions (collectively IDRs) and IDR variants.

What can GOOSE do?
--------------------
The main functionalities of GOOSE are:

- Generate synthetic IDRs where you can specify length and:

  #. Simultaneously specify average hydrophobicity, fraction of charged residues (FCR), net charge per residue (NCPR), and kappa (quantifies opposite charge distribution)  
  #. Fractions of amino acids (multiple fractions simultaneously)  
  #. Classes of amino acids (multiple classes simultaneously)  
  #. End-to-end distance (Re) or radius of gyration (Rg)  

- Generate IDR variants. There are over a dozen different kinds of sequence variants in GOOSE, and they are intended to change your IDR of interest in ways that let you test various hypotheses.  
- Generate sequences by sequence optimization.  In addition, you can define your own functions to design sequences!

  #. This approach offers greater flexibility and control over the design process.
  #. The optimizer comes with many built-in functions for optimization.  
  #. You can define your own functions to design sequences.

Importantly, GOOSE incorporates stochasticity into sequence generation, which allows you to create many sequences or sequence variants with the exact same overall properties but different primary sequences. 

What is new in V0.2.2?
-----------------------
The highlights of the new features in GOOSE V0.2.2 are:

- Complete overhaul of the SequenceOptimizer architecture:

  #. Better support for optimization of properties that have highly variably sensitivities and/or scales of values.
  #. Better logging and debugging information.
  #. Adaptive optimization that can help avoid getting 'stuck' in local minima.

- Rewrite of the SequenceOptimizer properties.

  #. Added functionality to set targets values to be minimum, maximum, or exact values. 
  #. Added linear profiles of properties in SequenceOptimizer so you can optimize across sliding windows of values across a sequence.
  #. Added the ability to optimize towards arbitrary vectors for linear profiles. 
  #. Added the ability to optimize towards arbitrary matrices for properties involving matrix calculations. 
  #. Introduction of numerous new SequenceOptimizer properties. 

- Update of demo notebooks in /demos to reflect changes in code and provide new demos. 
- Numerous bug fixes and performance improvements.


For full details, see the change log `on GitHub <https://github.com/idptools/goose>`_ 

How can I use GOOSE?
--------------------
You can use GOOSE from Python or from Google Colab. The Colab notebook can be found `in this link <https://colab.research.google.com/drive/1U9B-TfoNEZbbjhPUG5lrMPS0JL0nDB3o?usp=sharing>`_.


Installation - GOOSE takes flight!
===================================
Right now you can only install GOOSE through Github. 

A few notes on how to best use GOOSE.
* We strongly recommend using GOOSE in a virtual environment. This is not required, but it will help you avoid any issues with package dependencies. See conda or venv for more information on how to set up a virtual environment.
* GOOSE was tested largely in Python 3.11, so you will likely run into the fewest issues if you use that version. 
* As of GOOSE v0.2.0, GOOSE should work with Python 3.12 and above. However, if you run into issues, let us know and we will do our best to fix them.


As of GOOSE V0.2.0, you should be able to install GOOSE with a single command. 

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
GOOSE has some important limitations that users should be aware of. 

GOOSE makes sequences predicted to be disordered
-------------------------------------------------
GOOSE makes sequences **predicted** to be disordered based on the disorder predictor metapredict. Although modern disorder predictors have proven to be *quite good*, one should aways keep in mind that predicted disorder is **not** gaurenteed disorder. 

Allowed error in sequence properties
-------------------------------------
By default when using the ``create`` functionality, GOOSE allows a *small* amount of error in properties. This is a balance between accuracy and speed. The allowed error is:

* For hydropathy, the allowed error is 0.07. Override by specifying hydropathy_tolerance.
* For kappa, allowed error is 0.03. Override by specifying kappa_tolerance.
* If you choose an NCPR / FCR combination that is mathematically impossible, GOOSE will get as close as it can.

In addition, if you install GOOSE locally, you can go into goose/backend/parameters and modify these values globally.

Specified properties and speed
-------------------------------------
The more properties you specify, or the more constraints you put on sequence design, the more time it will take GOOSE to generate your sequence. 

Failed sequence generation
---------------------------
Sometimes GOOSE can't make your sequence. Here are some tips on getting around this:

- Run the code a few more times. GOOSE often will eventually make your sequence thanks to the inherent stochasticity in sequence generation. 
- If using the ``create`` functionality:

  #. Increase ``attempts``. Default is 100.
  #. Reduce the disorder cutoff value by specifying ``disorder_cutoff``. Default is 0.5.
  #. Increase the allowed error in properties by specifying ``hydropathy_tolerance`` and ``kappa_tolerance``. Default is 0.07 and 0.03, respectively.
  #. *Slightly* adjust your specified properties. 
  #. Try using the ``SequenceOptimizer`` instead as it offers more flexibility.

- If using the ``SequenceOptimizer`` functionality:

  #. Increase ``max_iterations``. Default is 1,000.
  #. Increase the tolerance allowed for each property by specifying the ``tolerance`` argument when defining each property. Default is 0.00.
  #. Try changing the weights of your specified properties. 


Limits on specifying sequence properties
-----------------------------------------
When using the ``create`` functionality, GOOSE will only return sequences with disorder values above the disorder threshold. Some sequence compositions are simply not predicted to be disordered. It is also important to note that some combinations of sequence properties are not mathematically possible. GOOSE uses a rescaled Kyte Doolittle hydropathy scale for calculating mean hydrophobicity. This scale goes from 0 to 9 where higher values are more hydrophobic. The charged residues have low hydrophobicity values (R = 0, K = 0.6, D = 1, E = 1). Therefore, if you have a sequence with too many charged residues, you limit how high the mean hydrophobicity can go. If you specify a high FCR and a high hydrophobicity, that sequence may be mathematically impossible to make. GOOSE will return an error if you do this. 

Best practices when using GOOSE
--------------------------------
It is best practice to double check that the sequences you make using GOOSE are what you intended. You can do this using the *analyze* module included with GOOSE. Although we have done extensive testing on GOOSE functionality, due to the massive sequence space that is possible when generating an IDR, you may encounter bugs. We would appreciate if you would report these bugs, and we will do our best to fix them as quickly as possible.
