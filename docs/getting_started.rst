What is GOOSE?
===============
**GOOSE : Generate disOrdered prOtiens Specifying propErties** is a python package developed to generate intrinsically disordered proteins or protein regions (collectively IDRs) and IDR variants. Basically, we want to make it easy for you to design IDRs or IDR variants for your research. My goal is to make this *as useful as possible to the protein disorder community*. If you have any feature requests, please email me at remenecker@wustl.edu with the subject line **GOOSE feature request** and I'll be happy to see if I can make it happen. 

Why was GOOSE made?
--------------------
It is difficult to know which properties of any given intrinsically disordered protein region (IDR) are important for function. Is it the hydropathy? Maybe the fraction of charged residues? Who knows! Further, generating IDR variants can be a tricky task. This is why we made GOOSE. GOOSE lets you **1.** generate *synthetic* IDRs to study how different sequence properties or fractions of amino acids impact some aspect of disordered proteins and **2.** create IDR variants based on an IDR that you are interested in so you can figure out which aspects of your IDR of interest are important for function. 

How can I use GOOSE?
--------------------
You can use GOOSE from Python or from a Google Colab notebook. The Colab notebook can be found at https://colab.research.google.com/drive/1U9B-TfoNEZbbjhPUG5lrMPS0JL0nDB3o?usp=sharing

Installation - GOOSE takes flight!
===================================
Right now you can only install GOOSE through Github. It will be on PyPi to allow for pip installation soon!

To install directly from the git repository simply run:

.. code-block:: bash

    $ pip install git+https://github.com/idptools/goose.git

Or to clone the GitHub repository and install locally run - 

.. code-block:: bash

    $ git clone https://github.com/idptools/goose.git
    $ cd goose
    $ pip install .

**Important note**

GOOSE also requires the package ``sparrow``. Sparrow should be downloaded automatically just by pip installing GOOSE, but if you have issues, try installing it by running:

.. code-block:: bash

    $ pip install git+https://github.com/holehouse-lab/sparrow.git

This will install SPARROW. If your attempted install of SPARROW fails, it may be because you do not have numpy or cython installed. I made them both required for installation of GOOSE, so if you install GOOSE first, you should be ok. If you still have problems, try running from terminal:

.. code-block:: bash

    $ pip install cython
    $ pip install numpy


What can GOOSE do?
===================
There are four main functionalities currently in GOOSE. 

**1.** Generatae IDRs where you can specify the *length*, *average hydrophobicity*, *fraction of charged residues (FCR)*, *net charge per residue (NCPR)*, and *kappa*, which is a property that defines opposite charge distribution in a sequence.

**2.** Generate IDR variants. The different variants I currently have implemented should allow you to address many questions you might be interested in by letting you choose which parts or properties of your IDR to hold constant and which to change. This can really narrow down *what makes your IDR of interest work*. 

**3.** Make sequence libraries. This includes libraries of sequences spanning different sequence properties or fractions of amino acids. 

**4.** Sequence analysis. Analysis includes sequence properties, fractions of amino acids, and various predicted sequence features. Predicted sequence features are all machine learning based and are generated using SPARROW (see https://github.com/idptools/sparrow). Predictions include phosphorylation sites, cellular localization, radius of gyration, end-to-end distance, and transcriptional activation domains. By using any predictions in GOOSE you acknowledge that they are just predictions and not ground truth and the accuracy of any predictions generated cannot be gauranteed.

Important limitations
======================
GOOSE has some important limitations that users should be aware of. First, GOOSE makes sequences **predicted** to be disordered based on the disorder predictor metapredict V2. Although modern disorder predictors have proven to be *quite good*, one should aways keep in mind that predicted disorder is **not** gaurenteed disorder. 

Allowed error in sequence properties
-------------------------------------
GOOSE by default allows a *small* amount of error between some user input properties and properties of returned sequences. For hydropathy, the allowed error is 0.07, which is honestly negligible. For kappa, allowed error is 0.03. This is a balance between accuracy and speed. If you install GOOSE locally, you can go into goose/backend/properties and modify these values globally. Finally, if you choose an NCPR / FCR combination that is mathematically impossible, GOOSE will get as close as it can.

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

