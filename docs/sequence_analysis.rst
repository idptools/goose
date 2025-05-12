
Sequence analysis in GOOSE
===========================

GOOSE provides powerful sequence analysis tools. These tools are provided by SPARROW - see https://github.com/idptools/sparrow/tree/main. The predictions are **all machine learning based**. They should be treated **as predictions** and NOT ground truth. They are **NOT** equivalent to experimental validation.
Because GOOSE relies entirely on SPARROW, we generally recommend just using SPARROW. However, we provide some analysis tools for your convenience.

To start using the analysis tools in GOOSE, first import **analyze**

.. code-block:: python

    from goose import analyze

and then you have a full suite of sequence anlysis tools. The analyze module includes tools for calculating and predicting the following characterstics of your protein of interest (which doesn't even need to be an IDR made in GOOSE!). 

1. Sequence length
2. Fraction of charged residues (FCR)
3. Net charge per residue (NCPR)
4. Average hydropathy
5. kappa - A measurement of charge asymmetry where the max value (1) is the greatest possible charge asymmetry for a given sequence and the min value (0) is the most symmetrical positions possible for oppositely charged residues.
6. The fractions of all amino acids in the sequence.
7. Predicted phosphosites for S, T, and Y phosphorylation sites
8. Predicted cellular localization signals including those for nuclear localization, nuclear export and mitochondrial targeting.
9. Predicted transcriptional activation domains
10. Predicted radius of gyration (Rg) and predicted end-to-end distance (Re)
11. Helical regions.


Using the analyze module in GOOSE
----------------------------------
**Analyzing general properties**

To analyze general properties in GOOSE, you can use the 'properties' function.
This function returns a dict containing all basic properties calculated including length, FCR, NCPR, hydropathy, kappa, and fractions of amino acids.

**Example**

.. code-block:: python

    test = 'GNGGNRAENRTERKGEQTHKSNHNDGARHTDRRRSHDKNAASRE'
    analyze.properties(test)
    {'length': 44, 'FCR': 0.4090909090909091, 'NCPR': 0.09090909090909091, 'hydropathy': 2.027272727272727, 'kappa': 0.050873519852757454, 'fractions': {'A': 0.09091, 'C': 0.0, 'D': 0.06818, 'E': 0.09091, 'F': 0.0, 'G': 0.11364, 'H': 0.09091, 'I': 0.0, 'K': 0.06818, 'L': 0.0, 'M': 0.0, 'N': 0.13636, 'P': 0.0, 'Q': 0.02273, 'R': 0.18182, 'S': 0.06818, 'T': 0.06818, 'V': 0.0, 'W': 0.0, 'Y': 0.0}}


**Predicting phosphosites**

GOOSE has 3 separate networks each trained on a different phosphorylation sites (S, T, and Y phosphosites). The phosphosites() function in analyze gives you information on all of them at once. 

**Example**

.. code-block:: python

    test = 'STAYELANSYPSSYGTTLYKSTSYTEDSPTGYYTTYVSRQ'
    analyze.phosphosites(test)
    {'S': [37], 'T': [21, 24, 29, 33, 34], 'Y': [3, 13, 18, 23]}

**It is extremely important to note that this predictor is not a gaurentee of a phosphorylation event!**. Protein phospohrylation is incredibly complex, this predictor should be used more as a way to check on something that you want to avoid being phosphorylated (although as with any 'predictor', nothing can be gaurenteed 100%).


**Predicting subcellular localization**

If you design an IDR and it ends up somewhere you don't want it to, that's a bad day in the lab. To try to mitigate this problem, we are working on machine learning-based predictors of cellular localization sequences to try to determine where a given protein might end up. As stated in the phosphosite section, this is not a perfect solution, but it is nonetheless better than nothing. 

**Example**

For a known mitochondrial protein...

.. code-block:: python

    test = 'MAAAAASLRGVVLGPRGAGLPGARARGLLCSARPGQLPLRTPQAVALSSKSGLSRGRKVMLSALGMLAAGGAGLAMALHS'
    analyze.cellular_localization(test)
    {'mitochondria': {'MAAAAASLRGVVLGPRGAGLPGARARGLLCSARPGQLPLRTPQAVALSSKSGLSRGRKVMLSAL': [1, 65]}, 'NES': 'No NES sequences predicted.', 'NLS': 'No NLS targeting sequences predicted.'}

In the returned dict, the key 'mitochondria' brings up the sequence predicted to be the targeting sequence as well as the coordinates of that sequence (where 1 is the first amino acid). For the other two keys, 'NES' for nuclear export sequences and 'NLS' for nuclear localization sequences, because none were detected the value for each of those key value pairs just states that none were predicted.


**Predicting transcriptional activation domains**

If you design an IDR, you might not to want it to inadvertantly have a transcriptional activation domain (TAD). To see if your sequence has a predicte TAD...

**Example**

For a subset of a protein with a known TAD...

.. code-block:: python

    test = 'PNNLNEKLRNQLNSDTNSYSNSISNSNSNSTGNLNSSYFNSLNIDSMLDDYVSSDLLLNDDDDDTNLSR'
    analyze.transcriptional_activation(test)
    {'TGNLNSSYFNSLNIDSML': [31, 49]}
    
If there is a TAD present, the function returns the TAD subsequence along with the coordinates for the TAD in the input sequence. 

**Predict everything**

If you just want a summary of... well basically everything we've covered so far from properties to all predicted features, that's pretty easy to do! 

**Example**

.. code-block:: python

    test = 'PNNLNEKLRNQLNSDTNSYSNSISNSNSNSTGNLNSSYFNSLNIDSMLDDYVSSDLLLNDDDDDTNLSR'
    analyze.everything(test)
    {'length': 69, 'FCR': 0.2029, 'NCPR': -0.11594, 'hydropathy': 3.41304, 'kappa': 0.4222, 'fractions': {'A': 0.0, 'C': 0.0, 'D': 0.14493, 'E': 0.01449, 'F': 0.01449, 'G': 0.01449, 'H': 0.0, 'I': 0.02899, 'K': 0.01449, 'L': 0.14493, 'M': 0.01449, 'N': 0.23188, 'P': 0.01449, 'Q': 0.01449, 'R': 0.02899, 'S': 0.21739, 'T': 0.04348, 'V': 0.01449, 'W': 0.0, 'Y': 0.04348}, 'helical regions': [[3, 13], [44, 50]], 'predicted phosphosites': {'S': [67], 'T': [30, 64], 'Y': [18, 37, 50]}, 'predicted cellular localization': {'mitochondria': 'No mitochondrial targeting sequences predicted.', 'NES': 'No NES sequences predicted.', 'NLS': 'No NLS targeting sequences predicted.'}, 'predicted transcriptional activation': {'TGNLNSSYFNSLNIDSML': [31, 49]}, 'predicted polymer properties': {'Rg': 23.2913, 'Re': 55.2608}, 'fraction aromatic': 0.05797, 'fraction polar': 0.50725, 'fraction aliphatic': 0.2029, 'sequence': 'PNNLNEKLRNQLNSDTNSYSNSISNSNSNSTGNLNSSYFNSLNIDSMLDDYVSSDLLLNDDDDDTNLSR'}

The analyze.everything() function will return a dictionary holding all of the information from sequence properties to predicted phosphosites, cellular localization, and transcriptional activation domains all from one simple function!

**Predict differences between sequences**

If you generate a sequence variant and want to see if you've broken or introduced any sequence features including TADs, cellular localization signalgs, and phosphosites, you can use the ``analyze.prediction_diffs()`` function. The function takes in two sequences as the input and then returns the predicted differences between those sequences. 

**Example**

.. code-block:: python

    test1 = 'PNNLNEKLRNQLNSDTNSYSNSISNSNSNSTGNLNSSYFNSLNIDSMLDDYVSSDLLLNDDDDDTNLSR'
    test2 = 'PTTITEKIKTTITTDTTTFTTTITTTTTTSTGNLNSSYFNSLNIDSMLDDYVSSDLLLNDDDDDTNLSR'
    analyze.prediction_diffs(test1, test2)
    {'S phosphorylation': 'No differences.', 'T phosphorylation': 'sequence 1: [30, 64], sequence 2: [64]', 'Y phosphorylation': 'sequence 1: [18, 37, 50], sequence 2: [37, 50]', 'NLS': 'No differences.', 'NES': 'No differences.', 'mitochondrial': 'No differences.', 'predicted transcriptional activation': ['Sequence 1 predicted TAD - TGNLNSSYFNSLNIDSML : [31, 49] not in sequence 2', 'Sequence 2 predicted TAD - GNLNSSYF : [32, 40] not in sequence 1']}


