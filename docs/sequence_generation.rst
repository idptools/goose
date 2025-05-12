Creating sequences and sequence variants with GOOSE
=====================================================
Both sequence generation and variant generation are in the ``create`` moduel. First import *create from goose*

.. code-block:: python

    from goose import create

Once *create* has been imported, you can start making sequences and sequence variants!

Generating sequences with specified properties
===============================================
GOOSE can generate sequences by either specifying sequence properties or fractions of amino acids. 

The ``create.sequence()`` function lets you create sequences predicted to be disordered with various specified properties. 

The only required argument is the length, which must be between 10 and 10,000. In addition, you can also specify several parameters.

1. ``hydropathy``: Average hydropathy, which must be between 0.1 and 6.1. It is called `hydropathy` instead of `hydrophobicity` because we are using a rescaled Kyte-Doolittle hydropathy scale.

2. ``FCR``: The fraction of charged residues (`FCR`), which must be between 0 and 1.

3. ``NCPR``: The net charge per residue (`NCPR`), which must be between -1 and +1.

4. ``kappa``: The kappa value is a charge asymmetry parameter where higher values mean greater charge asymmetry. This must be between 0 and 1. NOTE: kappa values of 0 or 1 can take a long time to make because there are very few sequences that satisfy these values. GOOSE is much faster at making sequences between 0.1 and 0.9. Also, specifying kappa **requires** that you do not have an FCR=0 and that you do not have FCR=NCPR (there needs to be some + and some - charged residues for kappa).

**Importantly, you do not need to specify all of these sequence properties simultaneously.** For example, if you specify FCR and hydropathy, GOOSE will return sequences that have varying NCPR and kappa values while making sure that the specified hydropathy and FCR values are what you input. In this way, you can generate many sequences that have the fixed properties that you want to stay fixed, while other properties vary.

In addition to the above parameters, you can specify:

1. ``cutoff``:  The disorder cutoff defines the value used to accept/reject sequences. A higher cutoff means higher confidence that the sequence is disordered. The default value is 0.5. We use Metapredict V3 in GOOSE (see https://github.com/idptools/metapredict).

2. ``attempts``: The number of attempts defines how many times GOOSE will try to generate a desired sequence. For more complex sequence compositions you may need to increase this from the default value. 

3. ``exclude``: A list of residues to exclude from your generated sequence. **NOTE** - If you specify FCR or NCPR, you CANNOT have any charged residues in the ``exclude`` list. Additionally, **you can not specify more than 10 residues to be excluded**. 

4. ``use_weighted_probabilities``: Whether to use weighted probabilities for generating sequences. Default is False. The TL;DR is that we have generated many thousands of sequences that fit different parameter combinations and created weighted probabilities for each residue to increase the chance of generating the objective sequence. This was necessary in older versions of GOOSE. However, we have implemented Numpy vectorization to speed up sequence generation, and the weighted probabilities are not necessary anymore (though can still speed things up). However, if you want to use them, you can set this to True.

5. ``strict_disorder``: Whether you want **all** residues to be above the cutoff value. By default, GOOSE lets a small number (at most 5 for very large sequences) of residues be below 0.5 because a single amino acid (or even a few) below the 0.5 threshold is not realistically going to be a folded domain.

6. ``return_all_sequences``: Whether to return all sequences generated. Default is False. This is new as of v0.2.0 because GOOSE generates and checks many sequences simultaneously, and it is now possible to return all of them that match the specified parameters.

7. ``custom_probabilities``: You can now specify your own custom probabilities for generating sequences. This is a dictionary where the keys are the amino acids and the values are the probabilities. The sum of all probabilities must equal 1.

8. ``metapredict_version``: You can specify the version of metapredict you want to use. The default is 3, but you can set it to 2 or 1 if you want to use the older versions.


Examples of sequence generation by properties
----------------------------------------------

Just specifying length:

.. code-block:: python

    create.sequence(40)
    'GDHNKAGQPPRKCSDQGGAGAPNPDCDPDTAPMDGDRMTN'


Specifying length and hydropathy:

.. code-block:: python

    create.sequence(100, hydropathy = 3)
    'MTSYGRDGSPETGEGSTGTNSSSSRSMMGSTHNWQQYNGGTTSGTSSTGDSHRTHGDHSAGETTSGGDSEGTDETSTTTNGRGSSSGHDGSTGQDTNTRR'

Hydropathy values can be between 0.0 and 6.1. Higher values can take slightly longer to generate. 

**Specifying length and fraction of charged residues (FCR):**

.. code-block:: python

    create.sequence(40, FCR = 0.3)
    'GDRPSEHGQGPRKEDGMDQDDVSTEGHEWSNNPCNQSNNP'

FCR values can be between 0 and 1

**Specifying length and the net charge per residue (NCPR):**

.. code-block:: python

    create.sequence(40, NCPR = -0.2)
    'MQKNDRAPDHKDREKDGPIKERPEECPDDEQSDDEECPSH'

NCPR values can be between -1 and 1.

 
**Specifying multiple properties**

GOOSE lets you combine different properties simultaneously. Importantly, any value you do not specify will just be random.

**Examples**

**FCR & NCPR**

.. code-block:: python

    create.sequence(100, FCR = 0.3, NCPR = -0.1)
    'TSNQDKEMPQQHSPRCQPGEKVSDPPRSSDNSTNGGARPQQDWRPPEHMNPNRYEPNTMHQNREGRESAGGKDWPNPTIDQNQDPHEDTDNQEEESDHPC'

You cannot have values for NCPR where the absolute value of NCPR is greater than the specified FCR value. 

**Important note on combining FCR and NCPR!** Whenever NCPR and FCR are combined, if the combinations of the length, NCPR, and FCR are not mathematically possible, GOOSE will get as close as it can. In addition, GOOSE  prioritizes NCPR over FCR, and the resulting sequence may deviate in terms of FCR as a result.

**FCR & Hydropathy**

.. code-block:: python

    create.sequence(100, FCR = 0.3, hydropathy = 3.2)
    'KVDSGTTSCSGERESDSGDLKSSKEGSSGSGSSSKSSKSKEATGSSTDTTAAAGGKGGGGGGDGGKGDGRGKGGGGGGEGRDGGGGGGEGGRGGGGRKRD'

**Note** - The higher FCR (or NCPR because the absolute value of NCPR must be at least equal to FCR) specified, the lower the maximum possible hydropathy because charged residues have a *very* low hydropathy value.

**NCPR & Hydropathy**

.. code-block:: python

    create.sequence(100, NCPR = -0.3, hydropathy = 2.4)
    'REARGDAKGERDRGGDAKDKGAESGKDDDGEEEGAGEEEGEEGDDEAEADRADKERAERDKGDRDRAEGRAEKGAAAAEGADEGADEADEEEDDDADDEE


**NCPR, FCR, & Hydropathy**

.. code-block:: python

    create.sequence(100, hydropathy = 2.65, NCPR = -0.3, FCR = 0.4)
    'NETPARPETHRDTASTSEGDETSEPEGTWSSNEADTDDDAETEHSPMSEDGERCESSKDAPPMRDEEGDDEDVEDTPDVSSSPDYEPGGHYSESNNDWPD'


**NCPR, FCR, Hydropathy, and kappa**

.. code-block:: python

    create.sequence(100, hydropathy = 2.65, NCPR = 0.0, FCR = 0.4, kappa=0.2)
    'GKDETATKRQKAPPVDRREAPAKHKRTTAGRRDRSPKEKETRMGQGGPEGESPSSGGDETEGIMARKASEDSTPGKMNSSRDRSDGEHGETPPVEPDPNH'


**Hydropathy, FCR, NCPR, excluding values, and increasing attempt number**

.. code-block:: python

    create.sequence(100, FCR=0.6, NCPR=0.6, hydropathy=3, exclude=['C'], attempts=1000)
    'VSKKLKAKIKSPKRKRKKKKLKVKARSRKRAKLSVVKRKRMSVKVAKRSKVRAFMVRRKKKPKPFKRKVKAVRKKKRRPKKKRIAKKRVKKVKRKRKKVI'

**Specifying custom probabilities**

.. code-block:: python

    create.sequence(100, hydropathy=2.5, custom_probabilities={'A':0.1, 'R':0.1, 'D':0.1, 'E':0.1, 'G':0.1, 'H':0.1, 'I':0.1, 'K':0.1, 'L':0.1, 'M':0.1})
    'IIPKLKDKIRKIGEEDGRKKKRAIRKRVKRRRKIRCRMKMEDERKRAARLRSKRKKDDGKHKTAKKKRERRKLKRYRRLELRGKDKGDDKIVRKKDMIDP'

**Specifying metapredict version**
.. code-block:: python

    create.sequence(100, hydropathy=2.5, metapredict_version=2)
    'RDRFSEYKNTKEQAFDSYQLERHKERESQTRKRHRPQREKQRPDGERHKHEFMEWKLERRRCTEDGDKEFRLQALGRCESPIGMQMHTPDIADPKRDRRN'


Generating Sequences specifying Fractions of Amino Acids
=========================================================

The ``create.seq_fractions()`` function lets you create sequences predicted to be disordered with specified fractions of various amino acids. With this function, you can specify multiple amino acids simultaneously. Each fraction should be specified using a decimal value (for example, if you want one-tenth of the amino acids to be alanine use ``A=0.1``).

For each amino acid, we had GOOSE attempt (at least 10,000 times for each value) to make sequences with increasing fractions of each amino acid until we identified the maximum possible fraction. The default maximum values for each amino acid are as follows - 

.. code-block:: python

    "A" - 0 : 0.95, 
    "R" - 0 : 1.0, 
    "N" - 0 : 1.0, 
    "D" - 0 : 1.0, 
    "C" - 0 : 1.0, 
    "Q" - 0 : 1.0, 
    "E" - 0 : 1.0, 
    "G" - 0 : 1.0, 
    "H" - 0 : 1.0, 
    "I" - 0 : 0.53, 
    "L" - 0 : 0.42, 
    "K" - 0 : 1.0, 
    "M" - 0 : 0.62, 
    "F" - 0 : 1.0, 
    "P" - 0 : 1.0, 
    "S" - 0 : 1.0, 
    "T" - 0 : 1.0, 
    "W" - 0 : 0.55, 
    "Y" - 0 : 0.99, 
    "V" - 0 : 0.71

Note that if you pass in requested fractions, those fractions cannot be greater than 1. Any values that are remaining will be randomly added. 

In addition to specifying the specific amino acid fractions, other parameters can be passed to the `create.seq_fractions()` function:

1. ``cutoff``:  The disorder cutoff used defines a threshold GOOSE uses to accept/reject sequences. 

2. ``attempts``: The number of attempts defines how many times GOOSE will try to generate a desired sequence. 

3. ``max_aa_fractions``: If you wish to generate sequences with extreme compositions it may be necessary to over-ride the default max fractional values. This can be achieved by passing a max_aa_fractions dictionary, which should specify key-value pairs for amino acid-max fraction information. 

4. ``strict_disorder``: Whether you want **all** residues to be above the cutoff value. By default, GOOSE lets a small number (at most 5 for very large sequences) of residues be below 0.5 because a single amino acid (or even a few) below the 0.5 threshold is not realistically going to be a folded domain.

5. ``return_all_sequences``: Whether to return all sequences generated. Default is False. This is new as of v0.2.0 because GOOSE generates and checks many sequences simultaneously, and it is now possible to return all of them that match the specified parameters.

6. ``metapredict_version``: You can specify the version of metapredict you want to use. The default is 3, but you can set it to 2 or 1 if you want to use the older versions.


Examples of Sequence Generation by Fractions
---------------------------------------------

**Specifying a single amino acid fraction:**

.. code-block:: python

    create.seq_fractions(100, Q=0.3)
    'QEQNGVDQQETTPRQDYPGNQQPNQQAEGQQMQSTKMHDQHDSVNEDQEQNQNPWGHQPHMKGESNSSAREAQSEDQQNQAQNQQQNHDSTQQQDGQMDQ'

**Specifying multiple amino acids:**

.. code-block:: python

    create.seq_fractions(100, Q=0.3, S=0.3, E=0.1)
    'QEQQSQKASQSQVESQDSSESSAPGSSQMHQQQSQSQEGMEQHQSSVGNSSSYPQSEQSEQQRQQSSQDQQQQSSSQTSEENSQSRQHDMSDTEMSGSQR'

**Note** - 
Some combinations of amino acids are simply not possible to make that are predicted to be disordered using the default settings. Specifically, specifying high fractions of multiple aliphatics or aromatics may not be predicted to be disordered using the default cutoff value.

**Excluding a specific amino acids:**
If you want to exclude an amino acid, you can set it equal to 0.

.. code-block:: python

    create.seq_fractions(50, A=0)
    'NKERPTGSWDEPPFDEGSSGMTNEDMGNKPYPTTDMQPEKWPQNDQQGST'
    

**Overriding default max fractions:**  

.. code-block:: python

    create.seq_fractions(100, Y=0.5, max_aa_fractions={'Y':1}) 
    'SSYYYYYSYSSYYSYSSGHYYSYSSYYYSSSYYSSYGGTYGYYSYSYGYYSSYYYSYSSNYYYYYYYYSSYGNSGYGGYYSYYSSSQHHYSSYYYSYYSY'
 



Generating Sequences specifying Ensemble Dimensions
=========================================================

The ``create.seq_rg()`` and ``create.seq_re()``functions let you create sequences with a specified length and a predicted radius of gyration (Rg) or end-to-end distance (Re). For these functions, you must specify the length and an objective Re or Rg. In addition you can also specify:

1. ``cutoff``:  The disorder cutoff used defines a threshold GOOSE uses to accept/reject sequences. 

2. ``attempts``: The number of attempts defines how many times GOOSE will try to generate a desired sequence. 

3. ``strict_disorder``: Whether you want **all** residues to be above the cutoff value. By default, GOOSE lets a small number (at most 5 for very large sequences) of residues be below 0.5 because a single amino acid (or even a few) below the 0.5 threshold is not realistically going to be a folded domain. 

4. ``exclude_aas``: A list of residues to exclude from your generated sequence. There are some limitations on excluding AAs, specifically you can't simultaneously exclude  W, Y, G, F, Q, and N or D, E, K, P, S, and T. 

5. ``allowed_error``: How far off from your desired Re/Rg in Å GOOSE can be before returning the sequence. A higher value here will decrease the time it takes GOOSE to make the sequence. 

6. ``reduce_pos_charged``: Whether to reduce positively charged amino acids in the sequence. Default is True. The reason for this is that in vivo data suggests that positively charged residues may not drive sequence expansion as much as was predicted by the model used here for predicted rg / re. Therefore, when set to True, this function will largely avoid high numbers of (+) charged residues if possible. 

Examples of generating sequences by specifying Rg or Re
----------------------------------------------------------

**Specifying a length and Rg:**

.. code-block:: python

    create.seq_rg(50, 20)
    'NSETSEFYNDPVNAQPGDDHNSENNSVTYDNTGTYSNEFPDTEPSDLHAP'


**Specifying a length and Re:**

.. code-block:: python

    create.seq_re(50, 20)
    'FGQQGGQWGQWGNGQWGYWQNFGYGGNGGWYFYQWYNWFQYNWWFWQWWF'

  
**Specifying a length and Rg, allowing positive charged residues:**

.. code-block:: python

    create.seq_rg(50, 20, reduce_pos_charged=False)
    'NQKDSPEIDKPKPGNASGKFQTIRGNNRRKQKGGQGYPEKTIGERHMSEA'
  


Copyright (c) 2023, Ryan Emenecker - Holehouse Lab