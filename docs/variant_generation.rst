  
Creating Sequence Variants in Python
=====================================

First import ``create`` from goose (if you haven't done so for sequence generation).

.. code-block:: python

    from goose import create

Once ``create`` has been imported, you can start making sequence variants!

Apart from simply generating sequences, GOOSE can help you make different types of sequence variants. In contrast to when you generate a sequence, the primary input for the sequence variant functions is your sequence of interest. 

*Disorder cutoffs when creating sequence variants*:

When making sequence variants, by default GOOSE will use the predicted disorder values of your input sequence as the threshold disorder values for the returned sequence. However, you can change this by setting ``strict_disorder=True``, which will make GOOSE use the cutoff disorder value across the entire sequence.

Types of sequence variants
---------------------------

``minimal_var()`` - Variant where the sequence is changed as little as possible while still changing the properties you specify.

``constant_class_var()`` - Variant with the same properties as the input variant as well as the same order and number of amino acids by class. GOOSE will try to change the sequence as much as possible within these constraints.

``new_seq_constant_class_var()`` - Variant where the sequence composition is new but the numbers of each residue from each class and the overall properties are the same.

``constant_properties_var()`` - Variant where **only the sequence properties** are constrained. There are no constraints on classes of amino acids. 

``constant_residue_var()`` - Variant where specific residues are held constant. The variant will have the same aggregate properties as the original sequence.

``region_shuffle_var()`` - Variant that will shuffle specific regions of an IDR. Multiple regions can be specified simultaneously.

``excluded_shuffle_var()`` - Variant where you can specifically shuffle a sequence *except for any specified residues.*

``targeted_shuffle_var()`` - Variant where you specify *which residues are shuffled*. Any residues not specified will not be shuffled. 

``asymmetry_var()`` - Variant where a class of residues (see below for classes) or a user-specified list of residues is changed to become more asymmetrically or less asymmetrically distributed throughout the sequence. Does NOT change sequence composition.

``hydro_class_var()`` - Like the ``constant_class_var()``, properties and the order / number of amino acids by class is held constant. However, hydropathy can be increased or decreased within this constraint. *Note* - because classes of residues are constraints, there are limits to how much you can increase or decrease the hydropathy of any specific sequence.

``fcr_class_var()`` - Function to make a sequence variant that adjusts the FCR while minimizing changes to the position and number of amino acids by class.

``ncpr_class_var()`` - Function to make a sequence variant that adjusts the NCPR while minimizing changes to the position and number of amino acids by class.

``kappa_var()`` - Variant where you can alter the charge asymmetry by changing the kappa value. Requires the presence of positively charged and negatively charged residues in the original sequence. Higher kappa values increase charge asymmetry, lower kappa values reduce charge asymmetry. Values can be between 0 and 1. As mentioned in the sequence generation section, specifying kappa **requires** that you do not have an FCR=0 and that you do not have FCR=NCPR (there needs to be some + and some - charged residues for kappa). Also, GOOSE is much faster at making sequences with kappa between 0.1 and 0.9. Values below 0.1 or above 0.9 may take longer. 

``all_props_class_var()`` - Function to make a sequence variant that adjusts the FCR, NCPR, hydropathy, and kappa values while minimizing changes to the position and number of amino acids by class. If you don't specify one of the values, GOOSE will keep it the same as it was in the input sequence.

``re_var()`` - Function to make a sequence variant that adjusts the Re while keeping amino acid composition constant.

``rg_var()`` - Function to make a sequence variant that adjusts the Rg while keeping amino acid composition constant.

``weighted_shuffle_var`` - Generate variants where you can specify residues or classes of residues to shuffle along with a weight to dictate how severely to shuffle the sequence.

**A note about FCR_class(), NCPR_class(), and all_props_class_var() variants** - 
For the ``fcr_class_var()``, ``ncpr_class_var()``, and ``all_props_class_var()`` variants, the changes to amino acid by class is **MINIMIZED** but not necessarily kept exactly the same. This is because if you (for example) change FCR in your sequence, it is IMPOSSIBLE to keep the order and number of all amino acids by class the same in the returned variant. Similarly, with the NCPR variant, if you change the NCPR to the extent that the FCR has to change as well, then it will change the order / number of amino acids by class.

For some variants, in addition to being able to specify residues using your own custom-defined list, you can specify amino acids by class. The classes are categorized as followed:

``aromatic`` : 'F', 'W', 'Y' 
``polar`` : 'Q', 'N', 'S', 'T' 
``positive`` : 'K', 'R' 
``negative`` : 'D', 'E' 
``hydrophobic``' : 'I', 'V', 'L', 'A', 'M'
``Special Cases`` : 'C', 'P', 'G', and 'H'
The ``Special Cases`` residues are, for any function that accounts for the class of a residue, not interchangable with any other residues. 

The constant_class_var()
------------------------

The ``constant_class_var()`` generates a variant with the same properties as the input variant as well as the same order and number of amino acids by class.

**Example**

.. code-block:: python

    test = 'QEQNGVDQQETTPRQDYPGNQQPNQQAEGQQMQ'
    create.constant_class_var(test)
    'QENQGADQQDQNPRNEWPGNNNPNQTADGNSAT'


The new_seq_constant_class_var()
---------------------------------

The ``new_seq_constant_class_var()`` makes a sequence where the sequence composition is new but the numbers of each residue from each class and the overall properties are the same.

**Example**

.. code-block:: python

    test = 'QGENNENPQDQGSREGPQNNAWAQNNQDAQTSP'
    create.new_seq_constant_class_var(test)
    'QNSAQNDGQNENYQPQGDNPDKNGTSQEAPQAN'


The constant_properties_var()
---------------------------------

The ``constant_properties_var()`` makes a sequence where **only the sequence properties** are constrained.

**Example**

.. code-block:: python

    test = 'QEQNGVDQQETTPRQDYPGNQQPNQQAEGQQMQ'
    create.constant_properties_var(test)
    'TDTGGPDPQDNPTKPENTNQNSGQTQSENSNSN'


The constant_residue_var()
----------------------------

``constant_residue_var()`` - Variant where specific residues are held constant. The variant will have the same aggregate properties as the original sequence. You can specify more than one residue to be held constant at once.

**Example**

.. code-block:: python

    test = 'QEQNGVDQQETTPRQDYPGNQQPNQQAEGQQMQ'
    create.constant_residue_var(test, constant=['T', 'Q'])
    'QDQSMNDQQETTGKQDNAGGQQHPQQPDAQQSQ'


The region_shuffle_var()
--------------------------

``region_shuffle_var()`` - Variant that will shuffle specific regions of an IDR. Multiple regions can be specified simultaneously.
**Note** - The region_shuffle_var does **NOT** use index values like you would normally in Python. For the region_shuffle_var, 1 = the first amino acid in the sequence **NOT 0**. 

**Example with one shuffled region**

.. code-block:: python

    test = 'QQQEEENNNDDDQQQEEENNNDDD'
    create.region_shuffle_var(test, shuffle=[3,9])
    'QQNNQENEEDDDQQQEEENNNDDD'

**Example with two residues constant**

.. code-block:: python

    test = 'QQQEEENNNDDDQQQEEENNNDDD'
    create.region_shuffle_var(test, shuffle=[[3,9], [15, 23]])
    'QQNENEQENDDDQQNEDDQNEEND'

**Notice that when you specify 2 regions, you use a list of lists (a nested list).**

The excluded_shuffle_var()
-----------------------------

``excluded_shuffle_var()`` - Variant where you can specifically shuffle a sequence *except for any specified residues.*

**Example**

.. code-block:: python

    test = 'QQQEEENNNDDDQQQEEENNNDDD'
    create.excluded_shuffle_var(test, exclude_aas=['N', 'D'])
    'EQEEEQNNNDDDQQEQQENNNDDD'

The targeted_shuffle_var()
---------------------------

``targeted_shuffle_var()`` - Variant where you specify *which residues are shuffled*. Any residues not specified will not be shuffled. 

**Example**

.. code-block:: python

    test = 'QQQEEENNNDDDQQQEEENNNDDD'
    create.targeted_shuffle_var(test, target_aas=['N', 'D'])
    'QQQEEENNDNNNQQQEEEDDNDDD'

The asymmetry_var()
---------------------

``asymmetry_var()`` - Variant where a class of residues or a user-specified list of residues is changed to become more asymmetrically or less asymmetrically distributed throughout the sequence. Does NOT change sequence composition.

**Example** - 

**Changing polar residues, no specification of changes property** - 

.. code-block:: python

    test = 'NSQSSQDSQDKSQGSQNQQEQSDSSEQTKQEEDGQTSSDSREQSQSHSQQ'
    create.asymmetry_var(test, 'decrease', 'polar')
    'NSQSQDSQDKSQGQNQQEQSDSSEQTSKQSEEDQGQTSSDSREQSQSHSQ'
    
**Example** - 

**Changing polar residues, increased number of changes** - 

.. code-block:: python

    test='NSQSSQDSQDKSQGSQNQQEQSDSSEQTKQEEDGQTSSDSREQSQSHSQQ'
    create.asymmetry_var(test, 'increase', 'polar', number_changes=30)
    'NQSTQQQSQQSNSTQSSQQQQSSQSSSSQSSSQQDDKGEDEKEEDGDREH'
    

**Changing polar residues, decrease asymmetry** - 

.. code-block:: python

    test='QELQAAAALQQPQTGKSASVQDSALSALQSLLARQSSLSL'
    create.asymmetry_var(test, 'decrease', 'aliphatic', number_changes=30)
    'QELQALQQPQLATGKSLASVAQDSASAQSLARLQSSLASL'
    

**Changing custom list, increase asymmetry** - 

.. code-block:: python

    test='RGNNLAGIVLGAAGAMNGRTEGRKGEQTHGKSGNDDRGHTGDRSHGNKNRGE'
    create.asymmetry_var(test, 'increase', ['G', 'T'], number_changes=20)
    'RNNLAIVLAAAMNRTERKEQHKSNDDRHGGTGGGGGGGGGGTGDRSHNKNRE'
    


The hydro_class_var()
----------------------

``hydro_class_var()`` - Like the ``constant_class_var()``, properties and the order / number of amino acids by class is held constant. However, hydropathy can be increased or decreased within this constraint. *Note* - because classes of residues are constraints, there are limits to how much you can increase or decrease the hydropathy of any specific sequence. If you go past the maximum change, GOOSE will raise an error (see below).

**Example decreasing hydropathy** - 
The starting hydropathy of the sequence below is  2.0272. Let's raise it to around 2.7.

.. code-block:: python

    test = 'GNGGNRAENRTERKGEQTHKSNHNDGARHTDRRRSHDKNAASRE'
    create.hydro_class_var(test, hydropathy=2.7)
    'GTGGTKMETKTEKKGESTHKTSHSDGLKHTDKKKTHDKTLASRE'

**Example where hydropathy is raised higher than possible**

.. code-block:: python

    test = 'GNGGNRAENRTERKGEQTHKSNHNDGARHTDRRRSHDKNAASRE'
    create.hydro_class_var(test, hydropathy=3.7)
    goose.goose_exceptions.GooseInputError:
    Unable to get to objective hydropathy without changing classes of residues.
    For this sequence the lowest possible hydrpathy is 1.611364.
    For this sequence the highest possible hydropathy is 2.834091.


The fcr_class_var()
--------------------

``fcr_class_var()`` - Function to make a sequence variant that adjusts the FCR while minimizing changes to the position and number of amino acids by class.

**Example** - 
The starting FCR of the sequence is 0.0. Let's increase to 0.2.

.. code-block:: python

    test = 'TTGGATSQAGGATHAQSHANSGTQSTSSPQTQGVNTTSANGQHGQATNQS'
    create.fcr_class_var(test, FCR=0.2)
    'TTGGATSQAGGATHAESHARSGTDSTSSPKTQGVETTSAKGDHGKATEKS'


The ncpr_class_var()
---------------------

``ncpr_class_var()`` - Function to make a sequence variant that adjusts the NCPR while minimizing changes to the position and number of amino acids by class.

**Example** - 
The starting NCPR of the sequence is 0.909. Let's lower it to 0.0.

.. code-block:: python

    test = 'GNGGNRAENRTERKGEQTHKSNHNDGARHTDRRRSHDKNAASRE'
    create.ncpr_class_var(test, NCPR=0)
    'GNEGERGENRAENRTDGKQDTKHESRNDHRNEGRAHRTSHNAAS'


The kappa_var()
----------------

``kappa_var()`` - Variant where you can alter the charge asymmetry by changing the kappa value. Requires the presence of positively charged and negatively charged residues in the original sequence. Higher kappa values increase charge asymmetry, lower kappa values reduce charge asymmetry. Values can be between 0 and 1. 

**Example** - 

First we can take something with very symmetrically positions oppositely charged amino acids and increase the kappa value. For reference, the starting kappa value for this 'test' sequence was 0.0012.

.. code-block:: python

    test = 'QNEKRDQNEKRDQNEKRDQNEKRDQNEKRDQN'
    create.kappa_var(test, kappa=0.9)
    'RKKKRRKQRKQNRQNQQNNNDNQDDEEDEDEE'

Now we can take this newly generated and make the charges more moderately symmetrical (something between what we started with and what we made in the previous example).

.. code-block:: python

    previous_variant = 'QNEKRDQNEKRDQNEKRDQNEKRDQNEKRDQN'
    create.kappa_var(previous_variant, kappa=0.15)
    'KRQKRDQREKRDNKEKNDQNEDRDQNENNEQQ'

**Note** - GOOSE will allow deviation from your input kappa value by up to 0.03. This is to keep GOOSE from being extremely slow. If you need something closer to your desired value, you can try generating a few variants. You'll likely quickly get the exact value you want within a few tries.


The all_props_class_var()
---------------------------

The ``all_props_class_var()`` makes a variant sequence that adjusts the FCR, NCPR, kappa, and mean hydropathy while minimizing changes to the order/number of amino acids *by class*. There is only a limited extent to which the NCPR or NCPR can be altered due to the fact that some FCR/hydropathy values are not compatible.

**Example changing all properties** - 
In this example we will change all 4 possible properties.

.. code-block:: python

    test = 'GNGGNRAENRTERKGEQTHKSNHNDGARHTDRRRSHDKNAASRE'
    create.all_props_class_var(test, hydropathy=2.5, FCR=0.23, NCPR=0, kappa=0.1)
    'GSGGTKIESRTEKSGQQTHDSNHNNGAEHTNNKDSHQNNAASQK'


**Example changing 2 properties** - 
In this example we will just change kappa and hydropathy.

.. code-block:: python

    test = 'GNGGNRAENRTERKGEQTHKSNHNDGARHTDRRRSHDKNAASRE'
    create.all_props_class_var(test, kappa=0.3, hydropathy=2.6)
    'KTGGTKRGSKTARKGKSTHTTKHDEGVRTHDRRLSHEENADSTE'


The re_var() and rg_var
---------------------------

The ``re_var()`` and ``rg_var()`` let you increase or decrease the Re / Rg of your sequence while holding amino acid composition constant. You can choose to just get something that maximally increases or decreases the Re / Rg or you can choose to get back a series of sequence that have increasingly altered Re / Rg from the starting sequence. You need to specify your sequence and ``decrease`` to decrease the Rg / Re or ``increase`` to increase the Rg / Re. The sequence we will start with for the examples below has an 'Rg' = 12.6429 and 'Re' = 19.8837.

**Example chagning Re** - 

.. code-block:: python

    test = 'FYFLGQGQQYYYYQQKQFFQFYYQQFFGFYGSNFQGGNYFGGYQQNQYFG'
    create.re_var(test, 'increase')
    {'QGQGKGFGQQYGQYYNFQFSYFFGYFGFFQYNYYFYFQQQQGYYNFGQQL': 26.680789338243116}

In this example, we increased the Re from 19.8837Å to 26.680789Å. Now let's decrease the Re. 

.. code-block:: python

    test = 'FYFLGQGQQYYYYQQKQFFQFYYQQFFGFYGSNFQGGNYFGGYQQNQYFG'
    create.re_var(test, 'decrease')
    {'YQYYYNYKQYFYGNQGQQFGGQYGYYLNFFGGFFFGQGQFQYQSQFFQQF': 15.126193729754828}

**Example chagning Rg** - 

.. code-block:: python

    test = 'FYFLGQGQQYYYYQQKQFFQFYYQQFFGFYGSNFQGGNYFGGYQQNQYFG'
    create.rg_var(test, 'increase')
    {'FFQKQFGNQGQYQGQQLQGYQYFQGGNGQFSFNQYGFYYYYQFYFFYYGF': 15.116732605201102}

In this example, we increased the Rg from 12.6429Å to 15.1167Å. Now let's decrease the Rg. 

.. code-block:: python

    test = 'FYFLGQGQQYYYYQQKQFFQFYYQQFFGFYGSNFQGGNYFGGYQQNQYFG'
    create.rg_var(test, 'decrease')
    {'LQYQYYYQFGSQYFFNYGQGFFFFQGQFGKQFGGYYGYYFQQNQGNFYQQ': 11.65882476485573}


**Additional usage**  

*Getting variants that span a range of Re or Rg values*  

In addition to getting a single sequence variant when increasing or decreasing the Rg / Re of your sequence, you can also get numerous sequences that span a range of Rg / Re values above or below your sequence. To do this, set ``return_all=True`` in the ``create.re_var()`` or ``create.rg_var()`` function. 

**Example getting variants spanning dimensions for Re or Rg** - 

.. code-block:: python

    test = 'FYFLGQGQQYYYYQQKQFFQFYYQQFFGFYGSNFQGGNYFGGYQQNQYFG'
    create.re_var(test, 'increase', return_all=True)
    {'FGGNNQFYFKGYYQFYGQFFQYYFSQQGLYGFQNGFYQYQGFYFGQQQYQ': 19.88582890272039, 'NQFGSYFFGFFFGQYQQYYGQQGGYGFQLQNGYYQGFQQFNFFYYKQYYQ': 20.09158138778087, 'FYQKFGFYYNQQQQLQGFQGGGGYYYFFQFQQQFQFYGYGYFSNFNQYYG': 20.29236054471707, 'GFQFQNYQYYGYYQQGFYYQFKQYGFQYGGFFLFGQYQGFNNQFQYQGSF': 20.493062151452005, 'NYFGQYQYYQQKGFYFSQYFGQFFFFQYFFQYYQLGQFQGNQGQGNGGYY': 20.69843194312801, 'FGQKGYNNQFGFFYYFGQYQLYGYGFNYYGQQGQGFQFSFQYQQQQFFYY': 20.899782611330014, 'QFYQFGYYQFQQFGFNGYQYYYYGGFNQFYKFNQQLQFQQGFYQSYFGGG': 21.10184471833488, 'QYYFQQFGGFQFKFQFLFQNFYQNQYYQYFQGGGQGYYGFYQQFYGYGNS': 21.302967793554927, 'YQYFQQYYFFGQQQQFFGQYGGYFQLYYFYFNKNNYGQYQGFFGQFSGQG': 21.503399660459397, 'YQYFFNQGYGNGQFFQGQFFGYFQQGYYYYYKGQFQFGSNYGQQYFQQFL': 21.708127134164094, 'GFGNYQGGKYFQFQFFQGYQFQQYQLGGNYGFNYFYQSYQFQFQFYQYYG': 21.913070399733165, 'GNYGFYYGQQFFYYQYNFQFQYGQGQFFFQGYQNSFFFQQGKLYYGYQQG': 22.11329490414296, 'GLFYQQYYGFGQQQYQFYQFSGFGNFYQFGGYNQYGFQFYNYQGQFFQYK': 22.313844782223256, 'LQFGYNKFSGFYQQYFYFYFYQFNGYFQQGGYYQYFGNGQGFQQYQQQFG': 22.51790127809966, 'FQGYGFQGGGYNSYFFQNKQGFQQYQFGQQYYYYYFFQFNFGYFQLYQQG': 22.718152756492504, 'FQGGFFFQYYFFGYFFGQFYGNGYYLQFQYGNQYQQQQQYFKNQYSYGQG': 22.923480401319996, 'QQKYGYQFQGYFFYYQFGFFFGQQNFFQNQYGLGYFNYFGQQQGYYYSGQ': 23.124772063934085, 'SYQYFNQYGGYGFNQQQFQNFGYQFKQQYFGFFYGFFYGYQYQGQFGYQL': 23.333739510598978, 'GSQNFYQYQLYYFYYGYGFYGNFNFQQQFYGFQQYQGFGFQYQKQFGQFG': 23.54112054988592, 'KQQGQFGGNFFYYQQYGFQYGYYQLQQGYFYNYGFSFQFQGYNQFYFFQG': 23.747530525783173, 'QKGFGQYQQYFGFFQYLNQYFFYQFYYFGNYNGGGQYFQQYFGFGQYSQQ': 23.98253798133204, 'QFSNFGQGQFYYQYFQFNFFGQKYFQQFYQYGFQGQYGGYNQYQFGYYGL': 24.190898513378375, 'GNFFKQGYQYGYYYGQGNFGFYQFYYQGFQFNQYFQYYQFFFQGGQQSLQ': 24.39256949762906, 'GQQNFQQYFFFQKSYFYFQGQYGFQFYQGFYYQGGFNYQYQGQYNGYGFL': 24.629380838294193, 'QYSQGQNQNFQFGQGYFYGQQYFYGGQFYKFYFYFGQYFGYFYNQFQGQL': 24.839275515626326, 'KQNFQFFNLGFYFYYYFFSYQGFQFQQGYGQYYYYQFGGFQQYNQQGGGQ': 25.04550510301195, 'GKFGYYGNYQQYFFFYYFYSQQFYQGFGFFGYQYNQGQGYQQFQFNQGLQ': 25.327400084592217, 'QQGKQFFQGGQQGNGQYYFFYQGYQQGYNYQGYQFQYFYYFFGFFYNSFL': 25.643400295987494, 'KQGGYGQGSQNQFQFQGQYGFYFGNFFNQFQFYQYYYYFQGGFFQQYYYL': 26.023861583369275, 'KNQGFGYGFFYYNYFQQGFFYQFFQGFSYFGYFYGYQQQYQYGQQQQLGN': 27.021238093522594}
    

Now let's do the same for Rg.

.. code-block:: python

    test = 'FYFLGQGQQYYYYQQKQFFQFYYQQFFGFYGSNFQGGNYFGGYQQNQYFG'
    create.rg_var(test, 'increase', return_all=True)
    {'FYFLGQGQQYYYYQQKQFFQFYYQQFFGFYGSNFQGGNYFGGYQQNQYFG': 12.642883713816577, 'YFGFFYQQGQFQQQFGYSNFQNGFLQYGQQKQGFYQFQYYYFGGGYYNFY': 12.843408303787768, 'YFGQQYFLQQYFQYFYFNFSYQYQKQQGGYFGGFFQGQQQGGFYYNGNYF': 13.047272610034963, 'QFGQLGGNYNQQFQYKQYQFGFYYQFFNYFYQSQGFFQYQYFGYGFQGYG': 13.249197318313685, 'GQFGKYYNGGGFNQQFYQQYYSGFQQQNGFYFGQFFYQQYQLGQFYYYFF': 13.452808743469815, 'NNYYQQNQGYGQYYQLFFQQQFGGFQFGGQFGYKGQQFFGYFSYYYYFQF': 13.663894490740866, 'GFQQGGFQQQYGGFSNQYKFQFFQYQQFQGNFGQYYNQYGYFGYLFYFYY': 13.868289004342328, 'QYGGQGQKNFYFGGYQQQQFQQYFFQFGNQSGQLFFFQFYYFYYYGYGYN': 14.118385874595063, 'QQYKQQNQYQQGGLGFYQGQQYYQFQGFFYFSQGGGYNFYFFNFYYFYFG': 14.330510120213422, 'FGFLGFFQGQKNQYQQFYGFQGQNFGQGQGQQFQGYSYYQFYYYYFYFNY': 14.801244899483716, 'NFQGNLQFYGQQQGQQNGGQQFQFGYKFQYQGGYGFYQYSYYYYFFFFFY': 15.209855224048953}


*Predicting the Rg / Re of your starting sequence simultaneously*  
GOOSE also lets you get the predicted Rg / Re of your starting sequence. When you set ``include_original=True`` for the ``create.re_var()`` or ``create.rg_var()`` functions, you will get back a dictionary where there with **original** and **variants** as the keys that correspond to values that are dictionaries with sequnce : Rg/Re pairs. You can also set ``return_all=True`` when using this to get back your original sequence and a series of variants with Rg / Re values.  


**Example including original sequence Re or Rg** - 


.. code-block:: python

    test = 'FYFLGQGQQYYYYQQKQFFQFYYQQFFGFYGSNFQGGNYFGGYQQNQYFG'
    create.re_var(test, 'increase', include_original=True)
    {'original': {'FYFLGQGQQYYYYQQKQFFQFYYQQFFGFYGSNFQGGNYFGGYQQNQYFG': 19.883686156942094}, 'variants': {'QKGFQGGQYQQQQFGFYFFYNYFYQNQFQNYYQFFYGFYYGGQGSGYQFL': 26.436047646562177}}


Now let's do the same for Rg. 

.. code-block:: python

    test = 'FYFLGQGQQYYYYQQKQFFQFYYQQFFGFYGSNFQGGNYFGGYQQNQYFG'
    create.rg_var(test, 'increase', include_original=True)
    {'original': {'FYFLGQGQQYYYYQQKQFFQFYYQQFFGFYGSNFQGGNYFGGYQQNQYFG': 12.642882870879607}, 'variants': {'QQNFYQYSGGFQFQKQQQQFNFGGQGFQGFGQFFGYQYYYNYYYLYGFYF': 14.667587127612029}}

The weighted_shuffle_var
---------------------------
``weighted_shuffle_var`` - Generate variants where you can specify residues or classes of residues to shuffle along with a weight to dictate how severely to shuffle the sequence.
**Example** -
In this example we will shuffle the sequence with a weight of 0.5.
.. code-block:: python

    test = 'QQQEEENNNDDDQQQEEENNNDDD'
    create.weighted_shuffle_var(test, shuffle_weight=0.5, target_aas=['Q', 'E'])
    'EQQQEENNNDDDEEQQEQNNNDDD'
