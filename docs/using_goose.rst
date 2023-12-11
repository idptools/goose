Leaving the nest - how to use GOOSE in Python
==============================================
To use GOOSE from in Python, first import *create from goose*

.. code-block:: python

    from goose import create

Once *create* has been imported, you can start making sequences and sequence variants!

Generating sequences with specified properties
===============================================
GOOSE can generate sequences by either specifying sequence properties or fractions of amino acids. 

The ``create.sequence()`` function lets you create sequences predicted to be disordered with various specified properties. 

The only required argument is the length, which must be between 10 and 10,000. In addition, you can also specify several parameters.

1. ``hydropathy``: Average hydropathy, which must be between 0.1 and 6.1. It is called `hydropathy` instead of `hydrophobicity` because we are using a rescaled Kyte-Doolittle hydropathy scale.

2. ``FCR``: The fraction of charged residues, which must be between 0 and 1.

3. ``NCPR``: The net charge per residue (`NCPR`), which must be between -1 and +1.

4. ``kappa``: The kappa value is a charge asymmetry parameter where higher values mean greater charge asymmetry. This must be between 0 and 1. NOTE: kappa values of 0 or 1 can take a long time to make because there are very few sequences that satisfy these values. GOOSE is much faster at making sequences between 0.1 and 0.9. Also, specifying kappa **requires** that you do not have an FCR=0 and that you do not have FCR=NCPR (there needs to be some + and some - charged residues for kappa).

**Importantly, you do not need to specify all of these sequence properties simultaneously.** For example, if you specify FCR and hydropathy, GOOSE will return sequences that have varying NCPR and kappa values while making sure that the specified hydropathy and FCR values are what you input. In this way, you can generate many sequences that have the fixed properties that you want to stay fixed, while other properties vary.

In addition to the above parameters, you can specify:

1. ``cutoff``:  The disorder cutoff defines the value used to accept/reject sequences. A higher cutoff means higher confidence that the sequence is disordered. The default value is 0.5. If you have difficulty making your sequence, you might want to try lowering the cutoff value. We use Metapredict V2-FF in GOOSE (see https://github.com/idptools/metapredict).

2. ``attempts``: The number of attempts defines how many times GOOSE will try to generate a desired sequence. This is relevant because certain parameter combinations of values are more challenging than others, and GOOSE implements a stochastic design algorithm such that different sequences are generated every time. For harder sequence compositions you may need to increase this from the default value. 

3. ``exclude``: A list of residues to exclude from your generated sequence. **NOTE** - If you specify FCR or NCPR, you CANNOT have any charged residues in the ``exclude`` list. Additionally, **you can not specify more than 10 residues to be excluded**. 

4. ``no_constraints``: A boolean value that lets you generate sequences without any constraints on sequence properties. Constraints on sequence properties are only added if ``kappa`` is specified and ``NCPR`` or ``FCR`` are not specified. The reason for these constraints is to make sure that there are oppositely charged residues in the sequence generated. In other words, FCR and NCPR are chosen to be random values that are result in at least some oppositely charged residues to be in the sequence generated (but the FCR and NCPR values are **still randomly chosen**, just within a more constrained space than 0-1). To bypass this, you can set ``no_constraints=True`` when generated a sequence; although, this is generally not recommended. 


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

Hydropathy values can be between 0.0 and 6.1. **Note**: the higher the hydropathy *over 5*, the longer it will take GOOSE to generate the sequence. Sequences that are very hydrophobic and disordered can be tricky to make. **Note**: whenever specifying hydropathy values, GOOSE will return a sequence within 0.07 of the specified value! This amount of error helps keep GOOSE fast (and a difference of 0.07 is typically negligible).

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


**Kappa with no_constraints set to True**

.. code-block:: python

    create.sequence(100, kappa=0.2, no_constraints=True)
    'THQQSEEAQASSQTTSEDGKQSHEGHEASGAKNESHGHHKQSNGTKHGDTRTHDTQNKTTAQSHRGDENRKKEGNDDEGAHAADDAHPAHSGTRQHQTKH'



**Note** - Generating this sequence fails frequently. To bypass this I increased the number of attempts by specifying ``attempts=1000``. It should be noted that creation of this sequence still fails occassionally even when attempts is increased to 1000. 

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

Note that if you pass in requested fractions, those fractions must be equal to or less than 1.

In addition to specifying the specific amino acid fractions, other parameters can be passed to the `create.seq_fractions()` function:

1. ``cutoff``:  The disorder cutoff used defines a threshold GOOSE uses to accept/reject sequences. 

2. ``attempts``: The number of attempts defines how many times GOOSE will try to generate a desired sequence. 

3. ``max_aa_fractions``: If you wish to generate sequences with extreme compositions it may be necessary to over-ride the default max fractional values. This can be achieved by passing a max_aa_fractions dictionary, which should specify key-value pairs for amino acid-max fraction information. 


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

Examples of generating sequences by specifying Rg or Re
----------------------------------------------------------

**Specifying a length and Rg:**

.. code-block:: python

    create.seq_rg(50, 20)
    'PDAESQCNTSRWIVSHPQSNTKYPDSRTDESASPQQEDPSHSQEKIHSRM'


**Specifying a length and Re:**

.. code-block:: python

    create.seq_re(50, 20)
    'FYFLGQGQQYYYYQQKQFFQFYYQQFFGFYGSNFQGGNYFGGYQQNQYFG'

Creating Sequence Variants in Python
=====================================

Apart from simply generating sequences, GOOSE can help you make different types of sequence variants. In contrast to when you generate a sequence, the primary input for the sequence variant functions is your sequence of interest. 

*Disorder cutoffs when creating sequence variants*:

When making sequence variants, by default GOOSE will use the predicted disorder values of your input sequence as the threshold disorder values for the returned sequence. However, you can change this by setting ``strict_disorder=True``, which will make GOOSE use the cutoff disorder value across the entire sequence.

Types of sequence variants
---------------------------

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

**A note about FCR_class(), NCPR_class(), and all_props_class_var() variants** - 
For the ``fcr_class_var()``, ``ncpr_class_var()``, and ``all_props_class_var()`` variants, the changes to amino acid by class is **MINIMIZED** but not necessarily kept exactly the same. This is because if you (for example) change FCR in your sequence, it is IMPOSSIBLE to keep the order and number of all amino acids by class the same in the returned variant. Similarly, with the NCPR variant, if you change the NCPR to the extent that the FCR has to change as well, then it will change the order / number of amino acids by class.

For all some variants, in addition to being able to specify residues using your own custom-defined list, you can specify amino acids by class. The classes are categorized as followed:

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
    NDNQGAENNDQNPRNEYPGQSSPNQNADGNNAS


The new_seq_constant_class_var()
---------------------------------

The ``new_seq_constant_class_var()`` makes a sequence where the sequence composition is new but the numbers of each residue from each class and the overall properties are the same.

**Example**

.. code-block:: python

    test = 'QEQNGVDQQETTPRQDYPGNQQPNQQAEGQQMQ'
    create.new_seq_constant_class_var(test)
    QNSAQNDGQNENYQPQGDNPDKNGTSQEAPQAN


The constant_properties_var()
---------------------------------

The ``constant_properties_var()`` makes a sequence where **only the sequence properties** are constrained.

**Example**

.. code-block:: python

    test = 'QEQNGVDQQETTPRQDYPGNQQPNQQAEGQQMQ'
    create.constant_properties_var(test)
    MEHPTHDQYDQNHKQEPTGSNPNGTPHETNPQP


The constant_residue_var()
----------------------------

``constant_residue_var()`` - Variant where specific residues are held constant. The variant will have the same aggregate properties as the original sequence. You can specify more than one residue to be held constant at once.

**Example**

.. code-block:: python

    test = 'QEQNGVDQQETTPRQDYPGNQQPNQQAEGQQMQ'
    create.constant_residue_var(test, constant=['T', 'Q'])
    QEQSANDQQETTPKQEAPSPQQASQQHEGQQPQ


The region_shuffle_var()
-------------------

``region_shuffle_var()`` - Variant that will shuffle specific regions of an IDR. Multiple regions can be specified simultaneously.
**Note** - The region_shuffle_var does **NOT** use index values like you would normally in Python. For the region_shuffle_var, 1 = the first amino acid in the sequence **NOT 0**. 

**Example with one shuffled region**

.. code-block:: python

    test = 'QQQEEENNNDDDQQQEEENNNDDD'
    create.region_shuffle_var(test, shuffle=[3,9])
    QQNQENNEEDDDQQQEEENNNDDD

**Example with two residues constant**

.. code-block:: python

    test = 'QQQEEENNNDDDQQQEEENNNDDD'
    create.region_shuffle_var(test, shuffle=[[3,9], [15, 23]])
    QQNNEEQNEDDDQQNDENNEDEQD

**Notice that when you specify 2 regions, you use a list of lists (a nested list).**

The excluded_shuffle_var()
-----------------------------

``excluded_shuffle_var()`` - Variant where you can specifically shuffle a sequence *except for any specified residues.*

**Example**

.. code-block:: python

    test = 'QQQEEENNNDDDQQQEEENNNDDD'
    create.excluded_shuffle_var(test, exclude_aas=['N', 'D'])
    QEQEEQNNNDDDQQEQEENNNDDD

The targeted_shuffle_var()
---------------------------

``targeted_shuffle_var()`` - Variant where you specify *which residues are shuffled*. Any residues not specified will not be shuffled. 

**Example**

.. code-block:: python

    test = 'QQQEEENNNDDDQQQEEENNNDDD'
    create.targeted_shuffle_var(test, target_aas=['N', 'D'])
    QQQEEENNDNNNQQQEEEDDNDDD

The asymmetry_var()
---------------------

``asymmetry_var()`` - Variant where a class of residues or a user-specified list of residues is changed to become more asymmetrically or less asymmetrically distributed throughout the sequence. Does NOT change sequence composition.

**Example** - 

**Changing polar residues, no specification of changes property** - 

.. code-block:: python

    test = 'NSQSSQDSQDKSQGSQNQQEQSDSSEQTKQEEDGQTSSDSREQSQSHSQQ'
    create.asymmetry_var(test, 'decrease', 'polar')
    NSQSSQDSQDKSQGSQNQEQSDSSEQTQKQEEDGQTSSDSREQSQSHSQQ
    
**Example** - 

**Changing polar residues, increased number of changes** - 

.. code-block:: python

    test='NSQSSQDSQDKSQGSQNQQEQSDSSEQTKQEEDGQTSSDSREQSQSHSQQ'
    create.asymmetry_var(test, 'increase', 'polar', number_changes=30)
    SNSQQQSSNTQSSSQSSQSSSTQQSSQQQQSQQQDDKGEDEKEEDGDREH
    

**Changing polar residues, decrease asymmetry** - 

.. code-block:: python

    test='SNSQQQSSNTQSSSQSSQSSSTQQSSQQQQSQQQDDKGEDEKEEDGDREH'
    create.asymmetry_var(test, 'decrease', 'aliphatic', number_changes=30)
    QNQSTSQQSDDQQKSTNGSSEDSQQESQKQESSEQSSDQGSQDSSQREQH
    

**Changing custom list, increase asymmetry** - 

.. code-block:: python

    test='RGNNLAGIVLGAAGAMNGRTEGRKGEQTHGKSGNDDRGHTGDRSHGNKNRGE'
    create.asymmetry_var(test, 'increase', ['G'], number_changes=20)
    RNNLAIVLGGGGGGGGGGGGGAAAMNRTERKEQTHKSNDDRHTDRSHNKNRE
    


The hydro_class_var()
----------------------

``hydro_class_var()`` - Like the ``constant_class_var()``, properties and the order / number of amino acids by class is held constant. However, hydropathy can be increased or decreased within this constraint. *Note* - because classes of residues are constraints, there are limits to how much you can increase or decrease the hydropathy of any specific sequence. If you go past the maximum change, GOOSE will raise an error where the error message specifies the minimum and maximum possible values for your sequence (see below).

**Example decreasing hydropathy** - 
The starting hydropathy of the sequence below is  2.0272. Let's raise it to around 2.7.

.. code-block:: python

    test = 'GNGGNRAENRTERKGEQTHKSNHNDGARHTDRRRSHDKNAASRE'
    create.hydro_class_var(test, hydropathy=2.7)
    GTGGTKIETKTEKKGETTHKTTHTDGLKHTDKKKTHDKSAASRE

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
    TTGGATSQAGGATHAESHARSGTDSTSSPKTQGVETTSAKGDHGKATEKS


The ncpr_class_var()
---------------------

``ncpr_class_var()`` - Function to make a sequence variant that adjusts the NCPR while minimizing changes to the position and number of amino acids by class.

**Example** - 
The starting NCPR of the sequence is 0.909. Let's lower it to 0.0.

.. code-block:: python

    test = 'GNGGNRAENRTERKGEQTHKSNHNDGARHTDRRRSHDKNAASRE'
    create.ncpr_class_var(test, NCPR=0)
    GNEGERGENRAENRTDGKQDTKHESRNDHRNEGRAHRTSHNAAS


The kappa_var()
----------------

``kappa_var()`` - Variant where you can alter the charge asymmetry by changing the kappa value. Requires the presence of positively charged and negatively charged residues in the original sequence. Higher kappa values increase charge asymmetry, lower kappa values reduce charge asymmetry. Values can be between 0 and 1. 

**Example** - 

First we can take something with very symmetrically positions oppositely charged amino acids and increase the kappa value. For reference, the starting kappa value for this 'test' sequence was 0.0012.

.. code-block:: python

    test = 'QNEKRDQNEKRDQNEKRDQNEKRDQNEKRDQN'
    create.kappa_var(test, kappa=0.9)
    KKRRKKRRRKQNQNQNQNEEEDQDENEDQDDN

Now we can take this newly generated and make the charges more moderately symmetrical (something between what we started with and what we made in the previous example).

.. code-block:: python

    previous_variant = 'KKRRKKRRRKQNQNQNQNEEEDQDENEDQDDN'
    create.kappa_var(previous_variant, kappa=0.15)
    QEEEDRDKEKERDRDRDKNQNQNQNQNKKRQN

**Note** - GOOSE will allow deviation from your input kappa value by up to 0.03. This is to keep GOOSE from being extremely slow. If you need something closer to your desired value, you can try generating a few variants. You'll likely quickly get the exact value you want within a few tries.


The all_props_class_var()
---------------------------

The ``all_props_class_var()`` makes a variant sequence that adjusts the FCR, NCPR, kappa, and mean hydropathy while minimizing changes to the order/number of amino acids *by class*. There is only a limited extent to which the NCPR or NCPR can be altered due to the fact that some FCR/hydropathy values are not compatible.

**Example changing all properties** - 
In this example we will change all 4 possible properties.

.. code-block:: python

    test = 'GNGGNRAENRTERKGEQTHKSNHNDGARHTDRRRSHDKNAASRE'
    create.all_props_class_var(test, hydropathy=2.5, FCR=0.23, NCPR=0, kappa=0.1)
    GSGGTKIESRTEKSGQQTHDSNHNNGAEHTNNKDSHQNNAASQK


**Example changing 2 properties** - 
In this example we will just change kappa and hydropathy.

.. code-block:: python

    test = 'GNGGNRAENRTERKGEQTHKSNHNDGARHTDRRRSHDKNAASRE'
    create.all_props_class_var(test, kappa=0.3, hydropathy=2.6)
    GTKGGSKMTKTKKGKESTHKRTSHKSRDEKGVHTDSHEDNAASE


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


Sequence analysis in GOOSE
===========================

GOOSE provides powerful sequence analysis tools. These tools are provided by SPARROW - see https://github.com/idptools/sparrow/tree/main. The predictions are **all machine learning based**. They should be treated **as predictions** and NOT ground truth. They are **NOT** equivalent to experimental validation.

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

    test = 'GNGGNRAENRTSSKSERKGEQTHKSNHNDGARHTDRRRSHYDKNAASRE'
    analyze.phosphosites(test)
    {'S': [11, 12, 14], 'T': [21, 33], 'Y': [40]}

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


Generating sequence libraries in GOOSE
========================================
Apart from generating individual sequences, you can also generate libraries of synthetic sequences with predicted disorder that span user-defined ranges of sequence properties or fractions of amino acids. 

Generating sequence libraries by properties
--------------------------------------------

To generate sequence libraries by specifying properties, use the ``create.seq_property_library()`` function. An important things to note **GOOSE automatically gets rid of sequences not possible to make**. For example, a sequence with an NCPR value where the absolute value of the NCPR value is greater than the FCR will not be included because it's not possible to generate. In addition, GOOSE will check the hydropathy value and see if it is possible based on the FCR. If it's not possible, GOOSE will not make it. However **for incompatible hydropathy / charge values**, GOOSE will print out the sequences it was not able to generate and print them unless you set ``silent_failed_seqs`` to True. In this function, you can specify the length, FCR, NCPR, hydropathy, kappa, and disorder cutoff values. **Note** - GOOSE will do it's best to get the kappa value spot on, but it does allow for some error. It will also also adjust the FCR to match the NCPR if the two values are not compatible. This function will return a dictionary of sequences where each sequence is named after the property specified for generating the sequence.

 
For this function, you can specify the property as a single value, a list of two values, or a list with three values where GOOSE will make all possible values from the first number of the list to the second number of the list at an interval equal to the third value of the list.

**Example**

.. code-block:: python

    create.seq_property_library(40, FCR=0.1, NCPR=[0, 0.1], hydropathy=[2, 4, 0.5])

    {'>FCR_0.1_NCPR_0_hydropathy_2.0': 'QKPSQNKNHTPTGQGNSHPQDHPEQQQQQPPQQQSTQNTP', '>FCR_0.1_NCPR_0_hydropathy_2.5': 'NSNSTSENNKQNNGPHSPGTSQPPNFQPSAPQENSGNGKH', '>FCR_0.1_NCPR_0_hydropathy_3.0': 'SAPTQDPQSHYTQGGNEQTGGSPTGPPGWSHAKRSPSGQG', '>FCR_0.1_NCPR_0_hydropathy_3.5': 'NRSSGSCAPLNSAGGTTPGNKEVADPPPPGSTGSWGHQTH', '>FCR_0.1_NCPR_0_hydropathy_4': 'PSTHSSAGPSDTSASSSARSVPSSDSAVKSSCGSGASTTS', '>FCR_0.1_NCPR_0.1_hydropathy_2.0': 'QPPSPHQPLSHHSQQHNGNTKKQKSHQPNKNNSHPNNHNQ', '>FCR_0.1_NCPR_0.1_hydropathy_2.5': 'SGHSQGQNTHTKQGRQRGHGHVSPNQQHSSTPQHMQSPKT', '>FCR_0.1_NCPR_0.1_hydropathy_3.0': 'TSPSNHPQKPGPTPAGMQTGGTPGKTKHPHHPGSKLQQYT', '>FCR_0.1_NCPR_0.1_hydropathy_3.5': 'SAMLNASAGNPSGGQQRNSANLGPSRTTQKTSAQARSPTG', '>FCR_0.1_NCPR_0.1_hydropathy_4': 'PAPKPGAKVVSTSALQRVAKSSPPACSPGTHPGSSPTTSS'}

In the example above, the first value ``40`` in the function is the length. This value must be a single integer value and is **required**. The second value, ``FCR=0.1`` is specified as a single value, so all sequences generated will have that value. The third value ``NCPR=[0, 0.1]`` is specified as a list, so GOOSE made a set of sequences where NCPR was equal to the first value in the list and then a set where NCPR was equal to the second value in the list. Finally, ``hydropathy=[2, 4, 0.5]`` was set equal to a list with 3 values. This means GOOSE was told to generate sequences with a range of the first value of the list ``2``, to the second value in the list, ``4``, at an interval equal to the third value of the list ``0.5``. This resulted in the creation of sequences with values 2.0, 2.5, 3.0, 3.5, and 4.0.
**If the third value, which is the interval value, cannot be equally used between the range of sequences, GOOSE will just use the maximum value as the last value.** For example, hydropathy = [2, 4, 1.1] would result in hydropathy values of 2, 3.1, and 4.

**Additional Usage** - 
In addition, you can add a random name to the end of each sequence name by setting ``random_name=True``. This way, if you combine multiple libraries (if you want replicates of sequences with the varying properties you specify), you don't need to worry about overwriting anything due to a shared name. 

**Example**

.. code-block:: python

    create.seq_property_library(40, FCR=0.1, NCPR=[0, 0.1], hydropathy=[2, 4, 0.5], random_name=True)

    {'>FCR_0.1_NCPR_0_hydropathy_2.0_F0K7D2N5N6': 'TPNHTQPHKNHDNNNPSHNHTGNNPTNPQEHKGSQTNQPT', '>FCR_0.1_NCPR_0_hydropathy_2.5_G5R2A2L8F8': 'THQNPEDTHTTHPSMSRSNNPQLQNNGQRPAPPSSPHGHN', '>FCR_0.1_NCPR_0_hydropathy_3.0_E8V4W5C2N2': 'GPGSEHPHAPGDSSTGNNTSGPTKPSTGGALSQNRQPQYP', '>FCR_0.1_NCPR_0_hydropathy_3.5_C2C3D1F5P4': 'PPTQQPNGQSMSGGARHTTAAAAEGSARMAELSQNHSNGG', '>FCR_0.1_NCPR_0_hydropathy_4_I3S8I0G9W4': 'APTKGVAPETRSTSPAASSGAGGGGSSPASSMSPSSGDGS', '>FCR_0.1_NCPR_0.1_hydropathy_2.0_N8E1Y2R5D0': 'PPPSTGHQKQNYSQNHHNNPPQHQRWHRNGPPRPNSHQSG', '>FCR_0.1_NCPR_0.1_hydropathy_2.5_H3H7I2S0W7': 'NSGGSKRSSSPGPTPNQPQNGRNMPMPQNRQNHTNFQNTP', '>FCR_0.1_NCPR_0.1_hydropathy_3.0_I2K0M6E2T0': 'TRSHQQPQMHGMPSTSPNGCQTLNSPSSMRKGPPPQSGKN', '>FCR_0.1_NCPR_0.1_hydropathy_3.5_H3Q4E1L7W8': 'QTSSPQTMGRSQTTTGSASMQSSGMASTSRPPRFSSRSTG', '>FCR_0.1_NCPR_0.1_hydropathy_4_A4D7R8W6V7': 'RTSPTSVKPPTSACKTAAGSTPMTRSPSSSTLAVNGPPAP'}

Generating sequence libraries by fractions
-------------------------------------------

To generate sequence libraries by specifying fractions, use the ``create.seq_fractions_library()`` function. An important thing to note ***GOOSE automatically gets rid of sequences not possible to make***. This includes where the fraction of amino acids is greater than 1 (for obvious reasons) and if the fraction of any amino acid specified is greater than the max_aa_fractions limit. The max_aa_fractions limit can be manually overridden (see below). 

For this function, you can specify the amino as a single value, a list of two values, or a list with three values where GOOSE will make all possible values from the first number of the list to the second number of the list at an interval equal to the third value of the list. If only two values are given, GOOSE just uses those two values and does not assume any interval.

**Example**

.. code-block:: python

    create.seq_fractions_library(100, A=0.1, D=[0.1, 0.2], K=[0.1, 0.2, 0.05])

    {'>A_0.1_D_0.1_K_0.1': 'SKSSSGTKEHSEAEDNGEDGAATNHDHNDEHGRATGADDKKNHGHTKEHGAQHQSSQGNNNHDKSNSTRDAHNGARSDKRARNKEKQHQKGQAGENDHGE', '>A_0.1_D_0.1_K_0.15': 'EAAGGHQHGKRQSGKSQSADENKGRKKDESTNKDNTSRQSRETASQGKAKKNNGGPKKAGNQDDAQDESEGSQRSSQQAKDAKNGDDQTKTDEGHSTKAQ', '>A_0.1_D_0.1_K_0.2': 'TGGAGKDASAGDATKDRAKSDNKGKTKKERAAKTQNKSHNQAQEKRTGESSHKEKRKDGENQAKSHSKNHKQRKPADTQKTEDERHEEHGHEDKKDEDEQ', '>A_0.1_D_0.2_K_0.1': 'DSSNATTSNDQDDKSHDSNTQAHDQREVGNSKDSNSNASDDENKAGQENTSAEEDNPDDHEEKDDNDRDGHAKKKTSADKDGDNDREAKASHNGKNAEEG', '>A_0.1_D_0.2_K_0.15': 'DQGDEKAQTDASKSKNDTGAADKHGKAKQTGKEEENQDDGKTDDHSGPTDGQDNRGDKKSEGTDDKAKDQQDDTDEQATTTTKRAGHAADEDSTNRKKRS', '>A_0.1_D_0.2_K_0.2': 'SDDKDRRDKQAHSNKHADAKSNEASHRKKHAGKHGQDTGKKDDGQNKDSADKTHKTKDGDSEQKAHDTSEADQAKKDGDHNGEDGDEDGDKAGKKQGNKN'}

In the example above, the first value ``100`` in the function is the length. This value must be a single integer value and is **required**. The second value, ``A=0.1`` is specified as a single value, so all sequences generated will have that value. The third value ``D=[0.1, 0.2]`` is specified as a two item list, so GOOSE made a set of sequences where the fraction of **D** is either 0.1 or 0.2. Finally, ``K=[0.1, 0.2, 0.05]`` was set equal to a list with 3 values. This means GOOSE was told to generate sequences with a range of the first value of the list ``0.1``, to the second value in the list, ``0.2``, at an interval equal to the third value of the list ``0.05``. This resulted in the creation of sequences with values 0.1, 0.15, 0.2.
**If the third value, which is the interval value, cannot be equally used between the range of values from the lowest to the highest, GOOSE will just use start at the lowest value and increase until it can't any more and then will add the maximum value.** For example, ``K = [0.1, 0.2, 0.08]`` would result in K fraction values of , 0.1, 0.18, and 0.2.
In addition, you can add a random name to the end of each sequence name by setting ``random_name=True``. 

**Additional usage**

Some other things you can specify are:

``warn_user`` - used to determine whether to warn the user of any problems with generating the sequnces. By default is set to True, but you can set it to False. 
``robust_warnings`` - used to return a warning message FOR EVERY SEQUENCE THAT HAS A FRACTION VALUE NOT EQUAL TO THE INPUT VALUE. This can be annoying for large library generation and is by default set to False. Set to True to get more information on the sequences that have errors and what those errors are. 
``max_aa_fractions`` - used to override the max amino acid fractions. Input as a dict. Example below:

**Example**

.. code-block:: python

    create.seq_fractions_library(100, A=0.1, D=[0.1, 0.2], K=[0.1, 0.2, 0.05], max_aa_fractions= {'K': 0.16})

    {'>A_0.1_D_0.1_K_0.1': 'TQHHDEKNRRAEANDSPNGDEHAQDGKHSAEKQRTQHAENSDRDHSEGAKGNNQHGRKKQENERRAGGGKQTHKNTSQDHGDRNAKDDAQNGQQHHNKHA', '>A_0.1_D_0.1_K_0.15': 'QRKGNSANSGERADHTGDHDTQNAATTGKRKDNKEKKNDKHSARAQNTDKKAHTEKGSKHATQNAHNESQPGGDETNSKKHASTGQGKGNNDSRKGRDDN', '>A_0.1_D_0.2_K_0.1': 'AQQRQDQGGDAKDADDRTDARKDTETSPAKQEQAGRSDDKGPDDTDQKAESPTESNERDQQGQETGDDKQQKKGSEDAHSQDQDGKGPDQGKDAHQAGSR', '>A_0.1_D_0.2_K_0.15': 'EGPQATSTDDEDDHHKSKDESADEGAKSGKRTEENDRAATDTAHHATKDHDDHHKTDGPEKTKDETKKADKEGHHKKDTAEKQEDNANSSTDDTQPSKDD'}

In the above example, we manually overrode the max fraction for K and set it to 0.16. This eliminated sequences where the K fraction was 0.2 like in the example above where the max_aa_Fractions were left as default.


Copyright (c) 2023, Ryan Emenecker - Holehouse Lab