Creating Sequence Variants in Python
=====================================

First import ``create`` from goose (if you haven't done so for sequence generation).

.. code-block:: python

    from goose import create

Once ``create`` has been imported, you can start making sequence variants!

Apart from simply generating sequences, GOOSE can help you make different types of sequence variants. The primary input for sequence variant generation is your sequence of interest, and you specify the type of variant you want to create.

Overview of the variant() function
----------------------------------

GOOSE provides a unified interface for generating sequence variants through the ``variant()`` function. This function takes your input sequence, a variant type, and any additional parameters needed for that specific variant type.

.. code-block:: python

    variant_sequence = create.variant(sequence, variant_type, **kwargs)

*Disorder cutoffs when creating sequence variants*:

When making sequence variants, by default GOOSE will use the predicted disorder values of your input sequence as the threshold disorder values for the returned sequence. However, you can change this by setting ``strict_disorder=True``, which will make GOOSE use the cutoff disorder value across the entire sequence.

Types of sequence variants
---------------------------

The ``variant()`` function supports multiple variant types, each with specific parameters and behaviors:

**Shuffling methods:**
- ``'shuffle_specific_regions'`` - Shuffle only specified regions
- ``'shuffle_except_specific_regions'`` - Shuffle all except specified regions  
- ``'shuffle_specific_residues'`` - Shuffle only specific residue types
- ``'shuffle_except_specific_residues'`` - Shuffle all except specific residue types
- ``'weighted_shuffle_specific_residues'`` - Weighted shuffle of specific residues
- ``'targeted_reposition_specific_residues'`` - Reposition specific residues

**Residue asymmetry methods:**
- ``'change_residue_asymmetry'`` - Change residue asymmetry patterns

**Property methods:**
- ``'constant_properties'`` - Generate variant with constant properties (NCPR, FCR, hydropathy, and kappa)
- ``'constant_residues_and_properties'`` - Keep specified residues and properties constant. The sequence generated will have the same properties as the input sequence, but with specified residues kept constant. 
- ``'constant_properties_and_class'`` - Generate variant with constant properties and the number of amino acids by each amino acid class
- ``'constant_properties_and_class_by_order'`` - Generate variant with constant properties and the number and order of amino acids by class constant
**Property modification methods:**
- ``'change_hydropathy_constant_class'`` - Change hydropathy while keeping class constant
- ``'change_fcr_minimize_class_changes'`` - Change FCR while minimizing changes to amino acid classes. Prioritizes keeping aromatics constant then H, C, and P, then aliphatics, then polar.
- ``'change_ncpr_constant_class'`` - Change NCPR while keeping class constant
- ``'change_kappa'`` - Change kappa value. Sequence composition stays constant. 
- ``'change_properties_minimize_differences'`` - Change properties while minimizing differences. This function is a little bit slower because it tries to change the fewest residues possible to achieve the desired properties.
- ``'change_any_properties'`` - Change any combination of properties. Similar to change_properties_minimize_differences, but changes are not necessarily minimized.
- ``'change_dimensions'`` - Change sequence dimensions (Rg/Re). This allows changes in the sequence including the amino acids by class.

Common parameters
-----------------

Most variant types support these common parameters:

- ``num_attempts`` (int): Number of attempts to generate variant (default: 100)
- ``strict_disorder`` (bool): Whether to use strict disorder checking (default: False)
- ``disorder_cutoff`` (float): Disorder cutoff threshold (default: from parameters)
- ``metapredict_version`` (int): MetaPredict version to use (default: 3)
- ``hydropathy_tolerance`` (float): Hydropathy tolerance (default: from parameters) (only if hydropathy is a factor)
- ``kappa_tolerance`` (float): Kappa tolerance (default: from parameters) (only if kappa is a factor)

For some variants, you can specify amino acids by class. The classes are categorized as follows:

- ``aromatic``: 'F', 'W', 'Y' 
- ``polar``: 'Q', 'N', 'S', 'T' 
- ``positive``: 'K', 'R' 
- ``negative``: 'D', 'E' 
- ``hydrophobic``: 'I', 'V', 'L', 'A', 'M'
- ``cystine``: 'C'
- ``proline``: 'P'
- ``glycine``: 'G'
- ``histidine``: 'H'

The ``Special Cases`` residues are, for any function that accounts for the class of a residue, not interchangeable with any other residues.

Shuffling variants
------------------

Shuffle specific regions
~~~~~~~~~~~~~~~~~~~~~~~~

The ``'shuffle_specific_regions'`` variant type shuffles only specified regions of the sequence.

**Parameters:**
- ``shuffle_regions`` (list): List of tuples specifying (start, end) positions to shuffle

**Example:**

.. code-block:: python

    test = 'QQQEEENNNDDDQQQEEENNNDDD'
    variant_seq = create.variant(test, 'shuffle_specific_regions', 
                                shuffle_regions=[(2, 9), (14, 22)])
    print(variant_seq)
    # Output: 'QQEEQENNNDDDQQNQNENEDEDD'

**Note:** Region specifications use 0-based indexing where (start, end) includes positions from start to end-1, following Python slice conventions.

Shuffle except specific regions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``'shuffle_except_specific_regions'`` variant type shuffles all regions except those specified.

**Parameters:**
- ``excluded_regions`` (list): List of tuples specifying (start, end) positions to exclude from shuffling

**Example:**

.. code-block:: python

    test = 'QQQEEENNNDDDQQQEEENNNDDD'
    variant_seq = create.variant(test, 'shuffle_except_specific_regions',
                                excluded_regions=[(0, 5), (18, 24)])
    print(variant_seq)
    # Output: 'QQQEENQEDENQDENDEQNNNDDD'

Shuffle specific residues
~~~~~~~~~~~~~~~~~~~~~~~~~

The ``'shuffle_specific_residues'`` variant type shuffles only specific residue types.

**Parameters:**
- ``target_residues`` (list): List of residue types to shuffle

**Example:**

.. code-block:: python

    test = 'QQQEEENNNDDDQQQEEENNNDDD'
    variant_seq = create.variant(test, 'shuffle_specific_residues',
                                target_residues=['N', 'D'])
    print(variant_seq)
    # Output: 'QQQEEENNNDDDQQQEEENNNDDD'

Shuffle except specific residues
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``'shuffle_except_specific_residues'`` variant type shuffles all residues except those specified.

**Parameters:**
- ``excluded_residues`` (list): List of residue types to exclude from shuffling

**Example:**

.. code-block:: python

    test = 'QQQEEENNNDDDQQQEEENNNDDD'
    variant_seq = create.variant(test, 'shuffle_except_specific_residues',
                                excluded_residues=['N', 'D'])
    print(variant_seq)
    # Output: 'QQQEEENNNDDDQQQEEENNNDDD'

Weighted shuffle specific residues
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``'weighted_shuffle_specific_residues'`` variant type performs weighted shuffling of specific residues.

**Parameters:**
- ``target_residues`` (list): List of residue types to shuffle
- ``shuffle_weight`` (float): Weight for shuffling operations (0.0 to 1.0)

**Example:**

.. code-block:: python

    test = 'QQQEEENNNDDDQQQEEENNNDDD'
    variant_seq = create.variant(test, 'weighted_shuffle_specific_residues',
                                target_residues=['Q', 'E'],
                                shuffle_weight=0.5)
    print(variant_seq)
    # Output: 'QQQEEENNNDDDQQQEEENNNDDD'

Targeted reposition specific residues
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``'targeted_reposition_specific_residues'`` variant type repositions specific residues within the sequence.

**Parameters:**
- ``target_residues`` (list): List of residue types to reposition

**Example:**

.. code-block:: python

    test = 'QQQEEENNNDDDQQQEEENNNDDD'
    variant_seq = create.variant(test, 'targeted_reposition_specific_residues',
                                target_residues=['E'])
    print(variant_seq)
    # Output: 'QQQEEENNNDDDQQQEEENNNDDD'

Property-based variants
-----------------------

Constant properties
~~~~~~~~~~~~~~~~~~~

The ``'constant_properties'`` variant type generates a variant where only the sequence properties are constrained.

**Parameters:**
- ``exclude_residues`` (list, optional): List of residue types to exclude from the variant

**Example:**

.. code-block:: python

    test = 'QEQNGVDQQETTPRQDYPGNQQPNQQAEGQQMQ'
    variant_seq = create.variant(test, 'constant_properties')
    print(variant_seq)
    # Output: 'QEQNGVDQQETTPRQDYPGNQQPNQQAEGQQMQ'

Constant residues and properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``'constant_residues_and_properties'`` variant type keeps specified residues constant while maintaining properties.

**Parameters:**
- ``constant_residues`` (list): List of residue types to keep constant

**Example:**

.. code-block:: python

    test = 'QEQNGVDQQETTPRQDYPGNQQPNQQAEGQQMQ'
    variant_seq = create.variant(test, 'constant_residues_and_properties',
                                constant_residues=['T', 'Q'])
    print(variant_seq)
    # Output: 'QDQSMNDQQETTGKQDNAGGQQHPQQPDAQQSQ'

Constant properties and class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``'constant_properties_and_class'`` variant type generates a variant with the same properties and amino acid class distribution.

**Example:**

.. code-block:: python

    test = 'QEQNGVDQQETTPRQDYPGNQQPNQQAEGQQMQ'
    variant_seq = create.variant(test, 'constant_properties_and_class')
    print(variant_seq)
    # Output: 'QENQGADQQDQNPRNEWPGNNNPNQTADGNSAT'

Constant properties and class by order
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``'constant_properties_and_class_by_order'`` variant type generates a variant with the same properties and maintains the order of amino acid classes.

**Example:**

.. code-block:: python

    test = 'QGENNENPQDQGSREGPQNNAWAQNNQDAQTSP'
    variant_seq = create.variant(test, 'constant_properties_and_class_by_order')
    print(variant_seq)
    # Output: 'QGDNQDNPNEQGQRDGPNTSAYAQQNNELQNNP'

Property modification variants
------------------------------

Change hydropathy constant class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``'change_hydropathy_constant_class'`` variant type changes hydropathy while keeping amino acid classes constant.

**Parameters:**
- ``target_hydropathy`` (float): Target hydropathy value

**Example:**

.. code-block:: python

    test = 'GNGGNRAENRTERKGEQTHKSNHNDGARHTDRRRSHDKNAASRE'
    variant_seq = create.variant(test, 'change_hydropathy_constant_class',
                                target_hydropathy=2.7)
    print(variant_seq)
    # Output: 'GTGGTKIETKTEKKGETTHKTTHTDGLKHTDRKKTHDKSVMTKE'

**Note:** Due to class constraints, there are limits to how much you can increase or decrease the hydropathy of any specific sequence. GOOSE will raise an error if you exceed these limits.

Change FCR minimize class changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``'change_fcr_minimize_class_changes'`` variant type adjusts FCR while minimizing changes to amino acid classes.

**Parameters:**
- ``target_FCR`` (float): Target FCR value

**Example:**

.. code-block:: python

    test = 'TTGGATSQAGGATHAQSHANSGTQSTSSPQTQGVNTTSANGQHGQATNQS'
    variant_seq = create.variant(test, 'change_fcr_minimize_class_changes',
                                target_FCR=0.2)
    print(variant_seq)
    # Output: 'TTGGMTSDAGGATHMKSHANSKGTKSTSSPKTEGINTTTIDGDHGKMTDKT'

Change NCPR constant class
~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``'change_ncpr_constant_class'`` variant type adjusts NCPR while keeping amino acid classes constant.

**Parameters:**
- ``target_NCPR`` (float): Target NCPR value

**Example:**

.. code-block:: python

    test = 'GNGGNRAENRTERKGEQTHKSNHNDGARHTDRRRSHDKNAASRE'
    variant_seq = create.variant(test, 'change_ncpr_constant_class',
                                target_NCPR=0.0)
    print(variant_seq)
    # Output: 'GNGGNRAENRTEEKGEQTHKSNHNDGARHTDDRRSHDKNAASRE'

Change kappa
~~~~~~~~~~~~

The ``'change_kappa'`` variant type alters charge asymmetry by changing the kappa value.

**Parameters:**
- ``target_kappa`` (float): Target kappa value (0.0 to 1.0)

**Example:**

.. code-block:: python

    test = 'QNEKRDQNEKRDQNEKRDQNEKRDQNEKRDQN'
    variant_seq = create.variant(test, 'change_kappa', target_kappa=0.9)
    print(variant_seq)
    # Output: 'KQRKRKRKRKRNQNQNQNQNEDEDQNEDEDED'

**Note:** GOOSE allows deviation from your input kappa value by up to 0.03 to maintain performance. Higher kappa values increase charge asymmetry, lower values reduce it.

Change any properties
~~~~~~~~~~~~~~~~~~~~~

The ``'change_any_properties'`` variant type adjusts multiple properties simultaneously.

**Parameters:**
- ``target_FCR`` (float): Target FCR value
- ``target_NCPR`` (float): Target NCPR value
- ``target_kappa`` (float): Target kappa value
- ``target_hydropathy`` (float): Target hydropathy value

**Example:**

.. code-block:: python

    test = 'GNGGNRAENRTERKGEQTHKSNHNDGARHTDRRRSHDKNAASRE'
    variant_seq = create.variant(test, 'change_any_properties',
                                target_hydropathy=2.5,
                                target_FCR=0.23,
                                target_NCPR=0.0,
                                target_kappa=0.1)
    print(variant_seq)
    # Output: 'GNGGQNAEQRNTKEGNESHTSTHTGDRAHQKSNNHQTNLERVSN'

Change properties minimize differences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``'change_properties_minimize_differences'`` variant type changes properties while minimizing differences from the original.

**Parameters (all optional):**
- ``target_hydropathy`` (float): Target hydropathy value
- ``target_FCR`` (float): Target FCR value
- ``target_NCPR`` (float): Target NCPR value
- ``target_kappa`` (float): Target kappa value

**Example:**

.. code-block:: python

    test = 'GNGGNRAENRTERKGEQTHKSNHNDGARHTDRRRSHDKNAASRE'
    variant_seq = create.variant(test, 'change_properties_minimize_differences',
                                target_kappa=0.3,
                                target_hydropathy=2.6)
    print(variant_seq)
    # Output: 'KTGGTKRGSKTARKGKSTHTTKHDEGVRTHDRRLSHEENADSTE'

Asymmetry variants
------------------

Change residue asymmetry
~~~~~~~~~~~~~~~~~~~~~~~~

The ``'change_residue_asymmetry'`` variant type changes the asymmetry of specific residues without changing sequence composition.

**Parameters:**
- ``target_residues`` (list): List of residue types or classes to modify
- ``num_changes`` (int, optional): Number of changes to make
- ``increase_or_decrease`` (str, optional): Whether to 'increase' or 'decrease' asymmetry

**Example - decreasing polar residue asymmetry:**

.. code-block:: python

    test = 'NSQSSQDSQDKSQGSQNQQEQSDSSEQTKQEEDGQTSSDSREQSQSHSQQ'
    variant_seq = create.variant(test, 'change_residue_asymmetry',
                                target_residues=['polar'],
                                increase_or_decrease='decrease',
                                num_changes=5)
    print(variant_seq)
    # Output: 'NSQDSSDQSQKSQGSQENQDQEKQSESSEQDGTQDQTSRSSEQSQSHSQQ'

**Example - increasing asymmetry with custom residue list:**

.. code-block:: python

    test = 'RGNNLAGIVLGAAGAMNGRTEGRKGEQTHGKSGNDDRGHTGDRSHGNKNRGE'
    variant_seq = create.variant(test, 'change_residue_asymmetry',
                                target_residues=['G', 'T'],
                                increase_or_decrease='increase',
                                num_changes=20)
    print(variant_seq)
    # Output: GGGGGTGGTGGGTGGGRNNLAIVLAAAMNRERKEQHKSNDDRHDRSHNKNRE

Dimensional variants
--------------------

Change dimensions
~~~~~~~~~~~~~~~~~

The ``'change_dimensions'`` variant type adjusts sequence dimensions (Rg or Re) while keeping amino acid composition constant.

**Parameters:**
- ``increase_or_decrease`` (str): Whether to 'increase' or 'decrease' the dimension
- ``rg_or_re`` (str): Whether to optimize 'rg' or 're'
- ``num_dim_attempts`` (int, optional): Number of dimensional optimization attempts
- ``allowed_error`` (float, optional): Allowed error for dimensional constraints
- ``reduce_pos_charged`` (bool, optional): Whether to reduce positive charges
- ``exclude_aas`` (list, optional): Amino acids to exclude from generation

**Example - increasing Re:**

.. code-block:: python

    test = 'FYFLGQGQQYYYYQQKQFFQFYYQQFFGFYGSNFQGGNYFGGYQQNQYFG'
    variant_seq = create.variant(test, 'change_dimensions',
                                increase_or_decrease='increase',
                                rg_or_re='re')
    print(variant_seq)

**Example - decreasing Rg:**

.. code-block:: python

    test = 'FYFLGQGQQYYYYQQKQFFQFYYQQFFGFYGSNFQGGNYFGGYQQNQYFG'
    variant_seq = create.variant(test, 'change_dimensions',
                                increase_or_decrease='decrease',
                                rg_or_re='rg')
    print(variant_seq)

Error handling and troubleshooting
-----------------------------------

The ``variant()`` function provides comprehensive error handling:

**Common errors:**

1. **Invalid variant type:** Ensure the variant_type is one of the supported types listed above.
2. **Missing required parameters:** Each variant type has specific required parameters.
3. **Invalid parameter values:** Check that parameter values are within valid ranges.
4. **Variant generation failure:** If generation fails, try increasing ``num_attempts`` or adjusting target values.

**Example error handling:**

.. code-block:: python

    try:
        variant_seq = create.variant(sequence, 'change_kappa', target_kappa=0.5)
    except goose.goose_exceptions.GooseInputError as e:
        print(f"Input error: {e}")
    except goose.goose_exceptions.GooseFail as e:
        print(f"Generation failed: {e}")

**Tips for successful variant generation:**

- Start with moderate changes to properties
- Use higher ``num_attempts`` for difficult targets
- Check that your sequence has the necessary residue types for the variant
- For kappa variants, ensure your sequence has both positive and negative charges
- For class-based variants, remember that some property changes may not be possible due to class constraints

Function selection guide
------------------------

**Choose variant type based on your needs:**

- **Shuffling sequences:** Use shuffling variants to rearrange existing residues
- **Maintaining properties:** Use constant property variants to keep sequence characteristics
- **Changing specific properties:** Use property modification variants for targeted changes
- **Adjusting dimensions:** Use dimensional variants to change IDR dimensions 
- **Changing asymmetry:** Use asymmetry variants to modify residue distribution patterns

**Performance considerations:**

- Shuffling variants are generally fastest
- Property modification variants may require more attempts
- Dimensional variants can be computationally intensive
- Kappa variants work best with values between 0.1 and 0.9

Backward compatibility notes
----------------------------

The unified ``variant()`` function replaces many individual functions from previous versions:

- ``constant_class_var()`` → ``variant(seq, 'constant_properties_and_class')``
- ``constant_properties_var()`` → ``variant(seq, 'constant_properties')``
- ``region_shuffle_var()`` → ``variant(seq, 'shuffle_specific_regions')``
- ``targeted_shuffle_var()`` → ``variant(seq, 'shuffle_specific_residues')``
- ``excluded_shuffle_var()`` → ``variant(seq, 'shuffle_except_specific_residues')``
- ``kappa_var()`` → ``variant(seq, 'change_kappa')``
- ``hydro_class_var()`` → ``variant(seq, 'change_hydropathy_constant_class')``
- ``fcr_class_var()`` → ``variant(seq, 'change_fcr_minimize_class_changes')``
- ``ncpr_class_var()`` → ``variant(seq, 'change_ncpr_constant_class')``
- ``all_props_class_var()`` → ``variant(seq, 'change_any_properties')``
- ``re_var()`` / ``rg_var()`` → ``variant(seq, 'change_dimensions')``
- ``weighted_shuffle_var()`` → ``variant(seq, 'weighted_shuffle_specific_residues')``
- ``asymmetry_var()`` → ``variant(seq, 'change_residue_asymmetry')``

The new interface provides more consistent parameter names and improved error handling while maintaining all the functionality of the original functions.
