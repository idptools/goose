Creating sequences and sequence variants with GOOSE
=====================================================

GOOSE provides several functions for generating intrinsically disordered protein sequences with specific properties. This guide covers the main sequence generation functions:

- ``create.sequence()``: Generate sequences with specified physicochemical properties
- ``create.seq_by_fractions()``: Generate sequences with specified amino acid fractions  
- ``create.seq_by_classes()``: Generate sequences with specified amino acid class fractions
- ``create.seq_by_rg()``: Generate sequences with specified radius of gyration
- ``create.seq_by_re()``: Generate sequences with specified end-to-end distance

Getting Started
===============

First import ``create`` from goose.

.. code-block:: python

    from goose import create

Once ``create`` has been imported, you can start making sequences.

**Available Organism-Specific Probabilities**

GOOSE provides amino acid probabilities derived from intrinsically disordered regions (IDRs) of various organisms. These can be used in the ``custom_probabilities`` parameter for ``create.sequence()`` and ``remaining_probabilities`` parameter for ``create.seq_by_fractions()`` and ``create.seq_by_classes()``:

- **Organism-specific**: ``'mouse'``, ``'human'``, ``'fly'``, ``'neurospora'``, ``'yeast'``, ``'arabidopsis'``, ``'e_coli'``, ``'worm'``, ``'zebrafish'``, ``'frog'``, ``'dictyostelium'``
- **Special sets**: ``'unbiased'`` (equal probabilities for all amino acids), ``'all'`` (weighted average across all organisms)

These probabilities are based on MetaPredict V3 predictions of IDRs across different proteomes and can help generate sequences that are more representative of natural disordered regions in specific organisms.
Specific organisms are:
- **Mouse**: Mus musculus
- **Human**: Homo sapiens
- **Fly**: Drosophila melanogaster
- **Neurospora**: Neurospora crassa
- **Yeast**: Saccharomyces cerevisiae
- **Arabidopsis**: Arabidopsis thaliana
- **E. coli**: Escherichia coli
- **Worm**: Caenorhabditis elegans
- **Zebrafish**: Danio rerio
- **Frog**: Xenopus laevis
- **Dictyostelium**: Dictyostelium discoideum

Generating sequences with specified properties
===============================================

The ``create.sequence()`` function lets you create sequences predicted to be disordered with various specified properties. 

The only required argument is the length, which must be between the minimum and maximum allowed lengths as defined in the parameters module. In addition, you can also specify several parameters.

**Primary physicochemical properties:**

1. ``hydropathy``: Average hydropathy, which must be between 0 and 6.1. It is called `hydropathy` instead of `hydrophobicity` because we are using a rescaled Kyte-Doolittle hydropathy scale. Higher values indicate more hydrophobic sequences.

   .. note::
      The legacy parameter name ``hydrophobicity`` is still supported for backward compatibility but ``hydropathy`` is recommended.

2. ``FCR``: The fraction of charged residues (FCR), which must be between 0 and 1. This includes both positively and negatively charged residues.

3. ``NCPR``: The net charge per residue (NCPR), which must be between -1 and +1. Positive values indicate net positive charge, negative values indicate net negative charge.

4. ``kappa``: The kappa value is a charge patterning parameter where higher values mean more even charge distribution. This must be between 0 and 1. NOTE: kappa values of 0 or 1 can take a long time to make because there are very few sequences that satisfy these values. GOOSE is much faster at making sequences between 0.1 and 0.9. Also, specifying kappa **requires** that you do not have an FCR=0 and that you do not have FCR=NCPR (there needs to be some + and some - charged residues for kappa).

**Importantly, you do not need to specify all of these sequence properties simultaneously.** For example, if you specify FCR and hydropathy, GOOSE will return sequences that have varying NCPR and kappa values while making sure that the specified hydropathy and FCR values are what you input. In this way, you can generate many sequences that have the fixed properties that you want to stay fixed, while other properties vary.

**Additional parameters:**

1. ``disorder_cutoff``: The disorder threshold for sequence validation. Sequences must have disorder scores above this threshold. Default value comes from the parameters module. We use MetaPredict V3 in GOOSE (see https://github.com/idptools/metapredict).

   .. note::
      The legacy parameter name ``cutoff`` is still supported for backward compatibility but ``disorder_cutoff`` is recommended.

2. ``attempts``: The number of attempts defines how many times GOOSE will try to generate a desired sequence. Default is 20. Higher values increase success probability but take longer. For more complex sequence compositions you may need to increase this from the default value.

3. ``exclude``: A list of residues to exclude from your generated sequence. **NOTE** - If you specify FCR or NCPR, you CANNOT have any charged residues in the ``exclude`` list. Additionally, **you can not specify more than 10 residues to be excluded**.

4. ``use_weighted_probabilities``: Whether to use weighted probabilities for generating sequences. Default is False. The TL;DR is that we have generated many thousands of sequences that fit different parameter combinations and created weighted probabilities for each residue to increase the chance of generating the objective sequence. This was necessary in older versions of GOOSE. However, we have implemented Numpy vectorization to speed up sequence generation, and the weighted probabilities are not necessary anymore (though can still speed things up). However, if you want to use them, you can set this to True.

5. ``strict_disorder``: Whether you want **all** residues to be above the cutoff value. By default, GOOSE lets a small number (at most 5 for very large sequences) of residues be below the disorder threshold because a single amino acid (or even a few) below the threshold is not realistically going to be a folded domain.

6. ``return_all_sequences``: Whether to return all sequences generated. Default is False. This is new as of v0.2.0 because GOOSE generates and checks many sequences simultaneously, and it is now possible to return all of them that match the specified parameters.

7. ``custom_probabilities``: You can now specify your own custom probabilities for generating sequences. This can be either:
   
   - A dictionary where the keys are amino acids (single-letter codes) and the values are probabilities (sum must equal 1)
   - A string specifying organism-specific probabilities or predefined probability sets:
     
     - Organism-specific: ``'mouse'``, ``'human'``, ``'fly'``, ``'neurospora'``, ``'yeast'``, ``'arabidopsis'``, ``'e_coli'``, ``'worm'``, ``'zebrafish'``, ``'frog'``, ``'dictyostelium'``
     - ``'unbiased'``: Equal probabilities (5%) for all 20 amino acids
     - ``'all'``: Weighted average probabilities across all organisms

8. ``metapredict_version``: You can specify the version of MetaPredict you want to use. The default is 3, but you can set it to 2 or 1 if you want to use the older versions.

9. ``max_consecutive_ordered``: Maximum number of consecutive ordered residues allowed in the sequence. Default value comes from the parameters module.

10. ``max_total_ordered``: Maximum fraction of ordered residues allowed in the sequence. Default value comes from the parameters module.

11. ``batch_size``: Number of sequences to generate in each batch. Default value comes from the parameters module.

12. ``hydropathy_tolerance``: Tolerance for hydropathy matching. Default value comes from the parameters module.

13. ``kappa_tolerance``: Tolerance for kappa matching. Default value comes from the parameters module.


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

    # Using a custom probability dictionary
    custom_probs = {'A':0.1, 'R':0.1, 'D':0.1, 'E':0.1, 'G':0.1, 'H':0.1, 'I':0.1, 'K':0.1, 'L':0.1, 'M':0.1}
    create.sequence(100, hydropathy=2.5, custom_probabilities=custom_probs)
    
    # Using organism-specific probabilities
    create.sequence(100, hydropathy=2.5, custom_probabilities='human')
    
    # Using unbiased probabilities (equal for all amino acids)
    create.sequence(100, hydropathy=2.5, custom_probabilities='unbiased')
    
    # Using averaged probabilities across all organisms
    create.sequence(100, hydropathy=2.5, custom_probabilities='all')

**Specifying metapredict version:**

.. code-block:: python

    create.sequence(100, hydropathy=2.5, metapredict_version=2)
    'RDRFSEYKNTKEQAFDSYQLERHKERESQTRKRHRPQREKQRPDGERHKHEFMEWKLERRRCTEDGDKEFRLQALGRCESPIGMQMHTPDIADPKRDRRN'


**Specifying additional parameters:**

.. code-block:: python

    # Using custom disorder cutoff and batch size
    create.sequence(100, hydropathy=2.5, disorder_cutoff=0.6, batch_size=50)
    
    # Using tolerance parameters for better matching
    create.sequence(100, hydropathy=3.0, kappa=0.5, hydropathy_tolerance=0.1, kappa_tolerance=0.05)
    
    # Limiting ordered residues
    create.sequence(100, FCR=0.3, max_consecutive_ordered=3, max_total_ordered=0.1)


Error Handling
==============

GOOSE provides informative error messages when sequence generation fails or invalid parameters are provided:

**GooseInputError**: Raised when invalid parameters are provided, such as:
- Invalid sequence length
- Parameter values outside allowed ranges  
- Invalid parameter combinations
- Missing required parameters

**GooseFail**: Raised when sequence generation fails after all attempts, such as:
- Unable to generate sequence with specified properties
- Conflicting parameter constraints
- Insufficient attempts for complex sequences

**Tips for successful sequence generation:**

1. **Increase attempts**: For complex parameter combinations, increase the ``attempts`` parameter
2. **Adjust tolerances**: Use ``hydropathy_tolerance`` and ``kappa_tolerance`` for more flexible matching
3. **Check parameter ranges**: Ensure all parameters are within valid ranges
4. **Use batch generation**: Set ``return_all_sequences=True`` to get multiple sequences
5. **Optimize disorder settings**: Adjust ``disorder_cutoff`` and ``strict_disorder`` if needed


Generating Sequences specifying Fractions of Amino Acids
=========================================================

The ``create.seq_by_fractions()`` function lets you create sequences predicted to be disordered with specified fractions of various amino acids. This function provides fine-grained control over sequence composition by allowing you to specify the exact fraction of each amino acid type. With this function, you can specify multiple amino acids simultaneously. Each fraction should be specified using a decimal value (for example, if you want one-tenth of the amino acids to be alanine use ``A=0.1``).

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

Note that if you pass in requested fractions, those fractions cannot be greater than 1, and the sum of all specified fractions should not exceed 1. Any values that are remaining will be randomly added based on the remaining probabilities. 

In addition to specifying the specific amino acid fractions, other parameters can be passed to the ``create.seq_by_fractions()`` function:

1. ``disorder_cutoff``: The disorder threshold for sequence validation. Default is 0.6.

2. ``attempts``: The number of attempts defines how many times GOOSE will try to generate a desired sequence. Default is 100.

3. ``max_aa_fractions``: If you wish to generate sequences with extreme compositions it may be necessary to over-ride the default max fractional values. This can be achieved by passing a max_aa_fractions dictionary, which should specify key-value pairs for amino acid-max fraction information.

4. ``strict_disorder``: Whether you want **all** residues to be above the cutoff value. By default, GOOSE lets a small number (at most 5 for very large sequences) of residues be below the disorder threshold because a single amino acid (or even a few) below the threshold is not realistically going to be a folded domain.

5. ``remaining_probabilities``: Custom probabilities for amino acids not explicitly specified in fractions. Keys should be amino acid codes, values should be probabilities.

6. ``return_all_sequences``: Whether to return all sequences generated. Default is False. This is new as of v0.2.0 because GOOSE generates and checks many sequences simultaneously, and it is now possible to return all of them that match the specified parameters.

7. ``metapredict_version``: You can specify the version of MetaPredict you want to use. The default is 3, but you can set it to 2 or 1 if you want to use the older versions.

8. ``max_consecutive_ordered``: Maximum number of consecutive ordered residues allowed in the sequence.

9. ``max_total_ordered``: Maximum fraction of ordered residues allowed in the sequence.

10. ``batch_size``: Number of sequences to generate in each batch.


Examples of Sequence Generation by Fractions
---------------------------------------------

**Specifying a single amino acid fraction:**

.. code-block:: python

    create.seq_by_fractions(100, Q=0.3)
    'QEQNGVDQQETTPRQDYPGNQQPNQQAEGQQMQSTKMHDQHDSVNEDQEQNQNPWGHQPHMKGESNSSAREAQSEDQQNQAQNQQQNHDSTQQQDGQMDQ'

**Specifying multiple amino acids:**

.. code-block:: python

    create.seq_by_fractions(100, Q=0.3, S=0.3, E=0.1)
    'QEQQSQKASQSQVESQDSSESSAPGSSQMHQQQSQSQEGMEQHQSSVGNSSSYPQSEQSEQQRQQSSQDQQQQSSSQTSEENSQSRQHDMSDTEMSGSQR'

**Using organism-specific probabilities:**

.. code-block:: python

    # Use human-specific probabilities for remaining amino acids
    create.seq_by_fractions(100, Q=0.3, S=0.3, remaining_probabilities='human')
    
    # Use mouse-specific probabilities for remaining amino acids
    create.seq_by_fractions(100, Q=0.3, S=0.3, remaining_probabilities='mouse')
    
    # Use unbiased probabilities for remaining amino acids
    create.seq_by_fractions(100, Q=0.3, S=0.3, remaining_probabilities='unbiased')

**Note** - 
Some combinations of amino acids are simply not possible to make that are predicted to be disordered using the default settings. Specifically, specifying high fractions of multiple aliphatics or aromatics may not be predicted to be disordered using the default cutoff value.

**Excluding a specific amino acids:**
If you want to exclude an amino acid, you can set it equal to 0.

.. code-block:: python

    create.seq_by_fractions(50, A=0)
    'NKERPTGSWDEPPFDEGSSGMTNEDMGNKPYPTTDMQPEKWPQNDQQGST'
    

**Overriding default max fractions:**  

.. code-block:: python

    create.seq_by_fractions(100, Y=0.5, max_aa_fractions={'Y':1}) 
    'SSYYYYYSYSSYYSYSSGHYYSYSSYYYSSSYYSSYGGTYGYYSYSYGYYSSYYYSYSSNYYYYYYYYSSYGNSGYGGYYSYYSSSQHHYSSYYYSYYSY'


Generating Sequences specifying Amino Acid Classes
==================================================

The ``create.seq_by_classes()`` function lets you create sequences with specified fractions of amino acid classes rather than individual amino acids. This provides a higher-level approach to sequence composition control.

**Amino acid classes:**

- **Aromatic**: F, W, Y (phenylalanine, tryptophan, tyrosine)
- **Aliphatic**: A, I, L, V (alanine, isoleucine, leucine, valine)  
- **Polar**: N, Q, S, T (asparagine, glutamine, serine, threonine)
- **Positive**: K, R (lysine, arginine)
- **Negative**: D, E (aspartate, glutamate)
- **Glycine**: G (glycine)
- **Proline**: P (proline)
- **Cysteine**: C (cysteine)
- **Histidine**: H (histidine)

**Parameters:**

All class fractions should be between 0 and 1. Additional parameters include:

- ``num_attempts``: Number of attempts to generate the sequence (default: 10)
- ``strict_disorder``: Whether to use strict disorder checking (default: False)
- ``disorder_cutoff``: Disorder threshold for sequence validation (default: from parameters)
- ``metapredict_version``: Version of MetaPredict to use (default: 3)
- ``max_consecutive_ordered``: Maximum consecutive ordered residues allowed
- ``max_total_ordered``: Maximum fraction of ordered residues allowed
- ``remaining_probabilities``: Custom amino acid probabilities for sequence generation. This controls the probabilities for amino acids not covered by class specifications. This can be either:
  
  - A dictionary where the keys are amino acids (single-letter codes) and the values are probabilities (sum must equal 1)
  - A string specifying organism-specific probabilities or predefined probability sets:
    
    - Organism-specific: ``'mouse'``, ``'human'``, ``'fly'``, ``'neurospora'``, ``'yeast'``, ``'arabidopsis'``, ``'e_coli'``, ``'worm'``, ``'zebrafish'``, ``'frog'``, ``'dictyostelium'``
    - ``'unbiased'``: Equal probabilities (5%) for all 20 amino acids
    - ``'all'``: Weighted average probabilities across all organisms

**Examples:**

.. code-block:: python

    # Generate sequence with 20% aromatic and 10% positive residues
    create.seq_by_classes(100, aromatic=0.2, positive=0.1)
    
    # Generate sequence with multiple class constraints
    create.seq_by_classes(75, aromatic=0.15, polar=0.25, glycine=0.1)
    
    # Use organism-specific probabilities for remaining amino acids
    create.seq_by_classes(100, aromatic=0.2, positive=0.1, remaining_probabilities='human')
    
    # Use unbiased probabilities for remaining amino acids
    create.seq_by_classes(100, aromatic=0.2, positive=0.1, remaining_probabilities='unbiased')



Generating Sequences specifying Ensemble Dimensions
=========================================================

The ``create.seq_by_rg()`` and ``create.seq_by_re()`` functions let you create sequences with a specified length and a predicted radius of gyration (Rg) or end-to-end distance (Re). For these functions, you must specify the length and an objective Re or Rg. In addition you can also specify:

1. ``disorder_cutoff``: The disorder threshold for sequence validation. Default value comes from the parameters module.

2. ``attempts``: The number of attempts defines how many times GOOSE will try to generate a desired sequence. Default is 20. 

3. ``strict_disorder``: Whether you want **all** residues to be above the cutoff value. By default, GOOSE lets a small number (at most 5 for very large sequences) of residues be below the disorder threshold because a single amino acid (or even a few) below the threshold is not realistically going to be a folded domain.

4. ``exclude_aas``: A list of residues to exclude from your generated sequence. There are some limitations on excluding AAs, specifically you can't simultaneously exclude  W, Y, G, F, Q, and N or D, E, K, P, S, and T.

5. ``allowed_error``: How far off from your desired Re/Rg in Å GOOSE can be before returning the sequence. A higher value here will decrease the time it takes GOOSE to make the sequence. Default value comes from the parameters module.

6. ``reduce_pos_charged``: Whether to reduce positively charged amino acids in the sequence. Default is False. The reason for this is that in vivo data suggests that positively charged residues may not drive sequence expansion as much as was predicted by the model used here for predicted rg / re. Therefore, when set to True, this function will largely avoid high numbers of (+) charged residues if possible.

7. ``metapredict_version``: Version of MetaPredict to use for disorder prediction. Default is 3.

8. ``max_consecutive_ordered``: Maximum number of consecutive ordered residues allowed in the sequence. Default value comes from the parameters module.

9. ``max_total_ordered``: Maximum fraction of ordered residues allowed in the sequence. Default value comes from the parameters module.

Examples of generating sequences by specifying Rg or Re
----------------------------------------------------------

**Specifying a length and Rg:**

.. code-block:: python

    create.seq_by_rg(50, 20)
    'NSETSEFYNDPVNAQPGDDHNSENNSVTYDNTGTYSNEFPDTEPSDLHAP'


**Specifying a length and Re:**

.. code-block:: python

    create.seq_by_re(50, 20)
    'FGQQGGQWGQWGNGQWGYWQNFGYGGNGGWYFYQWYNWFQYNWWFWQWWF'

  
**Specifying a length and Rg, allowing positive charged residues:**

.. code-block:: python

    create.seq_by_rg(50, 20, reduce_pos_charged=True)
    'NQKDSPEIDKPKPGNASGKFQTIRGNNRRKQKGGQGYPEKTIGERHMSEA'


**Specifying a length and Re with custom error tolerance:**

.. code-block:: python

    create.seq_by_re(75, 40.0, allowed_error=2.0)
    'NQKDSPEIDKPKPGNASGKFQTIRGNNRRKQKGGQGYPEKTIGERHMSEA'


**Specifying a length and Rg with excluded amino acids:**

.. code-block:: python

    create.seq_by_rg(100, 25.0, exclude_aas=['C', 'M'])
    'NQKDSPEIDKPKPGNASGKFQTIRGNNRRKQKGGQGYPEKTIGERHMSEA'


Function Selection Guide
========================

Choose the appropriate function based on your needs:

**create.sequence()**: 
- **Best for**: Specifying physicochemical properties (charge, hydropathy, etc.)
- **Use when**: You want to control FCR, NCPR, hydropathy, or kappa values
- **Flexibility**: High - can combine multiple properties

**create.seq_by_fractions()**: 
- **Best for**: Precise amino acid composition control
- **Use when**: You need exact percentages of specific amino acids
- **Flexibility**: High - can specify any combination of amino acids

**create.seq_by_classes()**: 
- **Best for**: Controlling amino acid classes (aromatic, charged, etc.)
- **Use when**: You want broad compositional control without specifying individual amino acids
- **Flexibility**: Medium - works with predefined amino acid groups

**create.seq_by_rg()**: 
- **Best for**: Controlling sequence compactness
- **Use when**: You need a specific radius of gyration value
- **Flexibility**: Low - focused on dimensional properties

**create.seq_by_re()**: 
- **Best for**: Controlling end-to-end distance
- **Use when**: You need a specific end-to-end distance value  
- **Flexibility**: Low - focused on dimensional properties

**Parameter Compatibility**:
- All functions support disorder-related parameters (``disorder_cutoff``, ``strict_disorder``, ``metapredict_version``)
- Most functions support sequence generation parameters (``attempts``, ``return_all_sequences``, ``batch_size``)
- Dimensional functions (``seq_by_rg``, ``seq_by_re``) have specialized parameters for size control


Copyright (c) 2023, Ryan Emenecker - Holehouse Lab