Sequence Optimization with GOOSE
================================

**SUPER DUPER IMPORTANT NOTE**: GOOSE is an IDR design tool. **HOWEVER**, when using SequenceOptimizer, you can design anything you want. Thus, sequences are not guaranteed to be disordered unless you specify the ``FractionDisorder`` property. You can also use the ``MatchSequenceDisorder`` property if you'd like to make a variant that maintains the disorder profile of your original sequence (either exactly or as a minimum threshold).

GOOSE's ``SequenceOptimizer`` is a flexible tool for designing protein sequences that match user-defined biophysical or biochemical property targets. It uses stochastic optimization, mutating an initial sequence to minimize the difference between calculated and target property values. Importantly, you can simultaneously specify arbitrary numbers of properties to optimize towards. Further, if GOOSE lacks properties that you want to optimize towards, you can easily implement it by following our implementation examples below.

Critical Differences between SequenceOptimizer and Create Functionality
-----------------------------------------------------------------------

We consider ``SequenceOptimizer`` to be a distinct approach to sequence / sequence variant generation from our original functionality as in the ``create`` module. The key differences are as follows:

* ``SequenceOptimizer`` is not fully optimized for every property it *can* optimize towards. What we mean by this is that sequences are basically generated by minimizing the error in the sequence through stochastic optimization. While most GOOSE functionality implements optimization like this to some extent, the code in the ``create`` module is highly optimized. The cost of this is that implementing new features in ``create`` is painful. In contrast, implementing new functionality in ``SequenceOptimizer`` is absolutely trivial. Furthermore, ``SequenceOptimizer`` is far, far, far more flexible.
* ``SequenceOptimizer`` will return you a sequence that is not exactly what you want it to be. It will return the *most optimized sequence*. In contrast, ``create`` will raise an Exception if it is unable to generate your desired sequence. Therefore, we recommend double checking that the sequence you make using ``SequenceOptimizer`` is what you expect it to be (we recommend this with ``create`` as well, just to be safe).

.. contents:: Table of Contents
   :local:
   :depth: 2

Quick Start Example
-------------------

Design a sequence of length 50 with a target hydrophobicity:

.. code-block:: python

    import goose
    from sparrow import Protein

    optimizer = goose.SequenceOptimizer(target_length=50)
    optimizer.add_property(goose.Hydrophobicity, target_value=0.5, weight=1.0)
    optimizer.set_optimization_params(max_iterations=1000, tolerance=0.01)
    print("Starting optimization...")
    optimized_sequence = optimizer.run()
    print(f"Optimized Sequence: {optimized_sequence}")
    final_protein = Protein(optimized_sequence)
    print(f"Final Hydrophobicity: {final_protein.hydrophobicity:.2f}")

**Explanation:**
- ``SequenceOptimizer(target_length=50)``: Sets the sequence length.
- ``add_property(...)``: Adds a property to optimize (here, hydrophobicity).
- ``set_optimization_params(...)``: Sets optimization parameters (optional).
- ``run()``: Runs the optimization and returns the best sequence.

Property Classes Overview
-------------------------

GOOSE provides many built-in property classes. Each can be added to the optimizer with its required arguments. Here is a summary:

+-------------------------------+-----------------------------------------------+------------------------------------------------+
| Property Class                | Description                                   | Key Arguments                                  |
+===============================+===============================================+================================================+
| Hydrophobicity                | Target average hydrophobicity                 | target_value, weight                           |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| FCR                           | Fraction of Charged Residues                  | target_value, weight                           |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| NCPR                          | Net Charge Per Residue                        | target_value, weight                           |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| Kappa                         | Charge patterning parameter                   | target_value, weight                           |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| SCD                           | Sequence Charge Decoration                    | target_value, weight                           |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| ComputeIWD                    | Inverse Weighted Distance for residues        | residues (tuple), target_value, weight         |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| RadiusOfGyration              | Predicted radius of gyration                  | target_value, weight                           |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| EndToEndDistance              | Predicted end-to-end distance                 | target_value, weight                           |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| FractionDisorder              | Fraction of disordered residues               | target_value, weight, disorder_cutoff          |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| TargetAminoAcidFractions      | Target amino acid composition                 | target_fractions (dict), weight                |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| MaxFractions, MinFractions    | Max/min fractions for specific amino acids    | max_fractions/min_fractions (dict), weight     |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| SelfEpsilon                   | Self-interaction potential                    | target_value, weight, model                    |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| MatchSequenceDisorder         | Match disorder profile of a target sequence   | target_sequence, weight, exact_match,          |
|                               |                                               | target_value                                   |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| Complexity                    | Wootton-Federhen (SEG) sequence complexity    | target_value, weight                           |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| MatchingResidues              | Number of residues matching a target sequence | target_sequence, target_value, weight          |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| MaxMatchingResidues           | Penalize if matches exceed target value       | target_sequence, target_value, weight          |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| EpsilonVectorBySequence       | Match epsilon interaction vectors to target   | original_sequence, target_interacting_sequence,|
|                               |                                               | weight, model                                  |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| EpsilonByValue                | Target total epsilon interaction value        | target_value, target_sequence, weight, model   |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| EpsilonBySequence             | Match total epsilon value to reference        | original_sequence, target_interacting_sequence,|
|                               |                                               | target_value, weight, model                    |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| FDSurfaceInteractionByValue   | Target surface repulsion/attraction (folded)  | repulsive_target, attractive_target, weight,   |
|                               |                                               | model, path_to_pdb, probe_radius,              |
|                               |                                               | surface_thresh, sasa_mode, fd_start, fd_end,   |
|                               |                                               | preloaded_fd                                   |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| ChemicalFingerprint           | Match chemical fingerprint to target sequence | target_sequence, target_value, weight, model   |
+-------------------------------+-----------------------------------------------+------------------------------------------------+

.. note::
   Some properties (e.g., ComputeIWD, TargetAminoAcidFractions) require extra arguments. See the API or docstrings for details.

Optimizer Initialization and Parameters
---------------------------------------

You can control the optimizer's behavior with several parameters:

.. code-block:: python

    optimizer = goose.SequenceOptimizer(
        target_length=100,
        kmer_dict_file=None,  # Path to a custom k-mer bias pickle file, or None for default amino acid frequencies
        verbose=True,         # Print logging information during optimization
        gap_to_report=100,    # How often to update progress (e.g., every 100 iterations)
        num_shuffles=5,       # Number of global/local shuffles to try at shuffle intervals
        just_shuffle=False    # If True, only shuffles the sequence without k-mer mutations
    )

- **kmer_dict_file**: Use a custom k-mer dictionary (see below for details).
- **just_shuffle**: If True, only shuffles the sequence (preserves composition).

Set optimization parameters at any time:

.. code-block:: python

    optimizer.set_optimization_params(
        max_iterations=50000,  # Increasing iterations can allow you to optimize towards very hard to make sequences.
        tolerance=1e-3,        # Setting tolerance very low can be useful but not necessary depending on the parameter.
        window_size=15,        # Window size for local shuffling
        shuffle_interval=50,   # Perform shuffles every 50 iterations
        just_shuffle=False     # Can also be set here
    )

- **window_size**: Size of sequence segments for local shuffling.
- **shuffle_interval**: How often to perform shuffles during optimization.

Multiple Properties and Weights
-------------------------------

You can optimize for several properties at once. The ``weight`` argument controls the importance of each property in the combined error function.

.. code-block:: python

    import goose
    from sparrow import Protein

    optimizer = goose.SequenceOptimizer(target_length=75, verbose=False)
    optimizer.add_property(goose.FCR, target_value=0.39, weight=1.0)
    optimizer.add_property(goose.NCPR, target_value=-0.1, weight=1.5) # NCPR is more important
    optimizer.set_optimization_params(max_iterations=10000)
    optimized_sequence = optimizer.run()

    print(f"Optimized Sequence: {optimized_sequence}")
    final_protein = Protein(optimized_sequence)
    print(f"Final FCR: {final_protein.FCR:.2f}")
    print(f"Final NCPR: {final_protein.NCPR:.2f}")

.. note::
   The optimizer minimizes a weighted sum of property errors: sum(weight * |calculated - target|).

Advanced Features
-----------------

**Using an Initial Sequence**

Start from a specific sequence (must match ``target_length``):

.. code-block:: python

    initial_seq = "M" * optimizer.target_length
    optimizer.set_initial_sequence(initial_seq)

**Setting Fixed Ranges**

Keep certain regions unchanged during optimization (0-indexed, inclusive):

.. code-block:: python

    optimizer.set_fixed_ranges([(0, 9), (20, 29)])

This preserves residues 0-9 and 20-29.

Custom Properties
-----------------

If you need a property not provided by GOOSE, define your own by subclassing ``goose.backend.optimizer_properties.ProteinProperty``:

.. code-block:: python

    import goose
    from goose.backend.optimizer_properties import ProteinProperty
    import sparrow

    class AlanineCount(ProteinProperty):
        def __init__(self, target_value: float, weight: float = 1.0):
            super().__init__(target_value, weight)
        def calculate(self, protein: 'sparrow.Protein') -> float:
            return float(protein.sequence.count('A'))

    custom_optimizer = goose.SequenceOptimizer(target_length=30)
    custom_optimizer.add_property(AlanineCount, target_value=5.0, weight=1.0)
    custom_optimizer.set_optimization_params(max_iterations=500)
    custom_sequence = custom_optimizer.run()
    print(f"Custom Optimized Sequence: {custom_sequence}")
    print(f"Alanine count: {custom_sequence.count('A')}")

K-mer Dictionaries
------------------

A k-mer dictionary controls the amino acid or k-mer composition during sequence generation. By default, GOOSE uses single amino acid frequencies from ``amino_acids.py``. You can provide a custom dictionary via the ``kmer_dict_file`` argument (must be a pickle file with the correct format).

- Use a custom k-mer dictionary to bias sequence generation toward specific motifs or patterns.
- See the API for details on the expected format.

How the Optimizer Works
-----------------------

1. **Initialization**: Builds a starting sequence (random or user-provided).
2. **Mutation**: At each iteration, mutates the sequence (by k-mer replacement or shuffling).
3. **Property Calculation**: Calculates all property values for the new sequence.
4. **Error Calculation**: Computes the weighted sum of errors between calculated and target values.
5. **Selection**: Keeps the best sequence found so far.
6. **Stopping**: Stops when the error is below ``tolerance`` or ``max_iterations`` is reached.

Troubleshooting and Tips
------------------------

**Optimization not converging?**
- Increase ``max_iterations``.
- Check if your property targets are physically possible.
- Increase ``num_shuffles`` to escape local minima.

**Slow optimization?**
- Decrease ``max_iterations`` for faster (but less optimal) results.
- Set a reasonable ``tolerance`` (too small = slow).
- Increase ``gap_to_report`` to reduce logging overhead.

**Multiple property conflicts?**
- Adjust ``weight`` parameters to prioritize properties.
- Ensure your targets are compatible (e.g., high hydrophobicity and high charge may conflict).

**Fixed range issues?**
- Don't over-constrain the sequence with too many fixed regions.
- Remember: fixed ranges are 0-indexed and inclusive.

Glossary
--------

- **k-mer**: A substring of length k (e.g., 3-mer = 3 amino acids).
- **window_size**: The length of sequence segments for local shuffling.
- **shuffle_interval**: How often shuffling is performed during optimization.
- **fixed ranges**: Sequence regions that are not mutated.
- **property**: A biophysical or biochemical feature to optimize (e.g., hydrophobicity).
- **weight**: The importance of a property in the optimization objective.

See Also
--------

- :doc:`api`
- :doc:`getting_started`
- :doc:`sequence_generation`
- :doc:`variant_generation`
- :doc:`sequence_library_generation`
- :doc:`sequence_analysis`

For more details, see the API documentation or the source code in ``goose/optimize.py`` and ``goose/backend/optimizer_properties.py``.
