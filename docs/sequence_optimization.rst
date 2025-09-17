Using the SequenceOptimizer
==============================

GOOSE's ``SequenceOptimizer`` is a flexible tool for designing protein sequences that match user-defined biophysical or biochemical property targets. It uses stochastic optimization with adaptive scaling to  explore sequence space and minimize the difference between calculated and target property values. You can simultaneously optimize towards arbitrary numbers of properties with individual weights, tolerances, and constraint types.

**SUPER DUPER IMPORTANT NOTE**: GOOSE is an IDR design tool. **HOWEVER**, when using SequenceOptimizer, you can design anything you want. Thus, sequences are not guaranteed to be disordered unless you specify the ``FractionDisorder`` property. You can also use the ``MatchSequenceDisorder`` property if you'd like to make a variant that maintains the disorder profile of your original sequence (either exactly or as a minimum threshold).

Key Features of the New SequenceOptimizer
-----------------------------------------

The ``SequenceOptimizer`` has been completely rewritten to provide:

* **Adaptive Property Scaling**: Automatically adjusts optimization focus based on property convergence patterns and error magnitudes. This makes it easier to optimizer towards properties with highly variable scales or difficult optimization landscapes.
* **Diverse Initial Sequences**: If you are generating a completely new sequence, you can specifiy the number of starting sequences to find better optimization.
* **Flexible Constraint Types**: Support for exact matching, minimum thresholds, and maximum constraints for each specified property
* **Per-Property Tolerances**: Set individual error tolerances for each property, allowing fine-grained control
* **Advanced Convergence Detection**: Multiple convergence criteria including error tolerance, trend analysis, and stagnation detection
* **Performance Optimization**: Comprehensive caching system for faster optimization
* **Arbitrary Number of Properties**: Optimizer towards multiple instances of the same property. This was not previously supported.
* **Easier Property Value Setting**: For many of the properties, you can now set the target value using a sequence of interest rather than a numeric value. 
* **Match to arbitrary interaction matrices**: You can now optimize sequences to match arbitrary interaction matrices.
* **Linear Profiles for Values**: You can now set ``can_be_linear_profile=True`` for some properties and provide a sequence or list of target values. The optimizer will then attempt to match the profile along the values. 

Critical Differences between SequenceOptimizer and Create Functionality
-----------------------------------------------------------------------

The ``SequenceOptimizer`` represents a fundamentally different approach to sequence generation compared to the ``create`` module:

* **Flexibility vs. Speed**: ``SequenceOptimizer`` prioritizes extreme flexibility and handles complex multi-property optimization scenarios that would be difficult ``create``. However, for simple, well-defined property targets, ``create`` functions are typically faster.

* **Approximate vs. Exact Solutions**: ``SequenceOptimizer`` returns the best possible sequence within the optimization constraints and may not achieve exact target values. In contrast, ``create`` functions either generate sequences that exactly meet specifications or fail completely.

* **Extensibility**: Adding new properties to ``SequenceOptimizer`` requires only implementing a simple property class. Adding new functionality to ``create`` requires significant backend overhead.

* **Multi-Property Optimization**: ``SequenceOptimizer`` excels at balancing multiple competing properties simultaneously, while ``create`` functions typically handle individual properties or simple property combinations.

.. contents:: Table of Contents
   :local:
   :depth: 2

Quick Start Example
-------------------

Design a sequence of length 50 with a target hydrophobicity:

.. code-block:: python

    import goose
    from sparrow import Protein

    # Initialize optimizer with basic parameters
    optimizer = goose.SequenceOptimizer(
        target_length=50,
        max_iterations=1000,
        verbose=True
    )
    
    # Add hydrophobicity property with a tolerance
    optimizer.add_property(
        goose.Hydrophobicity, 
        target_value=0.5, 
        weight=1.0,
        tolerance=0.05  # Allow 5% deviation
    )
    
    # Run optimization
    optimized_sequence = optimizer.run()
    
    # Analyze results
    final_protein = Protein(optimized_sequence)
    print(f"Optimized Sequence: {optimized_sequence}")
    print(f"Final Hydrophobicity: {final_protein.hydrophobicity:.3f}")
    print(f"Target Hydrophobicity: 0.5 ± 0.05")

**Explanation:**
- ``SequenceOptimizer(target_length=50, max_iterations=1000, verbose=True)``: Creates optimizer with sequence length, iteration limit, and progress reporting.
- ``add_property(..., tolerance=0.05)``: Adds hydrophobicity optimization with 5% error tolerance.
- ``run()``: Executes optimization with adaptive scaling and convergence detection.

**Advanced Quick Start with Multiple Properties:**

.. code-block:: python

    import goose

    optimizer = goose.SequenceOptimizer(target_length=100, verbose=True)
    
    # Exact hydrophobicity target
    optimizer.add_property(
        goose.Hydrophobicity, 
        target_value=2.4, 
        weight=1.0,
    )
    
    # Minimum disorder requirement
    optimizer.add_property(
        goose.FractionDisorder, 
        target_value=0.8, 
        weight=2.0,  # Higher weight = more important
        constraint_type='minimum',
        disorder_cutoff=0.5
    )
    
    # Maximum FCR constraint
    optimizer.add_property(
        goose.FCR, 
        target_value=0.3, 
        weight=1.5,
        constraint_type='maximum'
    )
    
    optimized_sequence = optimizer.run()

Property Classes Overview
-------------------------

All property classes support three constraint types and individual tolerances:

* **exact**: Minimize absolute difference from target (default)
* **minimum**: Penalize only when below target value
* **maximum**: Penalize only when above target value

To specify constraint type, use the ``constraint_type`` argument when adding a property:

.. code-block:: python
    
    # Exact target (default)
    optimizer.add_property(goose.Hydrophobicity, target_value=0.5, constraint_type='exact')
    
    # Minimum requirement
    optimizer.add_property(goose.FractionDisorder, target_value=0.8, constraint_type='minimum')

    # Maximum constraint
    optimizer.add_property(goose.FCR, target_value=0.3, constraint_type='maximum')

**Basic Properties**

+-------------------------------+-----------------------------------------------+------------------------------------------------+
| Property Class                | Description                                   | Key Arguments                                  |
+===============================+===============================================+================================================+
| Hydrophobicity                | Average hydrophobicity (0-9.0 scale)          | target_value, weight, constraint_type          |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| FCR                           | Fraction of Charged Residues (0-1)            | target_value, weight, constraint_type          |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| NCPR                          | Net Charge Per Residue (-1 to 1)              | target_value, weight, constraint_type          |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| Kappa                         | Charge patterning parameter (0-1)             | target_value, weight, constraint_type          |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| SCD                           | Sequence Charge Decoration                    | target_value, weight, constraint_type          |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| SHD                           | Sequence Hydropathy Decoration                | target_value, weight, constraint_type          |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| Complexity                    | Wootton-Federhen complexity                   | target_value, weight, constraint_type          |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| ComputeIWD                    | Inverse Weighted Distance                     | residues (tuple), target_value, weight,        |
|                               |                                               | constraint_type                                |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| AminoAcidFractions            | Target amino acid composition                 | target_fractions (dict), weight,               |
|                               |                                               | constraint_type                                |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| MatchingResidues              | Number of matching residues to target         | target_sequence, target_value, weight,         |
|                               |                                               | constraint_type                                |
+-------------------------------+-----------------------------------------------+------------------------------------------------+

**Ensemble Properties**

+-------------------------------+-----------------------------------------------+------------------------------------------------+
| Property Class                | Description                                   | Key Arguments                                  |
+===============================+===============================================+================================================+
| RadiusOfGyration              | Predicted radius of gyration (A)              | target_value, weight, constraint_type          |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| EndToEndDistance              | Predicted end-to-end distance (A)             | target_value, weight, constraint_type          |
+-------------------------------+-----------------------------------------------+------------------------------------------------+


**Disorder**

+-------------------------------+-----------------------------------------------+------------------------------------------------+
| Property Class                | Description                                   | Key Arguments                                  |
+===============================+===============================================+================================================+
| FractionDisorder              | Fraction of disordered residues (0-1)         | target_value, weight, constraint_type,         |
|                               |                                               | disorder_cutoff                                |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| MatchSequenceDisorder         | Match disorder profile of target sequence     | target_sequence, weight, constraint_type,      |
|                               |                                               | exact_match, target_value                      |
+-------------------------------+-----------------------------------------------+------------------------------------------------+


**Interaction Properties (Epsilon-based)**

+------------------------+----------------------------------------+-------------------------------------------+
| Property Class         | Description                            | Key Arguments                             |
+========================+========================================+===========================================+
|| MeanSelfEpsilon       || Self-interaction potential            || target_value, weight,                    |
||                       ||                                       || preloaded_model, constraint_type, model  |
+------------------------+----------------------------------------+-------------------------------------------+
|| MeanEpsilonWithTarget || Mean interaction with target sequence || target_value, target_sequence, weight,   |
||                       ||                                       || constraint_type, model, preloaded_model  |
+------------------------+----------------------------------------+-------------------------------------------+
|| ChemicalFingerprint   || Match chemical fingerprint to target  || target_sequence, target_value, weight,   |
||                       ||                                       || constraint_type, model, preloaded_model, |
||                       ||                                       || window_size                              |
+------------------------+----------------------------------------+-------------------------------------------+

**Matrix-based Interaction Properties**

+-------------------------------+-----------------------------------------------+------------------------------------------------+
| Property Class                | Description                                   | Key Arguments                                  |
+===============================+===============================================+================================================+
| MatchSelfIntermap             | Match self-interaction matrix                 | sequence, weight, constraint_type, model,      |
|                               |                                               | preloaded_model, inverse, window_size,         |
|                               |                                               | allow_matrix_resizing                          |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| MatchIntermap                 | Match interaction matrix with target          | sequence, target_sequence, weight,             |
|                               |                                               | constraint_type, model, preloaded_model,       |
|                               |                                               | window_size, allow_matrix_resizing             |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| ModifyAttractiveValues        | Modify attractive interactions                | sequence, target_sequence, multiplier,         |
|                               |                                               | weight, constraint_type, model,                |
|                               |                                               | preloaded_model, window_size                   |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| ModifyRepulsiveValues         | Modify repulsive interactions                 | interacting_sequence,                          |
|                               |                                               | target_interacting_sequence, multiplier,       |
|                               |                                               | weight, constraint_type, model,                |
|                               |                                               | preloaded_model, window_size                   |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| ModifyMatrixValues            | Modify both attractive and repulsive          | interacting_sequence,                          |
|                               |                                               | target_interacting_sequence,                   |
|                               |                                               | repulsive_multiplier, attractive_multiplier,   |
|                               |                                               | weight, constraint_type, model,                |
|                               |                                               | preloaded_model, window_size                   |
+-------------------------------+-----------------------------------------------+------------------------------------------------+

**Folded Domain Surface Properties**

+-------------------------------+-----------------------------------------------+------------------------------------------------+
| Property Class                | Description                                   | Key Arguments                                  |
+===============================+===============================================+================================================+
| FDMeanSurfaceEpsilon          | Mean surface epsilon for folded domains       | target_value, weight, constraint_type, model,  |
|                               |                                               | preloaded_model, path_to_pdb, probe_radius,    |
|                               |                                               | surface_thresh, sasa_mode, fd_start, fd_end,   |
|                               |                                               | preloaded_fd                                   |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| FDSurfaceEpsilon              | Surface epsilon interactions                  | repulsive_target, attractive_target, weight,   |
|                               |                                               | constraint_type, model, preloaded_model,       |
|                               |                                               | path_to_pdb, probe_radius, surface_thresh,     |
|                               |                                               | sasa_mode, fd_start, fd_end, preloaded_fd      |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| FDSurfacePatchInteractions    | Surface patch interaction analysis            | target_value, weight, constraint_type, model,  |
|                               |                                               | preloaded_model, path_to_pdb, probe_radius,    |
|                               |                                               | surface_thresh, sasa_mode, fd_start, fd_end,   |
|                               |                                               | preloaded_fd, patch_residues                   |
+-------------------------------+-----------------------------------------------+------------------------------------------------+

Optimizer Initialization and Basic Parameters
-------------------------------------------------

The ``SequenceOptimizer`` provides extensive control over the optimization process through initialization parameters. You can see additional parameters to change in the Advanced Optimizer Configuration section below.

**Basic Parameters:**

.. code-block:: python

    optimizer = goose.SequenceOptimizer(
        target_length=100,              # Required: target sequence length
        max_iterations=1000,            # Maximum optimization iterations
        verbose=True                   # Enable progress reporting
    )

**Mutation and Diversity Parameters:**

.. code-block:: python

    optimizer = goose.SequenceOptimizer(
        target_length=100,
        # Candidate generation
        num_candidates=5,               # Candidate sequences per iteration
        num_starting_candidates=100,    # Number of sequences to start with. 
        min_mutations=1,                # Minimum mutations per candidate
        max_mutations=15,               # Maximum mutations per candidate
        mutation_ratio=10,              # Length divisor for mutation calculation
        
        # Shuffling for diversity
        enable_shuffling=True,          # Enable sequence shuffling
        shuffle_frequency=50,           # Shuffle every N iterations
        global_shuffle_probability=0.4, # Probability of global vs local shuffle
        shuffle_window_size=15          # Window size for local shuffling
    )


**Setting Initial Sequences:**

.. code-block:: python

    # Start from a specific sequence
    initial_seq = "MGSWAEFKQRLAAIKTRLQALGSQAGKKDAE" * 3  # Must match target_length
    optimizer.set_initial_sequence(initial_seq)
    
    # The optimizer will automatically calculate normalization factors
    # based on the initial sequence for adaptive scaling

Multiple Properties, Weights, and Tolerances
--------------------------------------------

The optimizer excels at balancing multiple competing properties simultaneously. Each property can have individual weights, tolerances, and constraint types:

.. code-block:: python

    import goose
    from sparrow import Protein

    # Create optimizer with advanced parameters
    optimizer = goose.SequenceOptimizer(
        target_length=100, 
        max_iterations=2000,
        verbose=True
    )

    # Critical property - must be close to target
    optimizer.add_property(
        goose.FractionDisorder, 
        target_value=0.85, 
        weight=3.0,                    # High importance
        tolerance=0.02,                # Very strict tolerance (2%)
        constraint_type='minimum'  # Must be at least 85% disordered
    )

    # Important but flexible property
    optimizer.add_property(
        goose.FCR, 
        target_value=0.4, 
        weight=2.0,                    # Medium-high importance
        tolerance=0.05,                # 5% tolerance
    )

    # Secondary property - more flexible
    optimizer.add_property(
        goose.NCPR, 
        target_value=-0.1, 
        weight=1.0,                    # Lower importance
        tolerance=0.1                 # 10% tolerance - quite flexible
    )

    # Compositional constraint
    optimizer.add_property(
        goose.AminoAcidFractions,
        target_fractions={'G': 0.15, 'P': 0.10, 'S': 0.12},
        weight=1.5,
        tolerance=0.03                # 3% tolerance on each amino acid
    )

    # Run optimization
    optimized_sequence = optimizer.run()

    # Analyze results
    final_protein = Protein(optimized_sequence)
    print(f"Optimized Sequence: {optimized_sequence}")
    print(f"Final FCR: {final_protein.FCR:.3f} (target: 0.4 ± 0.05)")
    print(f"Final NCPR: {final_protein.NCPR:.3f} (target: -0.1 ± 0.1)")
    fracs=final_protein.amino_acid_fractions
    print(f"Final fractions: G = {fracs['G']:.3f}, P = {fracs['P']:.3f}, S = {fracs['S']:.3f},")    


Custom Properties
-----------------

Creating custom properties is straightforward by subclassing ``CustomProperty``. The new system supports all constraint types and tolerances automatically:

.. code-block:: python

    import goose
    from goose.backend.optimizer_properties import CustomProperty, ConstraintType
    import sparrow

    class AlanineCount(CustomProperty):
        """Count the number of alanine residues in the sequence."""
        
        def __init__(self, target_value: float, weight: float = 1.0, 
                     constraint_type: ConstraintType = ConstraintType.EXACT):
            super().__init__(target_value, weight, constraint_type)
        
        def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
            """Calculate the raw property value (before constraint application)."""
            return float(protein.sequence.count('A'))

    class MotifCount(CustomProperty):
        """Count occurrences of a specific motif in the sequence."""
        
        def __init__(self, motif: str, target_value: float, weight: float = 1.0,
                     constraint_type: ConstraintType = ConstraintType.EXACT):
            super().__init__(target_value, weight, constraint_type)
            self.motif = motif
        
        def get_init_args(self) -> dict:
            """Override to include motif parameter for serialization."""
            return {
                "motif": self.motif,
                "target_value": self.target_value,
                "weight": self.weight,
                "constraint_type": self.constraint_type.value
            }
        
        def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
            sequence = protein.sequence
            count = 0
            start = 0
            while True:
                pos = sequence.find(self.motif, start)
                if pos == -1:
                    break
                count += 1
                start = pos + 1
            return float(count)

**Using Custom Properties:**

.. code-block:: python

    # Create optimizer
    optimizer = goose.SequenceOptimizer(target_length=100, verbose=True)

    # Add custom properties with different constraint types
    optimizer.add_property(
        AlanineCount, 
        target_value=12.0, 
        weight=1.0,
        constraint_type='exact',
        tolerance=1.0  # Allow ±1 alanine
    )

    optimizer.add_property(
        MotifCount,
        motif="GPG",
        target_value=3.0,  # Want exactly 3 GPG motifs
        weight=2.0,
        constraint_type='exact',
        tolerance=0.0  # Must be exact
    )

    # Standard properties
    optimizer.add_property(
        goose.FractionDisorder,
        target_value=0.8,
        weight=3.0,
        constraint_type='minimum',
    )

    # Run optimization
    optimized_sequence = optimizer.run()

    # Analyze results
    final_protein = sparrow.Protein(optimized_sequence)
    print(f"Optimized Sequence: {optimized_sequence}")
    print(f"Alanine count: {optimized_sequence.count('A')}")
    print(f"GPG motifs: {optimized_sequence.count('GPG')}")


.. note::
   **Best Practices for Custom Properties:**
   
   - Always implement ``calculate_raw_value()`` instead of ``calculate()``
   - Use ``get_init_args()`` if your property has additional parameters
   - The base class automatically handles constraint types and tolerances


Advanced Optimizer Configuration    
---------------------------------

Below are additional parameters to customize the optimization process. You can set these during initialization or modify them later using dedicated methods.
The default parameter values are chosen to provide robust performance across a wide range of scenarios. However, you can adjust them to better suit your specific optimization needs.

**Convergence and Tolerance Controls:**

.. code-block:: python

    optimizer = goose.SequenceOptimizer(
        target_length=100,
        # Error tolerance stopping
        error_tolerance=1e-6,           # Stop when total error below this value
        enable_error_tolerance=True,    # Enable error tolerance early stopping
        
        # Convergence detection
        convergence_tolerance=1e-4,     # Convergence criterion for early stopping
        convergence_window=20,          # Number of recent iterations to check
        enable_early_convergence=False, # Enable early stopping on convergence
        convergence_patience=20,        # Wait iterations after convergence
        
        # Stagnation detection
        stagnation_threshold=25,        # Iterations before considering stagnant
        stagnation_improvement_threshold=0.005  # Minimum improvement to avoid stagnation
    )


**Adaptive Scaling Parameters:**

.. code-block:: python

    optimizer = goose.SequenceOptimizer(
        target_length=100,
        # Adaptive scaling control
        enable_adaptive_scaling=True,   # Enable adaptive property scaling
        max_distance_factor=3.0,        # Maximum scaling based on distance
        distance_offset=0.2,            # Offset for distance calculation
        boost_factor=2.0,               # Factor to boost underperforming properties
        scale_momentum=0.5,             # Momentum for scale smoothing (0-1)
        scale_learning_rate=0.5,        # Learning rate for scale updates (0-1)
        min_scale=0.1,                  # Minimum allowed property scale
        max_scale=8.0,                  # Maximum allowed property scale
        
        # Thresholds for adaptive behavior
        low_contribution_threshold=0.15, # Threshold for low-contributing properties
        high_error_threshold=0.05,      # Threshold for high-error properties
        stagnation_multiplier=1.0       # Multiplier for stagnation response
    )

**History and Analysis Parameters:**

.. code-block:: python

    optimizer = goose.SequenceOptimizer(
        target_length=100,
        # History tracking
        improvement_history_size=20,    # Recent improvements per property
        error_history_size=50,          # Recent error values to store
        
        # Analysis parameters
        min_analysis_samples=5,         # Minimum samples for analysis
        min_trend_samples=5,            # Minimum samples for trend calculation
        improvement_threshold=-0.001,   # Threshold for improvement detection
        stability_threshold=0.01,       # Variance threshold for stability
        
        # Progress reporting
        update_interval=10              # Update progress every N iterations
    )

**Dynamic Configuration Methods:**

You can modify convergence and error tolerance settings after initialization:

.. code-block:: python

    # Configure convergence detection
    optimizer.configure_convergence(
        tolerance=1e-5,                 # New convergence tolerance
        window=30,                      # New convergence window
        enable_early_stopping=True,     # Enable early stopping
        patience=15                     # New patience value
    )
    
    # Configure error tolerance
    optimizer.configure_error_tolerance(
        tolerance=1e-7,                 # New error tolerance
        enable=True                     # Enable/disable error tolerance stopping
    )
    
    # Get convergence information
    convergence_info = optimizer.get_convergence_info()
    print(f"Convergence status: {convergence_info}")



Troubleshooting and Optimization Tips
-------------------------------------

**Optimization Not Converging**

*Symptoms*: Error plateaus at high values, properties far from targets

*Solutions*:
- **Increase iterations**: ``max_iterations=5000`` or higher for complex problems
- **Enable adaptive scaling**: ``enable_adaptive_scaling=True`` (default)
- **Increase diversity**: ``shuffle_frequency=25``, ``num_candidates=10``
- **Check target compatibility**: Ensure properties don't fundamentally conflict
- **Use tolerances**: Set reasonable ``tolerance`` values for each property
- **Verify constraint types**: Make sure you're using appropriate constraints

**Slow Optimization Performance**

*Symptoms*: Optimization takes too long, high memory usage

*Solutions*:
- **Reduce candidates**: ``num_candidates=3`` for faster iterations (default is 5)
- **Disable expensive features**: ``enable_adaptive_scaling=False``, ``enable_shuffling=False``
- **Use stricter early stopping**: ``error_tolerance=1e-4``, ``enable_early_convergence=True``
- **Optimize caching**: Check cache hit rate with ``get_cache_statistics()``
- **Pre-load models**: Use ``preloaded_model`` for epsilon properties

**Property Conflicts and Balancing**

*Symptoms*: Some properties optimize while others get worse

*Solutions*:
- **Adjust weights**: Higher weight = higher priority
- **Use appropriate constraint types**: MINIMUM/MAXIMUM instead of EXACT when possible
- **Set generous tolerances**: Allow some flexibility in less critical properties
- **Check physical compatibility**: Some combinations may be impossible
- **Monitor individual properties**: Enable ``verbose=True`` to track individual progress

.. code-block:: python

    # Balanced multi-property optimization
    optimizer.add_property(goose.FractionDisorder, target_value=0.8, weight=3.0, 
                          constraint_type=ConstraintType.MINIMUM, tolerance=0.05)
    optimizer.add_property(goose.FCR, target_value=0.3, weight=1.0, 
                          constraint_type=ConstraintType.EXACT, tolerance=0.1)
    optimizer.add_property(goose.Hydrophobicity, target_value=0.4, weight=0.5, 
                          constraint_type=ConstraintType.EXACT, tolerance=0.2)

**Memory Issues with Large Sequences**

*Symptoms*: Out of memory errors, excessive RAM usage

*Solutions*:
- **Reduce history sizes**: ``improvement_history_size=5``, ``error_history_size=10``
- **Clear cache periodically**: Call ``optimizer._clear_evaluation_cache()`` if needed
- **Disable caching**: Set caching parameters conservatively
- **Use fewer candidates**: ``num_candidates=3`` for large sequences

.. code-block:: python

    # Memory-efficient settings for large sequences
    optimizer = goose.SequenceOptimizer(
        target_length=1000,
        improvement_history_size=5,
        error_history_size=10,
        num_candidates=3,
        debugging=False
    )

**Stagnation Issues**

*Symptoms*: Error doesn't improve for many iterations

*Solutions*:
- **Enable shuffling**: ``enable_shuffling=True`` with frequent shuffling
- **Adjust stagnation detection**: Lower ``stagnation_threshold=15``
- **Increase mutation diversity**: Higher ``max_mutations=20``
- **Check for impossible targets**: Some property combinations may be unachievable

Examples and Demo Notebooks
----------------------------

GOOSE includes comprehensive demo notebooks showcasing advanced ``SequenceOptimizer`` usage:

**Available Demos:**
- **Basic optimization**: see /demos/sequence_optimization.ipynb for basic usage. 
- **Custom properties**: Creating and implementing user-defined properties  
- **Interaction optimization**: Epsilon-based properties and matrix manipulations
- **Performance optimization**: Efficient settings for different use cases
- **Advanced features**: Adaptive scaling, convergence detection, and troubleshooting

**Demo Location:**
Check the ``demos`` directory for Jupyter notebooks with detailed examples and explanations.

API Reference
-------------

**Core Classes:**
- ``goose.SequenceOptimizer``: Main optimization engine
- ``goose.backend.optimizer_properties.ProteinProperty``: Base class for properties
- ``goose.backend.optimizer_properties.CustomProperty``: Base class for custom properties users can define

**Key Methods:**
- ``SequenceOptimizer.add_property()``: Add properties to optimize
- ``SequenceOptimizer.set_initial_sequence()``: Set starting sequence
- ``SequenceOptimizer.run()``: Execute optimization

See Also
--------

For complete API documentation, see ``goose/optimize.py`` and ``goose/backend/optimizer_properties.py``.

For implementation examples and advanced usage patterns, explore the demo notebooks in ``demos/``.
