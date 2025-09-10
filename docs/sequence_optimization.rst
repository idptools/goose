Using the SequenceOptimizer
==============================

**SUPER DUPER IMPORTANT NOTE**: GOOSE is an IDR design tool. **HOWEVER**, when using SequenceOptimizer, you can design anything you want. Thus, sequences are not guaranteed to be disordered unless you specify the ``FractionDisorder`` property. You can also use the ``MatchSequenceDisorder`` property if you'd like to make a variant that maintains the disorder profile of your original sequence (either exactly or as a minimum threshold).

GOOSE's ``SequenceOptimizer`` is a flexible tool for designing protein sequences that match user-defined biophysical or biochemical property targets. It uses advanced stochastic optimization with adaptive scaling, intelligent shuffling, and convergence detection to efficiently explore sequence space and minimize the difference between calculated and target property values. You can simultaneously optimize towards arbitrary numbers of properties with individual weights, tolerances, and constraint types.

Key Features of the New SequenceOptimizer
-----------------------------------------

The ``SequenceOptimizer`` has been completely rewritten to provide:

* **Adaptive Property Scaling**: Automatically adjusts optimization focus based on property convergence patterns and error magnitudes
* **Intelligent Shuffling**: Combines global and local sequence shuffling to escape local minima and maintain diversity
* **Flexible Constraint Types**: Support for exact matching, minimum thresholds, and maximum constraints per property
* **Per-Property Tolerances**: Set individual error tolerances for each property, allowing fine-grained control
* **Advanced Convergence Detection**: Multiple convergence criteria including error tolerance, trend analysis, and stagnation detection
* **Performance Optimization**: Comprehensive caching system and vectorized operations for faster optimization
* **Extensive Monitoring**: Detailed progress tracking, error analysis, and convergence diagnostics

Critical Differences between SequenceOptimizer and Create Functionality
-----------------------------------------------------------------------

The ``SequenceOptimizer`` represents a fundamentally different approach to sequence generation compared to the ``create`` module:

* **Flexibility vs. Speed**: ``SequenceOptimizer`` prioritizes extreme flexibility and handles complex multi-property optimization scenarios that would be impossible with ``create``. However, for simple, well-defined property targets, ``create`` functions are typically faster.

* **Approximate vs. Exact Solutions**: ``SequenceOptimizer`` returns the best possible sequence within the optimization constraints and may not achieve exact target values. In contrast, ``create`` functions either generate sequences that exactly meet specifications or fail completely.

* **Extensibility**: Adding new properties to ``SequenceOptimizer`` requires only implementing a simple property class. Adding new functionality to ``create`` requires significant backend development.

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
    from goose.backend.optimizer_properties import ConstraintType

    optimizer = goose.SequenceOptimizer(target_length=100, verbose=True)
    
    # Exact hydrophobicity target
    optimizer.add_property(
        goose.Hydrophobicity, 
        target_value=0.4, 
        weight=1.0,
        constraint_type=ConstraintType.EXACT
    )
    
    # Minimum disorder requirement
    optimizer.add_property(
        goose.FractionDisorder, 
        target_value=0.8, 
        weight=2.0,  # Higher weight = more important
        constraint_type=ConstraintType.MINIMUM,
        disorder_cutoff=0.5
    )
    
    # Maximum FCR constraint
    optimizer.add_property(
        goose.FCR, 
        target_value=0.3, 
        weight=1.5,
        constraint_type=ConstraintType.MAXIMUM
    )
    
    optimized_sequence = optimizer.run()

Property Classes Overview
-------------------------

All property classes support three constraint types and individual tolerances:

* **ConstraintType.EXACT**: Minimize absolute difference from target (default)
* **ConstraintType.MINIMUM**: Penalize only when below target value
* **ConstraintType.MAXIMUM**: Penalize only when above target value

.. code-block:: python

    from goose.backend.optimizer_properties import ConstraintType
    
    # Exact target (default)
    optimizer.add_property(goose.Hydrophobicity, target_value=0.5, constraint_type=ConstraintType.EXACT)
    
    # Minimum requirement
    optimizer.add_property(goose.FractionDisorder, target_value=0.8, constraint_type=ConstraintType.MINIMUM)
    
    # Maximum constraint
    optimizer.add_property(goose.FCR, target_value=0.3, constraint_type=ConstraintType.MAXIMUM)

**Basic Physicochemical Properties**

+-------------------------------+-----------------------------------------------+------------------------------------------------+
| Property Class                | Description                                   | Key Arguments                                  |
+===============================+===============================================+================================================+
| Hydrophobicity                | Average hydrophobicity (0-6.6 scale)         | target_value, weight, constraint_type          |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| FCR                           | Fraction of Charged Residues (0-1)           | target_value, weight, constraint_type          |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| NCPR                          | Net Charge Per Residue (-1 to 1)             | target_value, weight, constraint_type          |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| Kappa                         | Charge patterning parameter (0-1)            | target_value, weight, constraint_type          |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| SCD                           | Sequence Charge Decoration                   | target_value, weight, constraint_type          |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| SHD                           | Sequence Hydropathy Decoration               | target_value, weight, constraint_type          |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| Complexity                    | Wootton-Federhen (SEG) complexity            | target_value, weight, constraint_type          |
+-------------------------------+-----------------------------------------------+------------------------------------------------+

**Structural Properties**

+-------------------------------+-----------------------------------------------+------------------------------------------------+
| Property Class                | Description                                   | Key Arguments                                  |
+===============================+===============================================+================================================+
| RadiusOfGyration              | Predicted radius of gyration (Å)             | target_value, weight, constraint_type          |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| EndToEndDistance              | Predicted end-to-end distance (Å)            | target_value, weight, constraint_type          |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| ComputeIWD                    | Inverse Weighted Distance                     | residues (tuple), target_value, weight,        |
|                               |                                               | constraint_type                                |
+-------------------------------+-----------------------------------------------+------------------------------------------------+

**Disorder and Composition Properties**

+-------------------------------+-----------------------------------------------+------------------------------------------------+
| Property Class                | Description                                   | Key Arguments                                  |
+===============================+===============================================+================================================+
| FractionDisorder              | Fraction of disordered residues (0-1)        | target_value, weight, constraint_type,         |
|                               |                                               | disorder_cutoff                                |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| AminoAcidFractions            | Target amino acid composition                 | target_fractions (dict), weight,               |
|                               |                                               | constraint_type                                |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| MatchSequenceDisorder         | Match disorder profile of target sequence    | target_sequence, weight, constraint_type,      |
|                               |                                               | exact_match, target_value                      |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| MatchingResidues              | Number of matching residues to target        | target_sequence, target_value, weight,         |
|                               |                                               | constraint_type                                |
+-------------------------------+-----------------------------------------------+------------------------------------------------+

**Interaction Properties (Epsilon-based)**

+-------------------------------+-----------------------------------------------+------------------------------------------------+
| Property Class                | Description                                   | Key Arguments                                  |
+===============================+===============================================+================================================+
| MeanSelfEpsilon               | Self-interaction potential                    | target_value, weight, constraint_type, model, |
|                               |                                               | preloaded_model                                |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| MeanEpsilonWithTarget         | Mean interaction with target sequence        | target_value, target_sequence, weight,         |
|                               |                                               | constraint_type, model, preloaded_model       |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| ChemicalFingerprint           | Match chemical fingerprint to target         | target_sequence, target_value, weight,         |
|                               |                                               | constraint_type, model, preloaded_model,      |
|                               |                                               | window_size                                    |
+-------------------------------+-----------------------------------------------+------------------------------------------------+

**Matrix-based Interaction Properties**

+-------------------------------+-----------------------------------------------+------------------------------------------------+
| Property Class                | Description                                   | Key Arguments                                  |
+===============================+===============================================+================================================+
| MatchSelfIntermap             | Match self-interaction matrix                 | sequence, weight, constraint_type, model,      |
|                               |                                               | preloaded_model, inverse, window_size,        |
|                               |                                               | allow_matrix_resizing                          |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| MatchIntermap                 | Match interaction matrix with target         | sequence, target_sequence, weight,             |
|                               |                                               | constraint_type, model, preloaded_model,      |
|                               |                                               | window_size, allow_matrix_resizing             |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| ModifyAttractiveValues        | Modify attractive interactions                | sequence, target_sequence, multiplier,         |
|                               |                                               | weight, constraint_type, model,                |
|                               |                                               | preloaded_model, window_size                   |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| ModifyRepulsiveValues         | Modify repulsive interactions                 | interacting_sequence,                          |
|                               |                                               | target_interacting_sequence, multiplier,      |
|                               |                                               | weight, constraint_type, model,                |
|                               |                                               | preloaded_model, window_size                   |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| ModifyMatrixValues            | Modify both attractive and repulsive          | interacting_sequence,                          |
|                               |                                               | target_interacting_sequence,                   |
|                               |                                               | repulsive_multiplier, attractive_multiplier,  |
|                               |                                               | weight, constraint_type, model,                |
|                               |                                               | preloaded_model, window_size                   |
+-------------------------------+-----------------------------------------------+------------------------------------------------+

**Folded Domain Surface Properties**

+-------------------------------+-----------------------------------------------+------------------------------------------------+
| Property Class                | Description                                   | Key Arguments                                  |
+===============================+===============================================+================================================+
| FDMeanSurfaceEpsilon          | Mean surface epsilon for folded domains      | target_value, weight, constraint_type, model, |
|                               |                                               | preloaded_model, path_to_pdb, probe_radius,   |
|                               |                                               | surface_thresh, sasa_mode, fd_start, fd_end,  |
|                               |                                               | preloaded_fd                                   |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| FDSurfaceEpsilon              | Surface epsilon interactions                  | repulsive_target, attractive_target, weight,  |
|                               |                                               | constraint_type, model, preloaded_model,      |
|                               |                                               | path_to_pdb, probe_radius, surface_thresh,    |
|                               |                                               | sasa_mode, fd_start, fd_end, preloaded_fd     |
+-------------------------------+-----------------------------------------------+------------------------------------------------+
| FDSurfacePatchInteractions    | Surface patch interaction analysis           | target_value, weight, constraint_type, model, |
|                               |                                               | preloaded_model, path_to_pdb, probe_radius,   |
|                               |                                               | surface_thresh, sasa_mode, fd_start, fd_end,  |
|                               |                                               | preloaded_fd, patch_residues                   |
+-------------------------------+-----------------------------------------------+------------------------------------------------+

Optimizer Initialization and Parameters
---------------------------------------

The ``SequenceOptimizer`` provides extensive control over the optimization process through initialization parameters:

**Basic Parameters:**

.. code-block:: python

    optimizer = goose.SequenceOptimizer(
        target_length=100,              # Required: target sequence length
        max_iterations=1000,            # Maximum optimization iterations
        verbose=True,                   # Enable progress reporting
        debugging=False                 # Enable detailed debugging output
    )

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

**Mutation and Diversity Parameters:**

.. code-block:: python

    optimizer = goose.SequenceOptimizer(
        target_length=100,
        # Candidate generation
        num_candidates=5,               # Candidate sequences per iteration
        min_mutations=1,                # Minimum mutations per candidate
        max_mutations=15,               # Maximum mutations per candidate
        mutation_ratio=10,              # Length divisor for mutation calculation
        
        # Shuffling for diversity
        enable_shuffling=True,          # Enable sequence shuffling
        shuffle_frequency=50,           # Shuffle every N iterations
        global_shuffle_probability=0.4, # Probability of global vs local shuffle
        shuffle_window_size=15          # Window size for local shuffling
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
    from goose.backend.optimizer_properties import ConstraintType
    from sparrow import Protein

    # Create optimizer with advanced parameters
    optimizer = goose.SequenceOptimizer(
        target_length=100, 
        max_iterations=2000,
        verbose=True,
        enable_adaptive_scaling=True  # Automatically balance property importance
    )

    # Critical property - must be close to target
    optimizer.add_property(
        goose.FractionDisorder, 
        target_value=0.85, 
        weight=3.0,                    # High importance
        tolerance=0.02,                # Very strict tolerance (2%)
        constraint_type=ConstraintType.MINIMUM  # Must be at least 85% disordered
    )

    # Important but flexible property
    optimizer.add_property(
        goose.FCR, 
        target_value=0.4, 
        weight=2.0,                    # Medium-high importance
        tolerance=0.05,                # 5% tolerance
        constraint_type=ConstraintType.EXACT
    )

    # Secondary property - more flexible
    optimizer.add_property(
        goose.NCPR, 
        target_value=-0.1, 
        weight=1.0,                    # Lower importance
        tolerance=0.1,                 # 10% tolerance - quite flexible
        constraint_type=ConstraintType.EXACT
    )

    # Compositional constraint
    optimizer.add_property(
        goose.AminoAcidFractions,
        target_fractions={'G': 0.15, 'P': 0.10, 'S': 0.12},
        weight=1.5,
        tolerance=0.03,                # 3% tolerance on each amino acid
        constraint_type=ConstraintType.EXACT
    )

    # Run optimization
    optimized_sequence = optimizer.run()

    # Analyze results
    final_protein = Protein(optimized_sequence)
    print(f"Optimized Sequence: {optimized_sequence}")
    print(f"Final Disorder Fraction: {final_protein.disorder_fraction:.3f} (target: ≥0.85)")
    print(f"Final FCR: {final_protein.FCR:.3f} (target: 0.4 ± 0.05)")
    print(f"Final NCPR: {final_protein.NCPR:.3f} (target: -0.1 ± 0.1)")

**Understanding Error Calculation with Tolerances:**

.. code-block:: python

    # Property with tolerance=0.05 and target_value=0.4
    # Current value = 0.42
    # Raw error = |0.42 - 0.4| = 0.02
    # Since 0.02 <= 0.05 (tolerance), effective error = 0.0
    # This property is "satisfied" and contributes 0 to total error

    # Property with tolerance=0.0 (no tolerance) and target_value=0.3
    # Current value = 0.31  
    # Raw error = |0.31 - 0.3| = 0.01
    # Since tolerance=0.0, effective error = 0.01
    # This property contributes weight * 0.01 to total error

**Constraint Types in Detail:**

.. code-block:: python

    # EXACT: Minimize |current - target|
    optimizer.add_property(
        goose.Hydrophobicity, 
        target_value=0.5, 
        constraint_type=ConstraintType.EXACT
    )

    # MINIMUM: Penalize only when current < target
    optimizer.add_property(
        goose.FractionDisorder, 
        target_value=0.8, 
        constraint_type=ConstraintType.MINIMUM  # At least 80% disordered
    )

    # MAXIMUM: Penalize only when current > target  
    optimizer.add_property(
        goose.FCR, 
        target_value=0.3, 
        constraint_type=ConstraintType.MAXIMUM  # At most 30% charged
    )

.. note::
   **Adaptive Scaling**: When enabled, the optimizer automatically adjusts property weights based on convergence patterns, error magnitudes, and stagnation detection. Properties that are harder to optimize or contribute less to overall progress receive dynamic weight adjustments.

.. note::
   **Error Function**: Total error = sum(weight * max(0, effective_error - tolerance)) for each property, where effective_error depends on the constraint type.

Advanced Features and Optimization Strategies
--------------------------------------------

**Intelligent Initial Sequence Selection**

The optimizer can start from diverse initial sequences to find better optimization starting points:

.. code-block:: python

    # Method 1: Set a specific initial sequence
    initial_seq = "MGSWAEFKQRLAAIKTRLQALGSQAGKKDAE" * 3  # Must match target_length
    optimizer.set_initial_sequence(initial_seq)

    # Method 2: Let optimizer build diverse initial sequences internally
    # (This is done automatically in the optimization process)
    
    # The optimizer automatically calculates normalization factors
    # from the initial sequence for adaptive scaling

**Optimization Strategies for Different Scenarios**

**For Fast Optimization (Simple Properties):**

.. code-block:: python

    optimizer = goose.SequenceOptimizer(
        target_length=50,
        max_iterations=500,           # Fewer iterations
        num_candidates=3,             # Fewer candidates per iteration
        enable_adaptive_scaling=False, # Disable adaptive scaling
        enable_shuffling=False,       # Disable shuffling for speed
        verbose=False
    )

**For Challenging Multi-Property Optimization:**

.. code-block:: python

    optimizer = goose.SequenceOptimizer(
        target_length=200,
        max_iterations=5000,          # More iterations
        num_candidates=10,            # More candidates per iteration
        enable_adaptive_scaling=True, # Enable adaptive scaling
        enable_shuffling=True,        # Enable shuffling for diversity
        shuffle_frequency=25,         # More frequent shuffling
        stagnation_threshold=15,      # Detect stagnation sooner
        error_tolerance=1e-8,         # Very strict error tolerance
        verbose=True
    )

**For Epsilon-based Interaction Properties:**

.. code-block:: python

    # Pre-load epsilon model for efficiency
    import finches
    epsilon_model = finches.Load('mpipi')
    
    optimizer.add_property(
        goose.MeanSelfEpsilon,
        target_value=-2.5,
        weight=2.0,
        model='mpipi',
        preloaded_model=epsilon_model  # Avoid reloading model
    )

**Performance Monitoring and Diagnostics**

.. code-block:: python

    # Get cache statistics for performance analysis
    cache_stats = optimizer.get_cache_statistics()
    print(f"Cache hits: {cache_stats['cache_hits']}")
    print(f"Cache misses: {cache_stats['cache_misses']}")
    print(f"Hit rate: {cache_stats['hit_rate']:.2%}")

    # Get convergence information
    convergence_info = optimizer.get_convergence_info()
    print(f"Current error: {convergence_info['current_error']}")
    print(f"Error history: {convergence_info['error_history']}")
    print(f"Convergence detected: {convergence_info['converged']}")

**Working with Complex Properties**

.. code-block:: python

    # Chemical fingerprint matching
    target_sequence = "MGSWAEFKQRLAAIKTRLQALGSQAGKKDAE"
    optimizer.add_property(
        goose.ChemicalFingerprint,
        target_sequence=target_sequence,
        target_value=0.0,  # Perfect match
        weight=2.0,
        model='mpipi',
        window_size=15
    )

    # Matrix-based interaction matching
    optimizer.add_property(
        goose.MatchIntermap,
        sequence="MGSWAE" * 10,  # Original sequence
        target_sequence="FKQRLA" * 10,  # Target interaction partner
        weight=1.5,
        window_size=15,
        allow_matrix_resizing=True
    )

    # Folded domain surface interactions
    optimizer.add_property(
        goose.FDMeanSurfaceEpsilon,
        target_value=-1.5,
        weight=2.0,
        path_to_pdb="/path/to/structure.pdb",
        fd_start=50,  # Start of folded domain
        fd_end=150,   # End of folded domain
        surface_thresh=0.2
    )

**Memory and Performance Optimization**

.. code-block:: python

    # For large sequences or many properties, consider:
    optimizer = goose.SequenceOptimizer(
        target_length=500,
        # Reduce memory usage
        improvement_history_size=10,  # Smaller history
        error_history_size=25,        # Smaller error history
        
        # Optimize for performance
        update_interval=50,           # Less frequent progress updates
        debugging=False,              # Disable debug logging
        
        # Balance exploration vs exploitation
        num_candidates=8,             # Moderate candidate count
        max_mutations=20,             # Allow more mutations for large sequences
        mutation_ratio=15             # Adjust mutation ratio for sequence length
    )

Custom Properties
-----------------

Creating custom properties is straightforward by subclassing ``ProteinProperty``. The new system supports all constraint types and tolerances automatically:

.. code-block:: python

    import goose
    from goose.backend.optimizer_properties import ProteinProperty, ConstraintType
    import sparrow

    class AlanineCount(ProteinProperty):
        """Count the number of alanine residues in the sequence."""
        
        def __init__(self, target_value: float, weight: float = 1.0, 
                     constraint_type: ConstraintType = ConstraintType.EXACT):
            super().__init__(target_value, weight, constraint_type)
        
        def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
            """Calculate the raw property value (before constraint application)."""
            return float(protein.sequence.count('A'))

    class GlycineProlineRatio(ProteinProperty):
        """Calculate the ratio of glycine to proline residues."""
        
        def __init__(self, target_value: float, weight: float = 1.0,
                     constraint_type: ConstraintType = ConstraintType.EXACT):
            super().__init__(target_value, weight, constraint_type)
        
        def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
            sequence = protein.sequence
            g_count = sequence.count('G')
            p_count = sequence.count('P')
            
            # Handle division by zero
            if p_count == 0:
                return float('inf') if g_count > 0 else 0.0
            
            return g_count / p_count

    class MotifCount(ProteinProperty):
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
        constraint_type=ConstraintType.EXACT,
        tolerance=1.0  # Allow ±1 alanine
    )

    optimizer.add_property(
        GlycineProlineRatio, 
        target_value=1.5,  # 1.5 times more G than P
        weight=1.5,
        constraint_type=ConstraintType.MINIMUM,  # At least 1.5:1 ratio
        tolerance=0.1
    )

    optimizer.add_property(
        MotifCount,
        motif="GPG",
        target_value=3.0,  # Want exactly 3 GPG motifs
        weight=2.0,
        constraint_type=ConstraintType.EXACT,
        tolerance=0.0  # Must be exact
    )

    # Standard properties
    optimizer.add_property(
        goose.FractionDisorder,
        target_value=0.8,
        weight=3.0,
        constraint_type=ConstraintType.MINIMUM
    )

    # Run optimization
    optimized_sequence = optimizer.run()

    # Analyze results
    final_protein = sparrow.Protein(optimized_sequence)
    print(f"Optimized Sequence: {optimized_sequence}")
    print(f"Alanine count: {optimized_sequence.count('A')}")
    print(f"G:P ratio: {optimized_sequence.count('G')}/{optimized_sequence.count('P')}")
    print(f"GPG motifs: {optimized_sequence.count('GPG')}")
    print(f"Disorder fraction: {final_protein.disorder_fraction:.3f}")

**Advanced Custom Property with Caching:**

.. code-block:: python

    class ExpensiveProperty(ProteinProperty):
        """Example of a property that benefits from caching expensive calculations."""
        
        def __init__(self, target_value: float, weight: float = 1.0,
                     constraint_type: ConstraintType = ConstraintType.EXACT):
            super().__init__(target_value, weight, constraint_type)
            self._cache = {}  # Internal cache for expensive calculations
        
        def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
            sequence = protein.sequence
            
            # Check cache first
            if sequence in self._cache:
                return self._cache[sequence]
            
            # Expensive calculation here
            # (e.g., secondary structure prediction, contact prediction, etc.)
            result = self._expensive_calculation(sequence)
            
            # Cache the result
            self._cache[sequence] = result
            return result
        
        def _expensive_calculation(self, sequence: str) -> float:
            # Placeholder for expensive computation
            import time
            time.sleep(0.01)  # Simulate expensive calculation
            return len(sequence) * 0.1  # Dummy result

**Property with External Dependencies:**

.. code-block:: python

    class SecondaryStructureContent(ProteinProperty):
        """Calculate secondary structure content using external tool."""
        
        def __init__(self, structure_type: str, target_value: float, 
                     weight: float = 1.0, constraint_type: ConstraintType = ConstraintType.EXACT):
            super().__init__(target_value, weight, constraint_type)
            self.structure_type = structure_type.upper()  # 'H', 'E', 'C'
            
            if self.structure_type not in ['H', 'E', 'C']:
                raise ValueError("structure_type must be 'H' (helix), 'E' (strand), or 'C' (coil)")
        
        def get_init_args(self) -> dict:
            return {
                "structure_type": self.structure_type,
                "target_value": self.target_value,
                "weight": self.weight,
                "constraint_type": self.constraint_type.value
            }
        
        def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
            # Use SPARROW's secondary structure prediction
            ss_prediction = protein.secondary_structure_prediction
            
            # Count fraction of desired structure type
            count = ss_prediction.count(self.structure_type)
            return count / len(ss_prediction)

.. note::
   **Best Practices for Custom Properties:**
   
   - Always implement ``calculate_raw_value()`` instead of ``calculate()``
   - Use ``get_init_args()`` if your property has additional parameters
   - Consider caching for expensive calculations
   - Handle edge cases (division by zero, empty sequences, etc.)
   - The base class automatically handles constraint types and tolerances

Amino Acid Composition Control
-----------------------------

The optimizer uses standard amino acid frequencies for sequence generation by default. Unlike the older k-mer dictionary system, the new optimizer focuses on direct property optimization rather than biased sequence generation.

**Controlling Composition Through Properties:**

Instead of k-mer dictionaries, use composition-specific properties:

.. code-block:: python

    # Direct amino acid fraction control
    optimizer.add_property(
        goose.AminoAcidFractions,
        target_fractions={
            'G': 0.15,  # 15% glycine
            'P': 0.10,  # 10% proline  
            'S': 0.12,  # 12% serine
            'A': 0.08,  # 8% alanine
            'E': 0.10,  # 10% glutamate
            'K': 0.08   # 8% lysine
            # Other amino acids will be distributed among remaining positions
        },
        weight=2.0,
        tolerance=0.02  # Allow 2% deviation per amino acid
    )

    # Charge composition control
    optimizer.add_property(goose.FCR, target_value=0.25, weight=1.5)
    optimizer.add_property(goose.NCPR, target_value=0.1, weight=1.5)

    # Secondary structure amino acid preferences
    # (Disorder-promoting residues)
    disorder_promoting = ['G', 'P', 'S', 'Q', 'A', 'E', 'K', 'R']
    total_disorder_promoting = sum(
        optimizer.properties[0].target_fractions.get(aa, 0.0) 
        for aa in disorder_promoting
    )
    
    # Can also use custom properties for more complex composition rules
    class DisorderPromotingFraction(ProteinProperty):
        def __init__(self, target_value: float, weight: float = 1.0):
            super().__init__(target_value, weight)
        
        def calculate_raw_value(self, protein):
            sequence = protein.sequence
            disorder_aa = 'GPSQAEKR'
            count = sum(sequence.count(aa) for aa in disorder_aa)
            return count / len(sequence)

**Organism-Specific Composition:**

For organism-specific amino acid preferences, use the composition data and create appropriate ``AminoAcidFractions`` properties:

.. code-block:: python

    # Example: E. coli-like composition for disordered regions
    ecoli_idr_composition = {
        'A': 0.082, 'R': 0.045, 'N': 0.040, 'D': 0.062, 'C': 0.015,
        'Q': 0.042, 'E': 0.071, 'G': 0.089, 'H': 0.021, 'I': 0.051,
        'L': 0.087, 'K': 0.063, 'M': 0.024, 'F': 0.033, 'P': 0.055,
        'S': 0.083, 'T': 0.058, 'W': 0.012, 'Y': 0.028, 'V': 0.063
    }
    
    optimizer.add_property(
        goose.AminoAcidFractions,
        target_fractions=ecoli_idr_composition,
        weight=1.0,
        tolerance=0.015  # 1.5% tolerance per amino acid
    )

How the Advanced Optimizer Works
---------------------------------

The new ``SequenceOptimizer`` uses a sophisticated multi-stage optimization process:

**1. Initialization and Normalization**
   - Generates initial sequence (random or user-provided)
   - Calculates normalization factors based on initial property errors
   - Sets up adaptive scaling infrastructure and caching systems

**2. Candidate Generation**
   - Creates multiple candidate sequences per iteration using:
     - Point mutations (1 to max_mutations per candidate)
     - Local and global sequence shuffling
     - Intelligent diversity injection when stagnation is detected

**3. Property Evaluation with Caching**
   - Evaluates all properties for each candidate sequence
   - Uses comprehensive caching to avoid redundant calculations
   - Applies constraint types (EXACT/MINIMUM/MAXIMUM) and tolerances

**4. Adaptive Error Calculation**
   - Calculates weighted error: ``sum(weight * scale * max(0, effective_error - tolerance))``
   - ``effective_error`` depends on constraint type:
     - EXACT: ``|current - target|``
     - MINIMUM: ``max(0, target - current)``
     - MAXIMUM: ``max(0, current - target)``

**5. Adaptive Scaling (if enabled)**
   - Monitors property convergence patterns and error contributions
   - Dynamically adjusts property scales based on:
     - Distance from target values
     - Stagnation detection per property
     - Contribution to total error reduction
     - Improvement trends over recent history

**6. Selection and Best Tracking**
   - Selects candidate with lowest total error
   - Updates best sequence if improvement found
   - Maintains detailed error history and convergence statistics

**7. Convergence and Stagnation Detection**
   - **Error Tolerance**: Stops if total error < ``error_tolerance``
   - **Convergence Detection**: Stops if error change < ``convergence_tolerance`` for ``convergence_window`` iterations
   - **Stagnation Recovery**: Increases shuffling and mutation rates when no improvement for ``stagnation_threshold`` iterations
   - **Emergency Diversity**: Injects random sequences if severe stagnation detected

**8. Advanced Stopping Criteria**
   - Multiple convergence criteria can trigger early stopping
   - Patience mechanisms prevent premature termination
   - Comprehensive diagnostics available throughout optimization

**Optimization Flow Diagram:**

.. code-block:: text

    Initialize → Generate Candidates → Evaluate Properties → Calculate Errors
         ↑              ↓                      ↓                    ↓
    Update Best ← Select Best ← Apply Adaptive Scaling ← Check Convergence
         ↓              ↑                      ↑                    ↓
    Check Stopping ← Detect Stagnation ← Update Scaling ← Continue/Stop?
    
**Key Algorithmic Improvements:**

- **Vectorized Operations**: Batch processing of candidates for efficiency
- **Smart Caching**: Avoids recalculating identical sequences across iterations  
- **Dynamic Balancing**: Automatically focuses on difficult-to-optimize properties
- **Intelligent Diversity**: Context-aware shuffling and mutation strategies
- **Multi-Criteria Stopping**: Robust convergence detection with multiple fallbacks

Troubleshooting and Optimization Tips
------------------------------------

**Optimization Not Converging**

*Symptoms*: Error plateaus at high values, properties far from targets

*Solutions*:
- **Increase iterations**: ``max_iterations=5000`` or higher for complex problems
- **Enable adaptive scaling**: ``enable_adaptive_scaling=True`` (default)
- **Increase diversity**: ``shuffle_frequency=25``, ``num_candidates=10``
- **Check target compatibility**: Ensure properties don't fundamentally conflict
- **Use tolerances**: Set reasonable ``tolerance`` values for each property
- **Verify constraint types**: Make sure you're using appropriate constraints

.. code-block:: python

    # For difficult convergence
    optimizer = goose.SequenceOptimizer(
        target_length=100,
        max_iterations=10000,
        num_candidates=15,
        shuffle_frequency=20,
        stagnation_threshold=15,
        error_tolerance=1e-7,
        enable_adaptive_scaling=True
    )

**Slow Optimization Performance**

*Symptoms*: Optimization takes too long, high memory usage

*Solutions*:
- **Reduce candidates**: ``num_candidates=3`` for faster iterations
- **Disable expensive features**: ``enable_adaptive_scaling=False``, ``enable_shuffling=False``
- **Use stricter early stopping**: ``error_tolerance=1e-4``, ``enable_early_convergence=True``
- **Optimize caching**: Check cache hit rate with ``get_cache_statistics()``
- **Pre-load models**: Use ``preloaded_model`` for epsilon properties

.. code-block:: python

    # For faster optimization
    optimizer = goose.SequenceOptimizer(
        target_length=50,
        max_iterations=1000,
        num_candidates=3,
        enable_adaptive_scaling=False,
        enable_shuffling=False,
        update_interval=100,  # Less frequent progress updates
        verbose=False
    )

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
- **Use emergency diversity**: The optimizer does this automatically
- **Check for impossible targets**: Some property combinations may be unachievable

**Epsilon Property Optimization Issues**

*Symptoms*: Epsilon-based properties don't optimize well

*Solutions*:
- **Pre-load models**: Use ``preloaded_model`` to avoid reloading
- **Adjust window sizes**: Try different ``window_size`` values (10-20)
- **Use appropriate scales**: Epsilon properties often need higher weights
- **Enable matrix resizing**: ``allow_matrix_resizing=True`` for matrix properties

.. code-block:: python

    # Optimized epsilon property usage
    import finches
    model = finches.Load('mpipi')
    
    optimizer.add_property(
        goose.MeanSelfEpsilon,
        target_value=-2.0,
        weight=5.0,  # Higher weight for epsilon properties
        preloaded_model=model,
        tolerance=0.1
    )

**Debugging and Monitoring**

*Enable detailed diagnostics*:

.. code-block:: python

    optimizer = goose.SequenceOptimizer(
        target_length=100,
        verbose=True,
        debugging=True,  # Enable detailed logging
        update_interval=10  # Frequent progress updates
    )
    
    # Monitor convergence
    convergence_info = optimizer.get_convergence_info()
    print(f"Current error: {convergence_info['current_error']}")
    print(f"Error trend: {convergence_info['error_trend']}")
    
    # Check cache performance
    cache_stats = optimizer.get_cache_statistics()
    print(f"Cache hit rate: {cache_stats['hit_rate']:.2%}")

**Common Property Target Ranges**

.. code-block:: text

    Property                  Typical Range    Notes
    ----------------------   --------------   -------------------------
    Hydrophobicity           0.0 - 6.6        Higher = more hydrophobic
    FCR                      0.0 - 1.0        Fraction charged residues
    NCPR                     -1.0 - 1.0       Net charge per residue
    Kappa                    0.0 - 1.0        Charge patterning
    FractionDisorder         0.0 - 1.0        Use disorder_cutoff=0.5
    RadiusOfGyration         Variable         Depends on sequence length
    EndToEndDistance         Variable         Depends on sequence length
    Complexity               0.0 - ~4.3       SEG complexity score

**Performance Benchmarks**

.. code-block:: text

    Sequence Length    Simple Properties    Complex Properties    Epsilon Properties
    ---------------    -----------------    ------------------    ------------------
    50 residues        < 30 seconds         1-3 minutes          2-5 minutes
    100 residues       1-2 minutes          5-10 minutes         10-20 minutes  
    200 residues       5-10 minutes         20-30 minutes        30-60 minutes
    500+ residues      20+ minutes          60+ minutes          120+ minutes

*Note: Times are approximate and depend on target complexity, number of properties, and hardware.*

Glossary
--------

**Optimization Terms**
- **adaptive scaling**: Dynamic adjustment of property weights based on convergence patterns
- **constraint type**: How property errors are calculated (EXACT, MINIMUM, MAXIMUM)
- **convergence**: When optimization error stabilizes and stops improving significantly
- **property**: A biophysical or biochemical feature to optimize (e.g., hydrophobicity)
- **stagnation**: When optimization makes no progress for extended periods
- **tolerance**: Per-property error threshold below which the property is considered satisfied
- **weight**: The relative importance of a property in the total optimization objective

**Sequence Modification Terms**
- **candidate**: A trial sequence generated during each optimization iteration
- **mutation**: Random amino acid substitutions made to explore sequence space
- **shuffling**: Rearranging amino acids within sequence regions to escape local minima
- **window_size**: The length of sequence segments for local shuffling operations

**Property-Specific Terms**
- **disorder cutoff**: Threshold for classifying residues as disordered (default: 0.5)
- **epsilon**: Interaction energy parameters from FINCHES for protein-protein interactions
- **FCR**: Fraction of Charged Residues (charged amino acids / total amino acids)
- **intermap**: Matrix representation of pairwise interactions between sequence positions
- **NCPR**: Net Charge Per Residue (positive charges - negative charges) / total amino acids
- **preloaded model**: Pre-initialized interaction models to improve performance

**Algorithm Parameters**
- **cache hit rate**: Percentage of property calculations avoided through caching
- **convergence patience**: Number of iterations to wait after convergence before stopping
- **convergence window**: Number of recent iterations examined for convergence detection
- **error tolerance**: Global error threshold for early optimization termination
- **improvement history**: Recent optimization progress tracking for each property
- **normalization factors**: Scaling values calculated from initial sequence properties
- **shuffle frequency**: How often sequence shuffling is performed (every N iterations)
- **stagnation threshold**: Number of iterations without improvement before stagnation response

Examples and Demo Notebooks
--------------------------

GOOSE includes comprehensive demo notebooks showcasing advanced ``SequenceOptimizer`` usage:

**Available Demos:**
- **Basic optimization**: Single and multi-property examples
- **Custom properties**: Creating and implementing user-defined properties  
- **Interaction optimization**: Epsilon-based properties and matrix manipulations
- **Performance optimization**: Efficient settings for different use cases
- **Advanced features**: Adaptive scaling, convergence detection, and troubleshooting

**Demo Location:**
Check the ``goose/demos/`` directory for Jupyter notebooks with detailed examples and explanations.

API Reference
-------------

**Core Classes:**
- ``goose.SequenceOptimizer``: Main optimization engine
- ``goose.backend.optimizer_properties.ProteinProperty``: Base class for properties
- ``goose.backend.optimizer_properties.ConstraintType``: Constraint type enumeration

**Key Methods:**
- ``SequenceOptimizer.add_property()``: Add properties to optimize
- ``SequenceOptimizer.set_initial_sequence()``: Set starting sequence
- ``SequenceOptimizer.configure_convergence()``: Configure convergence detection
- ``SequenceOptimizer.configure_error_tolerance()``: Configure error tolerance
- ``SequenceOptimizer.run()``: Execute optimization
- ``SequenceOptimizer.get_convergence_info()``: Get convergence diagnostics
- ``SequenceOptimizer.get_cache_statistics()``: Get performance statistics

See Also
--------

- :doc:`api`
- :doc:`getting_started`
- :doc:`sequence_generation`
- :doc:`variant_generation`
- :doc:`sequence_analysis`

For complete API documentation, see ``goose/optimize.py`` and ``goose/backend/optimizer_properties.py``.

For implementation examples and advanced usage patterns, explore the demo notebooks in ``goose/demos/``.
