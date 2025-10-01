# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

from libc.stdlib cimport malloc, free, rand, srand, RAND_MAX
from libc.string cimport strcpy
from libc.time cimport time


# Seed C random number generator
srand(time(NULL))

cdef const unsigned char* AMINO_ACIDS = b"ACDEFGHIKLMNPQRSTVWY"
cdef int AMINO_ACIDS_LEN = 20


cpdef str cython_shuffle_sequence(str sequence, int window_size=10, float global_shuffle_probability=0.3):
    cdef int seq_len, i, j, start
    cdef char* seq_buffer
    cdef bytes py_bytes
    cdef char temp

    seq_len = len(sequence)
    if seq_len <= 1:
        return sequence

    seq_buffer = <char*>malloc(seq_len + 1)
    if not seq_buffer:
        raise MemoryError()

    # The `try...finally` block here ensures our C buffer is always freed 
    # even if errors occur we could also rewrite this with vecs which would be preferrable
    # but then we require clang or some other c++ compiler and lets keep this simple for now
    try:
        py_bytes = sequence.encode('utf-8')
        strcpy(seq_buffer, py_bytes)

        if (<float>rand() / RAND_MAX) < global_shuffle_probability:
            for i in range(seq_len - 1, 0, -1):
                j = rand() % (i + 1)
                temp = seq_buffer[i]
                seq_buffer[i] = seq_buffer[j]
                seq_buffer[j] = temp
        else:
            if seq_len > window_size:
                start = rand() % (seq_len - window_size + 1)
                for i in range(window_size - 1, 0, -1):
                    j = rand() % (i + 1)
                    temp = seq_buffer[start + i]
                    seq_buffer[start + i] = seq_buffer[start + j]
                    seq_buffer[start + j] = temp

        return seq_buffer[:seq_len].decode('utf-8')

    finally:
        # don't leak memory like a scrub with malloc
        free(seq_buffer)


cpdef str cython_make_point_mutation(str sequence, int position=-1):
    cdef int seq_len, pos
    cdef char* seq_buffer
    cdef bytes py_bytes
    cdef char old_aa, new_aa

    seq_len = len(sequence)
    if seq_len == 0:
        return sequence

    seq_buffer = <char*>malloc(seq_len + 1)
    if not seq_buffer:
        raise MemoryError()

    try:
        py_bytes = sequence.encode('utf-8')
        strcpy(seq_buffer, py_bytes)

        pos = position
        if pos == -1:
            pos = rand() % seq_len

        old_aa = seq_buffer[pos]


        while True:
            new_aa = AMINO_ACIDS[rand() % AMINO_ACIDS_LEN]
            if new_aa != old_aa:
                break
        
        seq_buffer[pos] = new_aa

        return seq_buffer[:seq_len].decode('utf-8')

    finally:
        free(seq_buffer)

cpdef str cython_make_multi_mutations(str sequence, int num_mutations):
    """
    Apply multiple point mutations to a sequence efficiently.
    
    Parameters
    ----------
    sequence : str
        The input amino acid sequence
    num_mutations : int
        Number of mutations to apply (must be <= sequence length)
    
    Returns
    -------
    str
        Mutated sequence with num_mutations point mutations at random positions
    
    Notes
    -----
    - Ensures no position is mutated twice (samples without replacement)
    - Ensures each mutation changes the amino acid (no silent mutations)
    - Much faster than calling cython_make_point_mutation in a loop
    """
    cdef int seq_len, i, j, pos, mutations_left
    cdef char* seq_buffer
    cdef char* positions_used
    cdef bytes py_bytes
    cdef char old_aa, new_aa
    cdef int attempts

    seq_len = len(sequence)
    if seq_len == 0 or num_mutations <= 0:
        return sequence
    
    # Clamp num_mutations to sequence length
    if num_mutations > seq_len:
        num_mutations = seq_len

    seq_buffer = <char*>malloc(seq_len + 1)
    if not seq_buffer:
        raise MemoryError()

    # Track which positions have been mutated (0 = unused, 1 = used)
    positions_used = <char*>malloc(seq_len)
    if not positions_used:
        free(seq_buffer)
        raise MemoryError()

    try:
        # Initialize
        py_bytes = sequence.encode('utf-8')
        strcpy(seq_buffer, py_bytes)
        
        # Mark all positions as unused
        for i in range(seq_len):
            positions_used[i] = 0

        mutations_left = num_mutations

        # Apply mutations
        while mutations_left > 0:
            # Find unused position (rejection sampling)
            attempts = 0
            while attempts < seq_len * 2:  # Prevent infinite loop
                pos = rand() % seq_len
                if positions_used[pos] == 0:
                    break
                attempts += 1
            
            # If we couldn't find unused position, fall back to sequential selection
            if attempts >= seq_len * 2:
                for i in range(seq_len):
                    if positions_used[i] == 0:
                        pos = i
                        break
            
            # Mark position as used
            positions_used[pos] = 1
            
            # Get old amino acid
            old_aa = seq_buffer[pos]
            
            # Select new amino acid (different from old)
            while True:
                new_aa = AMINO_ACIDS[rand() % AMINO_ACIDS_LEN]
                if new_aa != old_aa:
                    break
            
            # Apply mutation
            seq_buffer[pos] = new_aa
            mutations_left -= 1

        return seq_buffer[:seq_len].decode('utf-8')

    finally:
        free(seq_buffer)
        free(positions_used)
