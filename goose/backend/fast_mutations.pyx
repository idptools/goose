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