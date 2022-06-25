# various functions not integral to 
# folded region generation but were used in 
# optimizing it (to some extent...)

def gen_seq_list_beta(num_seqs=10000, seq_len=14, cutoff=0.7):
    successful_seqs = []
    for numseqs in range(0, num_seqs):
        curseq = gen_beta_starter(seq_len)
        if check_beta_strand(curseq, cutoff=cutoff):
            successful_seqs.append(curseq)
    print(successful_seqs)
    print(len(successful_seqs))



def gen_seq_list_helix(num_seqs=10000, seq_len=100, cutoff=0.7):
    successful_seqs = []
    for numseqs in range(0, num_seqs):
        curseq = gen_helix_starter(seq_len)
        if check_helicity(curseq, cutoff=cutoff):
            successful_seqs.append(curseq)
    print(successful_seqs)
    print(len(successful_seqs))


def gen_seq_list_coil(num_seqs=10000, seq_len=7, cutoff=0.7):
    successful_seqs = []
    for numseqs in range(0, num_seqs):
        curseq = gen_coil_starter(seq_len)
        if check_coil(curseq, cutoff=cutoff):
            successful_seqs.append(curseq)
    print(successful_seqs)
    print(len(successful_seqs))folded_region_backend_functions.py