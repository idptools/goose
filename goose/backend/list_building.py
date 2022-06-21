# code for builging weighted lists


def gen_new_list(beta_seqs):
    ''' 
    makes new weighted list for seq generation. old, poorly written, 
    very functional code.
    '''
    
    all_amino_acids = ''
    for i in beta_seqs:
        all_amino_acids += i

    final_list = []
    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    # now iterate through aas
    for aa in amino_acids:
        #count number times aa appeared
        cur_count = all_amino_acids.count(aa)
        total_aa = len(all_amino_acids)
        # now add fractional value to final list
        if cur_count != 0:
            cur_fraction = cur_count / total_aa
            cur_number = round(cur_fraction * 5000)
            if cur_number == 0:
                cur_number=1
            else:
                cur_number=cur_number
            for num_aa in range(0, cur_number):
                final_list.append(aa)
        # figure out fraction of aas for each aa
        cur_fraction = all_amino_acids.count(aa)

    # printing on purpose tbh
    print(final_list)