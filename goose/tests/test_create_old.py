import pytest
from goose import create
from sparrow import Protein
from goose import goose_exceptions
import metapredict as meta
import random


def test_create_lengths():

    with pytest.raises(goose_exceptions.GooseInputError):
        seq = create.sequence(1)


    with pytest.raises(goose_exceptions.GooseInputError):
        seq = create.sequence(9)
    
    assert len(create.sequence(10)) == 10
    assert len(create.sequence(100)) == 100
    assert len(create.sequence(1000)) == 1000

    pool = set([])
    doubles = 0
    for i in range(0,100):
        s = create.sequence(20)
        if s in pool:
            doubles = doubles + 1
        else:
            pool.add(s)

        if doubles > 1:
            raise Exception('Generated more than 1 duplicate sequence of 100 random seqs')
        

    # check sequences are disordered
    for i in range(0,20):
        length = random.randint(20,200)
        seq = create.sequence(length)
        if meta.percent_disorder(seq) < 90:
            raise Exception(f'Generated a disordered sequence that is less than 90% disorder\n:{seq}')
            
        
        
        
