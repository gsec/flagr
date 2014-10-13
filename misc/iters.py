import time, itertools, math, sys
import numpy as np
import matplotlib.pyplot as plt
def iters(mode='DIAGONAL'):
    #if mode == 'DIAGONAL':
    types = itertools.combinations_with_replacement([1,-1,0],3)
    print("Types:",types)
    perms = []
    #se = set()
    for i in types:
        t = set(itertools.permutations(i))
        while t:
            perms.append(t.pop())
        
    return np.pad(perms[:-1],(0,1), mode = 'constant')

nn = iters()
    
