# Test for consecutive indicies
from operator import itemgetter
from itertools import groupby
import  numpy as np

data = np.array([2, 3, 4, 5, 12, 13, 14, 15, 16, 17, 19])

startInd = np.array([], dtype=int)
endInd = np.array([], dtype=int)

for k, g in groupby(enumerate(data), lambda i: i[0]-i[1]):
    ind = list(map(itemgetter(1), g))
    startInd = np.append(startInd, ind[0])
    endInd = np.append(endInd, ind[-1])

for sI, eI in zip(startInd, endInd):
    print(data[sI:eI])