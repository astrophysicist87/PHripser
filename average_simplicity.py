import numpy as np
import sys

outfilename = sys.argv[1]
data = np.array([np.loadtxt(file, usecols=(1,2,3,4,5,6)) for file in sys.argv[2:]])

#print(data.shape)
data = np.mean(data, axis=0)

np.savetxt(outfilename, data)