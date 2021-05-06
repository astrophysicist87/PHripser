import numpy as np
import sys

def get_filenames(file):
    with open(file, "r") as a_file:
        for line in a_file:
            yield line.strip()

if len(sys.argv) > 3:
    outfilename = sys.argv[1]
    data = np.array([np.loadtxt(file, usecols=(1,2,3,4,5,6)) for file in sys.argv[2:]])
    #data = np.array([np.loadtxt(file) for file in sys.argv[2:]])
    #print(data.shape)
    data = np.std(data, axis=0)
    np.savetxt(outfilename, data)
else:
    outfilename = sys.argv[1]
    data = np.array([np.loadtxt(file, usecols=(1,2,3,4,5,6)) for file in get_filenames(sys.argv[2])])
    data = np.std(data, axis=0)
    np.savetxt(outfilename, data)

