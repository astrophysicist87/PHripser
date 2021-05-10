import numpy as np
from scipy import interpolate

def get_interpolation(x, y, xnew):
    return (interpolate.interp1d(x, y, kind='cubic'))( xnew )

path = "C:/Users/Christopher Plumberg/Desktop/Research/UIUC/PHripser/"
multiplicities = [50, 75, 100]

for multiplicity in multiplicities:
    filename = path + "persistence_pair_statistics_n" + str(multiplicity) + ".dat"
    outfilename = path + "interpolated_persistence_pair_statistics_n" + str(multiplicity) + ".dat"

    [woBEaverageFiltration, woBEaverage, woBEstddev, \
     wBEaverageFiltration, wBEaverage, wBEstddev]    \
        = ((np.loadtxt(filename))[1:,0:6]).T

    overlapMin = np.amax([np.amin(woBEaverageFiltration), np.amin(wBEaverageFiltration)])
    overlapMax = np.amin([np.amax(woBEaverageFiltration), np.amax(wBEaverageFiltration)])

    n=30000
    filtrationPoints = np.linspace( overlapMin, overlapMax, n)
    results = np.c_[
                filtrationPoints,
                get_interpolation(woBEaverageFiltration, woBEaverage, filtrationPoints),
                get_interpolation(woBEaverageFiltration, woBEstddev,  filtrationPoints),
                get_interpolation(wBEaverageFiltration,  wBEaverage,  filtrationPoints),
                get_interpolation(wBEaverageFiltration,  wBEstddev,   filtrationPoints)    
                ]

    np.savetxt(outfilename, results, fmt='%f')
    print('Saved to', outfilename)

