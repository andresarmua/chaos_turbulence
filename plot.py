#!/usr/bin/python
import sys
import numpy as np
import math
from scipy import stats
import matplotlib.pyplot as plt 
from matplotlib import pylab

def main():
    filename = sys.argv[1]

    arch = open(filename, 'r')
    ar = np.loadtxt(arch, usecols=(0, 1))
    ar=ar.transpose()

    x = ar[0]

    y = ar[1]

        
    
    pylab.plot(x,y)
    pylab.show()

if __name__ == '__main__':
    main()
