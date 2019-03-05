#!/usr/bin/python
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import skew
from scipy.stats import entropy
import matplotlib.patches as mpatches
import os
from cycler import cycler
def main():



    filename = sys.argv[1]
    (head,tail) = os.path.split(filename)
    split_filename = tail.split("_")
    forcing =split_filename[2]
    (method,HD_MHD) = split_filename[3].split(".")
    lat_size = int(split_filename[4])
    visc = float(split_filename[5])
    arch = open(filename,'r')
    a = np.loadtxt(arch, usecols=(1))
    L = len(a)
    
    print  'length of data = ', L-49
    a = a[49:L]
    lyapunov = np.mean(a)
    std_dev = np.std(a)
    hist, bin_edge= np.histogram(a, bins = 'auto', density = 'True')
    length_bin_edge = len(bin_edge)-1
    bin_centres = np.array([])
    for i in range(0,length_bin_edge):
        c = (1/2)*(bin_edge[i+1]-bin_edge[i])
        bin_centres = np.append(bin_centres,c)
    
   # print 'lyapunov = ', lyapunov
   # print 'sd = ', std_dev
   # print 'rel_sd = ', std_dev/lyapunov
    N=len(a)

    

    def f(t):
        return (1/np.sqrt(2*np.pi*std_dev**2))*np.exp(-((t-lyapunov)**2)/(2*(std_dev**2)))
   #Tests

    sk = skew(a)
    ku= stats.kurtosis(a)
   # print 'skewness = ', sk
    Q = f(bin_centres)
    S = entropy(hist,Q)
    W, p_value = stats.shapiro(a)
    SSKK, p_value2 = stats.normaltest(a)

    SSKK2, p_value4 = stats.mstats.normaltest(a)
    DD, p_value3 = stats.kstest(a,'norm')
   # print 'W_shapiro = ', W
   # print 'p_value_shapiro= ', p_value 
   # print 'Kullback-Leibler Entropy = ',S 
   # print 'p_value_K-S= ', p_value3
   # print 'p_value_DAgostino= ', p_value2

    #print 'p_value_DAgostino2= ', p_value4
    minimum = np.amin(a)
    maximum = np.amax(a)
    stepp = (maximum - minimum)/1000

    params = {'legend.fontsize':'x-large','axes.labelsize':'xx-large','xtick.labelsize':'x-large','ytick.labelsize':'x-large'}
    plt.rcParams.update(params)
    plt.rc('text',usetex=True)
    plt.rc('font',family='serif')
    plt.rcParams['axes.prop_cycle'] = cycler(color = 'krgcmyb')
   #plt.style.use('seaborn-muted')

    fig, ax = plt.subplots()
    b = np.arange(minimum,maximum,stepp) 
    
    
    gaussian = ax.plot(b,f(b),'k-', label = 'gaussian')
    
    histogram = ax.hist(a,bins='auto',density = 'True',label= 'histogram',color = '#bc5c47')
    ax.set(xlabel = '$\\tilde{\lambda}$', ylabel = 'Probability density')# title = ' %s %s \n $\\nu=%.3f$  $N^3$' %(method,forcing,visc,lat_size))
   
    #first_legend = ax.legend(loc=1)
   # axi = plt.gca().add_artist(first_legend)
    plt.text(0.78,0.75,'$n = %i$ ' %(N), horizontalalignment='left', verticalalignment='center', fontsize = 'x-large',transform=ax.transAxes, bbox=dict(edgecolor= 'none',facecolor='white'))
    plt.legend(frameon = False,loc = 1)
    plt.show()


if __name__ == '__main__':
    main()
        
