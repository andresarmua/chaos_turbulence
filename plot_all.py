#!/usr/bin/python
# This program plots Column 1 vs Column 2 of all files with a given extension in a given directory

import numpy as np
import matplotlib.pyplot as plt 
from glob import glob
from os import path

def main():
    print("This program plots column 1 vs. column 2 of all files in a directory, provided the extension")
    dir_name        = raw_input("Choose a directory(without initial /): ")      #asks for inputs
    extension       = raw_input("Choose an extension: ") 
    Column_1_label  = raw_input("column 1 label: ")
    Column_2_label  = raw_input("column 2 label: ") 
    t_i       = int(raw_input("choose initial time: "))
    
    Path = path.join(dir_name,'*.'+extension)                   #Creates entire path and 
    filename = glob(Path)                                       #sorts files numericallmpoy   
    filename.sort()
    no_files = len(filename)
    

    for i in range(no_files):                                   #plots all graphs
        arch = open(filename[i], 'r')
        ar = np.loadtxt(arch, usecols=(0, 1))
        ar=ar.transpose()
        

        lx = ar[0]
        ly = ar[1]

        x = lx[np.where(lx>t_i)]

        y = ly[np.where(lx>t_i)]
        

        plt.plot(x,y,'b.')
        plt.xlabel(Column_1_label)
        plt.ylabel(Column_2_label) 
        plt.title('$E$ vs. t for H0%i' %(i+1))
        plt.show()

    

if __name__ == '__main__':
    main()
