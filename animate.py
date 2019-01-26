#!/usr/bin/python
import sys
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt 
from matplotlib import animation

def main():

    filename = sys.argv[1]
    N_values_int = [2**j for j in range(12)]
    n   =   int(sys.argv[2])   
    N_values_str = []
    for i in range(len(N_values_int)):
        N_values_str.append(str(N_values_int[i]))
    t_i= float(raw_input('Set initial time: '))
    tf = float(raw_input('Set time between frames(ms): '))
    t_i = int(t_i*100)
    split_filename = filename.split('_')
    N = 0 
    for i in range(len(N_values_str)):
        if N_values_str[i] in split_filename:
            for j in range(len(split_filename)):
                if split_filename[j] == N_values_str[i]:
                    N = N_values_int[i]

    arch = open(filename, 'r')
    ar = np.loadtxt(arch, usecols=(0, 1))
    ar=ar.transpose()
    L = N/3 - 1
    j = int(L)

    fig = plt.figure()
    ax = plt.axes(xlim=(0,j+0.1*j),ylim=(1e-15,0.2))

    plt.grid(which = 'both',axis = 'both', color = '0.01', linestyle = '--', linewidth = 0.1)
    plt.yscale('log')
    plt.xlabel('k')
    plt.ylabel('$E_d(k)$')
    plt.title('Spectrum(t)  s \n AHOe')
    line, = ax.plot([],[],lw=2)
    

    plt.plot(ar[0,2999*j:3000*j], ar[1,(n*i+t_i)*j:(n*i+t_i+1)*j],'r')
    def init():

        line.set_data([],[]) 
        return line,

    def animate(i):  
        f= float(n*i+t_i)
        t = f/100
        plt.title('Spectrum t = %.2f' %t)

        x = ar[0,(n*i+t_i)*j:(n*i+t_i+1)*j]
        x = x.clip(1e-33)
        y = ar[1,(n*i+t_i)*j:(n*i+t_i+1)*j]    
        line.set_data(x,y)
        line.set_color('blue')
        return line,      
 
    anim=animation.FuncAnimation(fig,animate,init_func=init,frames = 2000, interval =tf)
    plt.show()

if __name__ == '__main__':
    main()
