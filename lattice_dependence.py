#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from os import path
from scipy import stats
from matplotlib import pylab
import sys


def main():
    
    print("--------------------------------------------------------\nThis program measures all Reynolds Numbers, T_0's, tau's and lattice sizes from .stats file and calculates lambda and error_lambda from .ftle files and plots error_lambda*tau against Re as well as lambda*T_0 vs Re \n --------------------------------------------------------\n\n")

    # Asks directory in which files are located
    # dir_name = raw_input("Directory where *stats and *ftle files are located: ")
    dir_name = sys.argv[1]
    #Defines arrays with all needed values 

    
    
    Reynolds_Numbers    = np.array([])
    Reynolds_Errors     = np.array([])

    Klmgrv_Microtimes   = np.array([])
    Klmgrv_Micro_Errors = np.array([])

    T_0s                = np.array([])
    T_0_Errors          = np.array([]) 
    
    Lyapunov_Exponents  = np.array([])
    Lyapunov_Errors     = np.array([])
    
    Lattice_Sizes       = np.array([])
    
    L                   = np.array([])
    L_error             = np.array([])
    U                   = np.array([])
    U_error             = np.array([])

    seeds               = np.array([])

    markers_seed        = np.array([])
    markers_lattice     = np.array([])
    leg                 = np.array([])
#---------------------------------stats files process-------------------------------------





    # Creates path and measures number of files

    Path = path.join(dir_name,'*.stats')    # Creates path
    stats_files = glob(Path)
    stats_files.sort()                      # Sorts files by numerical/alph order
    no_files = len(stats_files)             # Gets how many files are processes


    # Now data from the stats files is processesed and stored into arrays

    for i in range(no_files):
        
       
        # Read values from file name
        (head,tails) = path.split(stats_files[i])
        
        N_values = [2**j for j in range(12)]
        split_filename = tails.split("_")
        forcing = split_filename[2]
        method = split_filename[3]
        lat_size = int(split_filename[5])
        if lat_size not in N_values:
                exit()
        seed = int(list(split_filename[4])[1])
        (visc1,visc2,exten) = split_filename[6].split(".")
        visc = float(visc1 + "." + visc2)


        # Open file and read
        filename = open(stats_files[i],'r')
        
        
        
        ar = np.loadtxt(filename, usecols=(2,3,5,6))
        ar = ar.transpose()
        

        t_f = np.size(ar,1) - 1
        
        
    
        
        
        # Consider only steady state data
        
        diss   = ar[0,1500:t_f]  
        R_L    = ar[1,1500:t_f]
        u      = ar[2,1500:t_f]
        l      = ar[3,1500:t_f]
        
        # Calculates all necessary parameters
        
        
        epsilon         = np.mean(diss)
        epsilon_error   = np.std(diss)

        seeds           = np.append(seeds,seed)
        fmts    = ['*','v','s','o','D','X','p','P','<','^','H','>','8']
        for j in range(12):
            if 2**j == lat_size:
                k = j

        markers_seed    = np.append(markers_seed, fmts[seed])
        markers_lattice = np.append(markers_lattice, fmts[k])
        leg             = np.append(leg,N_values[k])
        L               = np.append(L,np.mean(l)) 
        L_error         = np.append(L_error, np.std(l))
        U               = np.append(U,np.mean(u))
        U_error         = np.append(U_error,np.std(u))

        Lattice_Sizes       = np.append(Lattice_Sizes, lat_size)
        Reynolds_Numbers    = np.append(Reynolds_Numbers,np.mean(R_L))    
        Reynolds_Errors     = np.append(Reynolds_Errors,np.std(R_L))
        
        

        Klmgrv_Microtimes   = np.append(Klmgrv_Microtimes,np.sqrt(visc/epsilon))
        Klmgrv_Micro_Errors = np.append(Klmgrv_Micro_Errors,(np.sqrt(visc)/2*(epsilon)**(3/2))*epsilon_error) 
        
        T_0s                = np.append(T_0s, L[i]/U[i])
        T_0_Errors          = np.append(T_0_Errors, (L[i]*U_error[i]/U[i]**2)+(L_error[i]/U[i]))
        




#---------------------------------ftle files process--------------------------------------
        




    # Creates path and measures number of files

    Path = path.join(dir_name,'*.ftle')    # Creates path
    ftle_files = glob(Path)
    ftle_files.sort()                      # Sorts files by numerical/alph order
    no_files = len(ftle_files)             # Gets how many files are processes

    
    


    # Now data from ftle's is processed and stored into arrays


    for i in range(no_files):
        
        # Open file and read
        filename = open(ftle_files[i],'r')
        ar = np.loadtxt(filename, usecols=(0,1))
        ar = ar.transpose()
        
        t_f = np.size(ar,1) - 1

        t_i = 300              #  ~ equivalent to 15s 

        Lyapunov_sample = ar[1,t_i:t_f]
        Lyapunov_Exponents  = np.append(Lyapunov_Exponents,np.mean(Lyapunov_sample))
        Lyapunov_Errors     = np.append(Lyapunov_Errors,np.std(Lyapunov_sample))




    
    params = {'legend.fontsize':'x-large','axes.labelsize':'xx-large','xtick.labelsize':'x-large','ytick.labelsize':'x-large'}
    
    plt.rcParams.update(params)
    plt.rc('text',usetex=True)
    plt.rc('font',family='serif')

    plt.style.use('seaborn-dark-palette')


    #plt.plot(Lyapunov_Exponents, Lyapunov_Errors, 'o', label = 'ErrorDependence')
    #plt.xlabel('<\lambda>')
    #plt.ylabel('\Delta \lambda')
    #plt.xlabel('$\lambda$')
    #plt.ylabel('$\Delta \lambda$')
    #plt.show()

    for x,y,er,m in zip(Lattice_Sizes,Lyapunov_Exponents, Lyapunov_Errors, markers_seed):
        plt.errorbar(x,y,yerr = er,label = 'Error Dependence',color = '#bc5c47',marker = m, ecolor = '0.6',capsize = 3, elinewidth = 1, capthick = 1)
    plt.xscale('log',basex=2)
    plt.grid(which ='both', axis = 'both', linewidth = 0.2, linestyle='--')
    plt.xlabel('$N$')
    plt.ylabel('$\lambda$')
    plt.show()

   # plt.plot(Lattice_Sizes,Reynolds_Numbers,'.',color = '#bc5c47',label = 'Error Dependence')
   # plt.xscale('log',basex = 2)
   #   #plt.yscale('log')
   # plt.xlabel('$N$')
   # plt.ylabel('$Re$')
   # plt.show()
    
    
    trigger = 14*[False]
    for x,y,m,l in zip(Reynolds_Numbers,Lyapunov_Exponents,markers_lattice,leg):
        h = int(np.log2(l))
        if trigger[h] == False:
            plt.plot(x,y,m,color = '#bc5c47',label='$N = %i$' %l)
            trigger[h] = True
        else:
            plt.plot(x,y,m,color = '#bc5c47')
    #plt.xscale('log')
    plt.grid(which ='both', axis = 'both', linewidth = 0.2, linestyle='--')
    plt.xlabel('$Re$')
    plt.ylabel('$\lambda$')
    plt.legend()
    plt.show()

    
#    for x,y,m in zip(Lattice_Sizes,Lyapunov_Errors*Klmgrv_Microtimes,markers_seed):
 #       plt.plot(x,y, m,label = 'Error Dependence',color = '#bc5c47')

 #   plt.xscale('log',basex=2)
    #plt.yscale('log')
 #   plt.xlabel('N')
#    plt.ylabel('$\Delta \lambda \\tau$')
 #   plt.grid(which ='both', axis = 'both', linewidth = 0.2, linestyle='--')
 #   plt.title('Error dependence on Lattice')
 #   plt.show()


    for x,y,m in zip(Lattice_Sizes,Lyapunov_Errors*T_0s,markers_seed):
        plt.plot(x,y, m,label = 'Error Dependence', color = '#bc5c47')

    #plt.errorbar(Lattice_Sizes,Lyapunov_Errors*T_0s, yerr = Lyapunov_Errors*T_0_Errors, fmt= markers, ecolor = 'g',capsize =3, elinewidth = 1, capthick = 1,label = 'Error Dependence')
    #plt.plot(Lattice_Sizes,Lyapunov_Errors*T_0s, 'o', label = 'Error Dependence')
    plt.xscale('log',basex=2)
    #plt.yscale('log')
    plt.xlabel('$N$')
    plt.ylabel('$\sigma_{\lambda}\, T_0$')
    plt.grid(which ='both', axis = 'both', linewidth = 0.2, linestyle='--')
    #plt.title('Error dependence on Lattice')
    plt.show()

   # for x,y,m in zip(Lattice_Sizes,Lyapunov_Errors/Lyapunov_Exponents,markers_seed):
   #     plt.plot(x,y,m,color= '#bc5c47',label = 'Error Dependence')

    #plt.plot(Lattice_Sizes,Lyapunov_Errors/Lyapunov_Exponents, 'o', label = 'Error Dependence')
  #  plt.xscale('log',basex = 2)
    #plt.yscale('log')
  #  plt.xlabel('$N$')
  #  plt.ylabel('$\\frac{\sigma_{\lambda}}{\lambda}$')
  #  plt.grid(which ='both', axis = 'both', linewidth = 0.2, linestyle='--')
    #plt.title('Relative error vs Lattice Size')
   # plt.show()





if __name__ == '__main__':
    main()
