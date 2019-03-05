#!/usr/bin/python3

import numpy as np
import sys
import matplotlib.pyplot as plt
from glob import glob
from os import path
from scipy import stats
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
#--------This program measures Lambda vs steptime vs Re and Delta_Lambda vs steptime --------------- 



def main():
    
#    print("--------------------------------------------------------\nThis program measures the dependence of the Lyapunov exponents and its fluctuations respect to the change in steptime in the FTLE method.\n --------------------------------------------------------\n\n")
    
    
    
    # Asks directory in which files are located
    #dir_name = raw_input("Directory where *.stats and *.ftle files are located: ")             
    dir_name = sys.argv[1]
    # Defines arrays with all needed values

   # Reynolds_Numbers    = np.array([])      
   # Reynolds_Errors     = np.array([])

    #Klmgrv_Microtimes   = np.array([])
    #Klmgrv_Micro_Errors = np.array([])

    #T_0s                = np.array([])
    #T_0_Errors          = np.array([]) 
    
    Lyapunov_Exponents  = np.array([])
    Lyapunov_Errors     = np.array([])
    dt                  = np.array([])
    #Urms                = np.array([])
    #Urms_Errors         = np.array([])

    #L                   = np.array([])
    #L_Errors            = np.array([])

#---------------------------------stats files process-------------------------------------





    # Creates path and measures number of files

    #Path = path.join(dir_name,'*.stats')    # Creates path
    #stats_files = glob(Path)
    #stats_files.sort()                      # Sorts files by numerical/alph order
    #no_files = len(stats_files)             # Gets how many files are processes


    
    # Now data from the stats files is processesed and stored into arrays

    #for i in range(no_files):
        
        # Read values from file name
      #  (head,tails) = path.split(stats_files[i])
        
     #   N_values = [2**j for j in range(12)]
     #   split_filename = tails.split("_")
     #   dt = float(split_filename[1])
    #    method,lat = split_filename[2].split('.')
    #    lat_size = int(lat)
    #    if lat_size not in N_values:
    #            print('wrong file reading')
    #            exit()
    #    visc = float(split_filename[3])

     #   # Open file and read
     #   filename = open(stats_files[i],'r')
     #   ar = np.loadtxt(filename, usecols=(2,8,16))  # Ekin, diss, L
     #   ar = ar.transpose()
        

     #   t_f = np.size(ar,1) - 1
        
        
    
        
        
        # Consider only steady state data
        
     #   Ekin   = ar[0,2000:t_f]  
     #   diss   = ar[1,2000:t_f]
     #   L      = ar[2,2000:t_f]
        
        # Calculates all necessary parameters
        
     #   epsilon             = np.mean(diss)
     #   epsilon_error       = np.std(diss)
     #   Li                  = np.mean(L) 
     #   Li_error            = np.std(L)
     #   U                   = np.mean(np.sqrt(2*Ekin/3))
     #   U_error             = np.std(np.sqrt(2*Ekin/3))
#
#        R_L                = U*L/visc  
 #       Reynolds_Numbers    = np.append(Reynolds_Numbers,np.mean(R_L))    
 #       Reynolds_Errors     = np.append(Reynolds_Errors,np.std(R_L))
 #       
 #       
#
#        Klmgrv_Microtimes   = np.append(Klmgrv_Microtimes,np.sqrt(visc/epsilon))
#        Klmgrv_Micro_Errors = np.append(Klmgrv_Micro_Errors,(np.sqrt(visc)/2*(epsilon)**(3/2))*epsilon_error) 
#        Log_tau             = np.log10(Klmgrv_Microtimes)    
#        T_0s                = np.append(T_0s, Li/U)
#        T_0_Errors          = np.append(T_0_Errors, (Li*U_error/U**2)+(Li_error/U))
#        
#        Urms                = np.append(Urms, U)
#        Urms_Errors         = np.append(Urms_Errors,U_error)
#
#        L                   = np.append(L, Li)
#        L_Errors            = np.append(L_Errors,Li_error)
#
#
#---------------------------------ftle files process--------------------------------------
        



    # Creates path and measures number of files

    Path = path.join(dir_name,'*.ftle')    # Creates path
    ftle_files = glob(Path)
    ftle_files.sort()                      # Sorts files by numerical/alph order
    no_files = len(ftle_files)             # Gets how many files are processes

    


    # Now data from ftle's is processed and stored into arrays


    for f in ftle_files:
        (head,tails) = path.split(f)
        N_values = [2**j for j in range(12)]
        split_filename = tails.split("_")
        dt = np.append(dt,float(split_filename[1]))
        method,lat = split_filename[2].split('.')
        lat_size = int(lat)
        if lat_size not in N_values:
                print('wrong file reading')
                exit()
        visc = float(split_filename[3])

        # Open file and read
        filename = open(f,'r')
        ar = np.loadtxt(filename, usecols=(0,1))
        ar = ar.transpose()
        t_f = np.size(ar[0]) - 1
        trig = False
        t_i = 0
        for k in range(len(ar[0])):
            if int(ar[0,k]) == 15 and trig == False:
                t_i = k
                trig == True
                  

        Lyapunov_sample     = ar[1,t_i:t_f]
        Lyapunov_Exponents  = np.append(Lyapunov_Exponents,np.mean(Lyapunov_sample))
        Lyapunov_Errors     = np.append(Lyapunov_Errors,np.std(Lyapunov_sample))
        Log_L = np.log10(Lyapunov_Exponents)
    
#--------------------------------------- set error bars -----------------------


    ''' Error_Lyap_T        = Lyapunov_Errors*T_0s + T_0_Errors*Lyapunov_Exponents
    Error_DeltaL_tau    = Lyapunov_Errors * Klmgrv_Micro_Errors
    Error_DeltaL_T      = Lyapunov_Errors * T_0_Errors 
    Error_DeltaLyap     = (Lyapunov_Errors ** 2)/(Lyapunov_Exponents ** 2)
    Error_Log_L         = Lyapunov_Errors/Lyapunov_Exponents
    ErrorFabiola Gianotti_Log_tau       = Klmgrv_Micro_Errors/Klmgrv_Microtimes
#-------------------------  define fit functions -----------------------------

    # fit lambda*T_0 ~ Re^alpha

    Log_R  = np.log10(Reynolds_Numbers)
    Log_Ly = np.log10(Lyapunov_Exponents*T_0s)
    
    Log_E_Ly = np.log10(Lyapunov_Errors * T_0s)


    # WLS fit 
    
    Log_R = sm.add_constant(Log_R)  # adds constant b in y = mx + b
    mod_wls = sm.WLS(Log_Ly, Log_R, weights = ((Lyapunov_Exponents*T_0s)/Error_Lyap_T)**2)
    res_wls = mod_wls.fit()
    intercept, slope = res_wls.params
    intercept_e,std_err = res_wls.bse
    r_value = res_wls.rsquared    

    # linear regression on Delta_Lyap vs Lyap
    
    slope2, intercept2, r_value2, p_value2 , std_err2 = stats.linregress(Lyapunov_Exponents, Lyapunov_Errors)
    

    # linear regression on Delta_Lyap*T_0 vs Re 

    Log_Ly_Er = np.log10(Lyapunov_Errors*T_0s)



    # WLS fit 
    
    mod_wls3 = sm.WLS(Log_Ly_Er, Log_R)
    res_wls3 = mod_wls3.fit()
    intercept3, slope3 = res_wls3.params
    intercept_e3,std_err3 = res_wls3.bse
    r_value3 = res_wls3.rsquared    
    # linear regression on Delta_Re vs Re


    slope4, intercept4, r_value4, p_value4, std_err4 = stats.linregress(Reynolds_Numbers, Reynolds_Errors)


    # linear regression on Delta_T vs T


    slope5, intercept5, r_value5, p_value5, std_err5 = stats.linregress(T_0s,T_0_Errors)

    #linear regression on Delta Lyap /Lyap vs Delta_Re/Re


    slope6, intercept6, r_value6, p_value6, std_err6 = stats.linregress(Reynolds_Errors/Reynolds_Numbers,Error_Lyap_T/(Lyapunov_Exponents*T_0s))



    #linear regression on Delta Lyap_T vs Lyap_T 


    slope8, intercept8, r_value8, p_value8, std_err8 = stats.linregress(Lyapunov_Exponents*T_0s,Error_Lyap_T)


    #linear regression on Urms vs Delta_Urms


    slope9, intercept9, r_value9, p_value9, std_err9 = stats.linregress(Urms,Urms_Errors)



    #linear regression on L vs Delta_L


#    slope10, intercept10, r_value10, p_value10, std_err10 = stats.linregress(L,L_Errors)


    #linear regression on Delta_Lyap/ Lyap vs Delta_Re/Re

    slope11, intercept11, r_value11, p_value11, std_err11 = stats.linregress(Reynolds_Errors/Reynolds_Numbers,Lyapunov_Errors/Lyapunov_Exponents)
    
    #WLS fit on Lyap vs tau 

    
    Log_tau = sm.add_constant(Log_tau)    # adds constant y = mx + b
    mod_wls12 = sm.WLS(Log_L, Log_tau, weights = Error_Log_L)
    res_wls12 = mod_wls12.fit()
    intercept12, slope12 = res_wls12.params
    intercept_e12,std_err12 = res_wls12.bse
    r_value12 = res_wls12.rsquared    
    

#-------------------------------- Define fit functions -------------------------- 

    # Set fit function to be plot
    def f(t):
        return 10**(intercept)*t**(slope)

    low = np.amin(Reynolds_Numbers)
    maxi = np.amax(Reynolds_Numbers)
    values = np.arange(low,maxi,0.1)


    def g(t):
        return slope2*t + intercept2 

    low2 = np.amin(Lyapunov_Exponents)
    maxi2 = np.amax(Lyapunov_Exponents)
    values2 = np.arange(low2,maxi2,0.001)
    

    def h(t):
        return 10**(intercept3)*t**(slope3)
    
   
    def i(t):
       return slope4*t + intercept4
   

    def j(t):
        return slope5*t + intercept5
    
    low5 = np.amin(T_0s)
    maxi5 = np.amax(T_0s)
    values5 = np.arange(low5,maxi5,0.001)
    
    
    def k(t):
        return slope6*t + intercept6

    low6 = np.amin(Reynolds_Errors/Reynolds_Numbers)
    maxi6 = np.amax(Reynolds_Errors/Reynolds_Numbers)
    values6 = np.arange(low6,maxi6,0.001)

    def l(t):
        return 0.1 + 0.53*0.08  + np.log(t)*0.02


    def m(t):
        return slope8*t + intercept8
    
    low8 = np.amin(Lyapunov_Exponents*T_0s)
    maxi8 = np.amax(Lyapunov_Exponents*T_0s)
    values8 = np.arange(low8,maxi8,0.001)



    def n(t):
        return slope9*t + intercept9
    
    low9 = np.amin(Urms)
    maxi9 = np.amax(Urms)
    values9 = np.arange(low9,maxi9,0.001)

    def o(t):
        return slope10*t + intercept10


    low10 = np.amin(L)
    maxi10 = np.amax(L)
    values10 = np.arange(low10,maxi10,0.001)


    def p(t):
        return slope11*t + intercept11


    def q(t):
        return slope12*t + intercept12


    low12 = np.amin(Klmgrv_Microtimes)
    maxi12 = np.amax(Klmgrv_Microtimes)
    values12 = np.arange(low12,maxi12,0.001)


    '''
#--------------------------- plot everything --------------------
    
    
    params = {'legend.fontsize':'x-large','axes.labelsize':'xx-large','xtick.labelsize':'x-large','ytick.labelsize':'x-large'}
    plt.rcParams.update(params)
    plt.rc('text',usetex=True)
    plt.rc('font',family='serif')

    #plt.style.use('seaborn-dark-palette')
    
    fig,[ax1,ax2] =  plt.subplots(2,1)
    
    # plot Lyap vs dt

   
    ax1.errorbar(dt, Lyapunov_Exponents, yerr = Lyapunov_Errors, color = '#bc5c47', marker = '.',markerfacecolor = 'None',linestyle = 'None',label ='Data', ecolor = '0.6', capsize = 3, elinewidth = 1, capthick = 1)
    #plt.xscale('log')
    #plt.yscale('log')
    ax1.set_xlabel('$\Delta t$')
    ax1.set_ylabel('$\lambda$')
    ax1.grid(which='both', axis='both',linewidth = 0.2, linestyle = '--')
    #ax1.set_title('$\lambda \,vs\,\Delta t$')
    
    ax2.plot(dt, Lyapunov_Errors, marker = '.', markerfacecolor = 'None', color = '#bc5c47',linestyle = 'None')
    
    ax2.set_xlabel('$\Delta t$')
    ax2.set_ylabel('$\sigma_{\lambda}$')
    ax2.grid(which='both', axis='both',linewidth = 0.2, linestyle = '--')

    
    plt.savefig('/home/andres/Desktop/Lambda_and_error_vs_steptime.png',bbox_inches='tight')
    plt.show()
   
    '''
    # plot Lyap vs Re and linear fit
    data = plt.errorbar(Reynolds_Numbers, Lyapunov_Exponents*T_0s, xerr = Reynolds_Errors, yerr = Error_Lyap_T, fmt ='o',label ='Data', ecolor = 'g', capsize = 3, elinewidth = 1, capthick = 1)
    fit = plt.plot(values,f(values),label = 'Linear fit')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Re')
    plt.ylabel('$\lambda T_0$')
    plt.grid(which='both', axis='both',linewidth = 0.2, linestyle = '--')
    plt.title('$\lambda \, T_0 vs Re$')
    plt.text(0.7,0.2,'$\lambda \,  T_0 = D Re^{\\alpha}$ \n $\\alpha= %.2f \pm %.2f$ \n $r = %.2f$ \n D = %.3f+%.3f'%(slope,std_err,r_value,10**intercept,10**intercept*np.log(10)*intercept_e), horizontalalignment = 'center', verticalalignment = 'center',transform = ax.transAxes, bbox = dict(facecolor = 'white'))
    plt.legend(shadow = 'True',loc=2)
    plt.show()
    
    # plot Delta_Lyap vs Lyap and linear fit
    

    data2 = plt.plot(Lyapunov_Exponents,Lyapunov_Errors, 'o',label = 'Data')
    fit2 = plt.plot(values2,g(values2),label = 'Linear fit')
    plt.xlabel('$\lambda$')
    plt.ylabel('$\Delta \lambda$')
    plt.grid(which='both', axis='both',linewidth = 0.2, linestyle = '--')
    plt.title('$\Delta \lambda \, vs \,  \lambda$')
    plt.text(0.7,0.2,'$\Delta \lambda  \simeq c_2 \, \lambda$ \n $c_2= %.2f \pm %.2f$ \n $r = %.2f$'%(slope2,std_err2,r_value2), horizontalalignment = 'center', verticalalignment = 'center',transform = ax.transAxes, bbox = dict(facecolor = 'white'))
    plt.legend(shadow = 'True',loc=2)
    plt.show()

    # plot Delta_Lyap*T_0s vs Re and linear fit
    
    
    
    
    plt.errorbar(Reynolds_Numbers,Lyapunov_Errors*T_0s, xerr = Reynolds_Errors,yerr = Error_DeltaL_T, fmt= 'o', ecolor = 'g', capsize= 3,elinewidth = 1, capthick = 1, label = 'Data ')
    fit3 = plt.plot(values,h(values), label = 'Linear fit')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Re')
    plt.ylabel('$\Delta \lambda \, T_0$')
    plt.grid(which='both', axis='both',linewidth = 0.2, linestyle = '--')
    plt.text(0.7,0.15,'$\Delta \lambda T_0  \propto Re^{\gamma}$ \n$\\gamma=%.2f \pm %.2f$ \n $r = %.2f$'%(slope3,std_err3,r_value3), horizontalalignment = 'center', verticalalignment = 'center',transform = ax.transAxes,bbox= dict(facecolor = 'white'))
    plt.title('$\Delta \lambda \, T_0$ vs. Re')
    plt.legend(shadow = 'True',loc=2)
    plt.show()

    
    # plot Delta_Lyap*tau vs Re

    plt.errorbar(Reynolds_Numbers,Lyapunov_Errors*Klmgrv_Microtimes, xerr = Reynolds_Errors, yerr = Error_DeltaL_tau, fmt= 'o', ecolor = 'g', capsize = 3, elinewidth = 1, capthick = 1, label = '$\Delta \lambda \, \\tau$ vs. Re ')
    plt.xlabel('Re')
    plt.ylabel('$\Delta \lambda \, \\tau$')
    plt.grid(which='both', axis='both',linewidth = 0.2, linestyle = '--')
    plt.title('$\Delta \lambda  \, \\tau$ vs. Re')
    plt.show()




    
    # plot Delta_Lyap T/Lyap T vs Re
#    fig7, ax7 = plt.subplots()
#    ax7.plot(np.log(Reynolds_Numbers), Error_Lyap_T/(Lyapunov_Exponents*T_0s), 'o', label= 'Data')
#    ax7.plot(np.log(values), l(values), label = 'curve fit')
    #plt.xscale('log')
    #plt.yscale('log')
#    plt.xlabel('log(Re)')
#    plt.ylabel('$\\frac{\Delta (\lambda T_0)}{\lambda T_0}$')
#    plt.grid(which='both', axis='both',linewidth = 0.2, linestyle = '--')
#    plt.title('$\\frac{\Delta (\lambda T_0)}{\lambda T_0}$ vs. Re')
#    plt.show()

    # plot Delta(Lyap T) vs Lyap T
    
    plt.plot(Lyapunov_Exponents*T_0s,Error_Lyap_T,'o', label = 'Data ')
    fit8 = plt.plot(values8,m(values8), label = 'linear fit')
    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel('$\lambda T_0$')
    plt.ylabel('$\Delta (\lambda \, T_0)$')
    plt.grid(which='both', axis='both',linewidth = 0.2, linestyle = '--')
    plt.text(0.7,0.2,'$\Delta (\lambda T_0) = A \lambda T_0$ \n A= $%.2f \pm %.2f$ \n $r = %.2f$'%(slope8,std_err8,r_value8), horizontalalignment = 'center', verticalalignment = 'center',transform = ax.transAxes,bbox= dict(facecolor = 'white'))
    plt.title('$\Delta (\lambda \, T_0)$ vs. $\lambda T_0$')
    plt.legend(shadow = 'True',loc=2)
    plt.show()



    
    # plot Delta_U vs U 
#    fig9, ax9 = plt.subplots()
#    ax9.plot(Urms, Urms_Errors,'o', label = 'Data')
#    fit9 = ax9.plot(values9,n(values9),label = 'linear fit')
#    plt.xlabel('$U_{rms}$')
#    plt.ylabel('$\Delta U_{rms}$')
#    plt.grid(which='both', axis='both',linewidth = 0.2, linestyle = '--')
#    plt.title('$\Delta U_{rms} vs U_{rms}$')
#    plt.text(0.5,0.2,'$\Delta U_{rms}  = m U_{rms} + b$ \n $ m = %.2f \pm %.2f$ \n $r= %.2f$ \n $b= %.2f$'%(slope9,std_err9,r_value9,intercept9), horizontalalignment = 'center', verticalalignment = 'center',transform = ax.transAxes, bbox = dict(facecolor = 'white'))
#    plt.legend(shadow = 'true',loc=2)
#    plt.show()


    # plot Delta_L vs L 
#    fig10, ax10 = plt.subplots()
#    ax10.plot(L, L_Errors,'o', label = 'Data')
#    fit10 = ax10.plot(values10,o(values10),label = 'linear fit')
#    plt.xlabel('$L$')
#    plt.ylabel('$\Delta L$')
#    plt.grid(which='both', axis='both',linewidth = 0.2, linestyle = '--')
#    plt.title('$\Delta L vs L$')
#    plt.text(0.5,0.1,'$\Delta L  = m L + b$ \n $ m = %.2f \pm %.2f$ \n $r= %.2f$ \n $b= %.2f$'%(slope10,std_err10,r_value10,intercept10), horizontalalignment = 'center', verticalalignment = 'center',transform = ax.transAxes, bbox = dict(facecolor = 'white'))
#    plt.legend(shadow = 'true',loc=2)
#    plt.show()
    
    # plot Delta_T vs T 
#    fig5, ax5 = plt.subplots()
#    ax5.plot(T_0s, T_0_Errors,'o', label = 'Data')
#    fit5 = ax5.plot(values5,j(values5),label = 'linear fit')
#    plt.xlabel('$T_0$')
#    plt.ylabel('$\Delta T_0$')
#    plt.grid(which='both', axis='both',linewidth = 0.2, linestyle = '--')
#    plt.title('$\Delta T_0 vs T_0$')
#    plt.text(0.5,0.1,'$\Delta T_0  = m T_0 + b$ \n $ m = %.2f \pm %.2f$ \n $r= %.2f$ \n $b= %.2f$'%(slope5,std_err5,r_value5,intercept5), horizontalalignment = 'center', verticalalignment = 'center',transform = ax.transAxes, bbox = dict(facecolor = 'white'))
#    plt.legend(shadow = 'true',loc=2)
#    plt.show()
    
    # plot Delta_Re vs Re 
    #fig4, ax4 = plt.subplots()
    plt.plot(Reynolds_Numbers, Reynolds_Errors,'o', label = 'Data')
    fit4 = plt.plot(values,i(values),label = 'Linear fit')
    plt.xlabel('$Re$')
    plt.ylabel('$\Delta Re$')
    plt.grid(which='both', axis='both',linewidth = 0.2, linestyle = '--')
    plt.title('$\Delta Re \, vs \, Re$')
    plt.text(0.7,0.2,'$\Delta Re  = c_1 Re$ \n $ c_1 = %.3f \pm %.3f$ \n $r = %.2f$'%(slope4,std_err4,r_value4), horizontalalignment = 'center', verticalalignment = 'center',transform = ax.transAxes, bbox = dict(facecolor = 'white'))
    plt.legend(shadow = 'true',loc=2)
    plt.show()

    # plot Delta_L_T/L_T vs Delta_R/R
    plt.plot(Reynolds_Errors/Reynolds_Numbers,Error_Lyap_T/(Lyapunov_Exponents*T_0s), 'o', label= 'Data')
    fit6 = plt.plot(values6,k(values6),label = 'Linear fit')
    plt.xlabel('$\\frac{\Delta Re}{Re}$')
    plt.ylabel('$\\frac{\Delta (\lambda T_0)}{\lambda T_0}$')
    plt.grid(which = 'both', axis = 'both', linewidth = 0.2, linestyle = '--')
    plt.title('$\\frac{\Delta\lambda T_0}{\lambda T_0} \, vs \, \\frac{\Delta Re}{Re}$')
    plt.text(0.7,0.23, '$\\frac{\Delta\lambda T_0}{\lambda T_0} = A \\frac{\Delta Re}{Re}+ B$ \n $A = %.2f \pm %.2f$ \n $r=%.2f$ \n $b = %.2f$'%(slope6, std_err6,r_value6,intercept6), horizontalalignment = 'center', verticalalignment = 'center',transform = ax.transAxes, bbox = dict(facecolor = 'white'))
    plt.legend(shadow = 'true',loc=2)
    plt.show()


    
    # plot Delta_L/L vs Delta_R/R
    plt.plot(Reynolds_Errors/Reynolds_Numbers,Lyapunov_Errors/Lyapunov_Exponents, 'o', label= 'Data')
    fit6 = plt.plot(values6,p(values6),label = 'Linear fit')
    plt.xlabel('$\\frac{\Delta Re}{Re}$')
    plt.ylabel('$\\frac{\Delta (\lambda)}{\lambda}$')
    plt.grid(which = 'both', axis = 'both', linewidth = 0.2, linestyle = '--')
    plt.title('$\\frac{\Delta\lambda}{\lambda} \, vs \, \\frac{\Delta Re}{Re}$')
    plt.text(0.7,0.23, '$\\frac{\Delta\lambda}{\lambda} = A \\frac{\Delta Re}{Re}+ B$ \n $A = %.2f \pm %.2f$ \n $r=%.2f$ \n $b = %.2f$'%(slope11, std_err11,r_value11,intercept11), horizontalalignment = 'center', verticalalignment = 'center',transform = ax.transAxes, bbox = dict(facecolor = 'white'))
    plt.legend(shadow = 'true',loc=2)
    plt.show()

    '''

if __name__ == '__main__':
    main()
