#!/usr/bin/python3

import numpy as np
import sys
import matplotlib.pyplot as plt
from glob import glob
from os import path
from scipy import stats
import math
import statsmodels.api as sm
from cycler import cycler

#--------This program measures Lambda*T_0 vs Re and Delta_Lambda*tau vs Re --------------- 



def main():
    
    print("--------------------------------------------------------\nThis program measures all Reynolds Numbers, T_0's and tau's from .stats file and calculates lambda and error_lambda from .ftle files and plots error_lambda*tau against Re as well as lambda*T_0 vs Re \n --------------------------------------------------------\n\n")
    
    
    
    # Asks directory in which files are located
    #dir_name = raw_input("Directory where *.stats and *.ftle files are located: ")             
    dir_name = sys.argv[1]
    # Defines arrays with all needed values

    Reynolds_Numbers    = np.array([])      
    Reynolds_Errors     = np.array([])

    Klmgrv_Microtimes   = np.array([])
    Klmgrv_Micro_Errors = np.array([])

    T_0s_b              = np.array([])
    T_0_Errors_b        = np.array([]) 
    
    T_0s                = np.array([])
    T_0_Errors          = np.array([]) 
    
    Lyapunov_Exponents  = np.array([])
    Lyapunov_Errors     = np.array([])

    Urms                = np.array([])
    Urms_Errors         = np.array([])

    L                   = np.array([])
    L_Errors            = np.array([])

#---------------------------------stats files process-------------------------------------





    # Creates path and measures number of files

    Path = path.join(dir_name,'*.stats')    # Creates path
    stats_files = glob(Path)
    stats_files.sort()                      # Sorts files by numerical/alph order
    no_files = len(stats_files)             # Gets how many files are processes


    
    # Now data from the stats files is processesed and stored into arrays

    for i in range(no_files):
        
        # Open file and read
        filename = open(stats_files[i],'r')
        ar = np.loadtxt(filename, usecols=(2,3,5,6,13,1))
        ar = ar.transpose()
        

        t_f = np.size(ar,1) - 1
        
        
    
        
        
        # Consider only steady state data
        
        diss   = ar[0,2000:t_f]  
        R_L    = ar[1,2000:t_f]
        u      = ar[2,2000:t_f]
        l      = ar[3,2000:t_f]
        e      = ar[5,2000:t_f]
        # Calculates all necessary parameters
        
        visc            = ar[4,2]
        epsilon         = np.mean(diss)
        epsilon_error   = np.std(diss)
        
        Li              = np.mean(l) 
        Li_error        = np.std(l)
        U               = np.mean(u)
        U_error         = np.std(u)
        E               = np.mean(e)
        E_error         = np.std(e)
        
        Reynolds_Numbers    = np.append(Reynolds_Numbers,np.mean(R_L))    
        Reynolds_Errors     = np.append(Reynolds_Errors,np.std(R_L))
        
        

        Klmgrv_Microtimes   = np.append(Klmgrv_Microtimes,np.sqrt(visc/epsilon))
        Klmgrv_Micro_Errors = np.append(Klmgrv_Micro_Errors,(np.sqrt(visc)/2*(epsilon)**(3/2))*epsilon_error) 
        
        T_0s_b              = np.append(T_0s_b, E/epsilon)
        T_0_Errors_b        = np.append(T_0_Errors_b, (E*epsilon_error/epsilon**2)+(E_error/epsilon))
        
        T_0s                = np.append(T_0s, Li/U)
        T_0_Errors          = np.append(T_0_Errors, (Li*U_error/U**2)+(Li_error/U))


        Urms                = np.append(Urms, U)
        Urms_Errors         = np.append(Urms_Errors,U_error)

        L                   = np.append(L, Li)
        L_Errors            = np.append(L_Errors,Li_error)


#---------------------------------ftle files process--------------------------------------
        



    # Creates path and measures number of files

    Path = path.join(dir_name,'*.ftle*')    # Creates path
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

        t_i = 51                  #  ~ equivalent to 15s, hardcoded, when the ftle's become stable 

        Lyapunov_sample     = ar[1,t_i:t_f]
        Lyapunov_Exponents  = np.append(Lyapunov_Exponents,np.mean(Lyapunov_sample))
        Lyapunov_Errors     = np.append(Lyapunov_Errors,np.std(Lyapunov_sample))



#--------------------------------------- set error bars -----------------------


    Error_Lyap_T        = Lyapunov_Errors*T_0s + T_0_Errors*Lyapunov_Exponents
    Error_DeltaL_tau    = Lyapunov_Errors * Klmgrv_Micro_Errors
    Error_DeltaL_T      = Lyapunov_Errors * T_0_Errors 
    Error_DeltaLyap     = (Lyapunov_Errors ** 2)/(Lyapunov_Exponents ** 2)


# T = E/eps


    Error_Lyap_T_b        = Lyapunov_Errors*T_0s_b + T_0_Errors_b*Lyapunov_Exponents
    Error_DeltaL_T_b      = Lyapunov_Errors * T_0_Errors_b 
#-------------------------  define fit functions -----------------------------

    # fit lambda*T_0 ~ Re^alpha

    Log_R  = np.log10(Reynolds_Numbers)
    Log_Ly = np.log10(Lyapunov_Exponents*T_0s)

    Log_Ly_b = np.log10(Lyapunov_Exponents*T_0s_b)
    Log_E_Ly_b = np.log10(Lyapunov_Errors * T_0s_b)


                        # L/U 
    
    Log_R = sm.add_constant(Log_R) 
    mod_wls = sm.WLS(Log_Ly, Log_R, weights = ((Lyapunov_Exponents*T_0s)/Error_Lyap_T)**2)
    res_wls = mod_wls.fit()
    intercept, slope = res_wls.params
    intercept_e,std_err = res_wls.bse
    r_value = res_wls.rsquared    

    print('alpha_(L/U) = %.3f +/- %.3f' %(slope,std_err))
    print('intercept = ', intercept)

                        # E/eps 
    
    mod_wls_b = sm.WLS(Log_Ly_b, Log_R, weights = ((Lyapunov_Exponents*T_0s_b)/Error_Lyap_T_b)**2)
    res_wls_b = mod_wls_b.fit()
    intercept_b, slope_b = res_wls_b.params
    intercept_e_b,std_err_b = res_wls_b.bse
    r_value_b = res_wls_b.rsquared    
    
    print('alpha_(E/eps) = %.3f +/- %.3f' %(slope_b,std_err_b))
    print('intercept = ',intercept_b)
    
    # linear regression on Delta_Lyap vs Lyap
    
    slope2, intercept2, r_value2, p_value2 , std_err2 = stats.linregress(Lyapunov_Exponents, Lyapunov_Errors)
    

    # linear regression on Delta_Lyap*T_0 vs Re 

    Log_Ly_Er = np.log10(Lyapunov_Errors*T_0s)


    Log_Ly_Er_b = np.log10(Lyapunov_Errors*T_0s_b)

    # WLS fit 
                #L/U
    mod_wls3 = sm.WLS(Log_Ly_Er, Log_R)
    res_wls3 = mod_wls3.fit()
    intercept3, slope3 = res_wls3.params
    intercept_e3,std_err3 = res_wls3.bse
    r_value3 = res_wls3.rsquared   
                #E/eps
    mod_wls3_b = sm.WLS(Log_Ly_Er_b, Log_R)
    res_wls3_b = mod_wls3_b.fit()
    intercept3_b, slope3_b = res_wls3_b.params
    intercept_e3_B,std_err3_b = res_wls3_b.bse
    r_value3_b = res_wls3_b.rsquared   

    print('gamma_L/U = %.3f +/- %.3f' %(slope3,std_err3))
    
    print('gamma_E/eps = %.3f +/- %.3f' %(slope3_b,std_err3_b))


    # linear regression on Delta_Re vs Re


    slope4, intercept4, r_value4, p_value4, std_err4 = stats.linregress(Reynolds_Numbers, Reynolds_Errors)


    # linear regression on Delta_T vs T


    slope5, intercept5, r_value5, p_value5, std_err5 = stats.linregress(T_0s,T_0_Errors)

    print('T_0_prec= %.3f +/- %.3f'%(slope5,std_err5))
    #linear regression on Delta Lyap /Lyap vs Delta_Re/Re


    slope6, intercept6, r_value6, p_value6, std_err6 = stats.linregress(Reynolds_Errors/Reynolds_Numbers,Error_Lyap_T/(Lyapunov_Exponents*T_0s))



    #linear regression on Delta Lyap_T vs Lyap_T 


    slope8, intercept8, r_value8, p_value8, std_err8 = stats.linregress(Lyapunov_Exponents*T_0s,Error_Lyap_T)


    #linear regression on Urms vs Delta_Urms


    slope9, intercept9, r_value9, p_value9, std_err9 = stats.linregress(Urms,Urms_Errors)



    #linear regression on L vs Delta_L


    slope10, intercept10, r_value10, p_value10, std_err10 = stats.linregress(L,L_Errors)


    #linear regression on Delta_Lyap/ Lyap vs Delta_Re/Re

    slope11, intercept11, r_value11, p_value11, std_err11 = stats.linregress(Reynolds_Errors/Reynolds_Numbers,Lyapunov_Errors/Lyapunov_Exponents)
    
    

#-------------------------------- Define fit functions -------------------------- 

    # Set fit function to be plot
    def f(t):
        return 10**(intercept)*t**(slope)

    def f_b(t):
        return 10**(intercept_b)*t**(slope_b)

    low = np.amin(Reynolds_Numbers)
    maxi = np.amax(Reynolds_Numbers)
    values = np.arange(low,maxi,0.1)
    intermediate = 10**((np.log10(maxi)+np.log10(low))/2)
    
    def f_Ruelle(t):
        return (f(intermediate))*(intermediate**(-0.5))*(t**(0.5))

    def f_Ruelle_b(t):
        return (f_b(intermediate))*(intermediate**(-0.5))*(t**(0.5))

    def g(t):
        return slope2*t + intercept2 

    low2 = np.amin(Lyapunov_Exponents)
    maxi2 = np.amax(Lyapunov_Exponents)
    values2 = np.arange(low2,maxi2,0.001)
    intermediatei2 = (maxi2 + low2)/2

    def h(t):
        return 10**(intercept3)*t**(slope3)
    
    def h_b(t):
        return 10**(intercept3_b)*t**(slope3_b)
   
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





#--------------------------- plot everything --------------------
    
    
    params = {'legend.fontsize':'x-large','axes.labelsize':'xx-large','xtick.labelsize':'x-large','ytick.labelsize':'x-large'}
    plt.rcParams.update(params)
    plt.rc('text',usetex=True)
    plt.rc('font',family='serif')
    plt.rcParams['axes.prop_cycle'] = cycler(color = 'krgcmyb')
   #plt.style.use('seaborn-muted')
    fig,ax =  plt.subplots()

    # plot Lyap vs Re and linear fit   T_0 = L/U
    data = plt.errorbar(Reynolds_Numbers, Lyapunov_Exponents*T_0s, xerr = Reynolds_Errors, yerr = Error_Lyap_T, marker ='.', markerfacecolor = 'None',linestyle = 'None',label ='Data',ecolor = '0.6', capsize = 3, elinewidth = 1, capthick = 1)
    
    fit_ = plt.plot(values, f(values), color = '#bc5c47',label = 'Linear fit')

    
    fit3_ = plt.plot(values, f_Ruelle(values), color = '0.7', linestyle = ':', label = 'Ruelle\'s prediction')


    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$Re$')
    plt.ylabel('$\lambda  \, T_0 \quad (T_0 = \\frac{L}{U})$')
    plt.grid(which='both', axis='both',linewidth = 0.2, linestyle = '--')
    plt.legend(loc=2,frameon = False)
    plt.show()
   

    # plot Lyap vs Re and linear fit  T_0 = E/eps
    data = plt.errorbar(Reynolds_Numbers, Lyapunov_Exponents*T_0s_b, xerr = Reynolds_Errors, yerr = Error_Lyap_T_b, marker ='.', markerfacecolor = 'None',linestyle = 'None',label ='Data',ecolor = '0.6', capsize = 3, elinewidth = 1, capthick = 1)
    

    fit2_ = plt.plot(values,f_b(values),color = '#bc5c47',label = 'Linear fit')
    
    fit4_ = plt.plot(values,f_Ruelle_b(values),color = '0.7', linestyle = ':', label = 'Ruelle\'s prediction')
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$Re$')
    plt.ylabel('$\lambda  \, T_0 \quad (T_0 = \\frac{E}{\\varepsilon})$')
    plt.grid(which='both', axis='both',linewidth = 0.2, linestyle = '--')
    plt.legend(loc=2,frameon = False)
    plt.show()

    # plot Delta_Lyap vs Lyap and linear fit
    

    data2 = plt.plot(Lyapunov_Exponents,Lyapunov_Errors, marker ='.', markerfacecolor = 'None',linestyle = 'None',label = 'Data')
    fit2 = plt.plot(values2,g(values2),color = '#bc5c47',label = 'Linear fit')
    plt.xlabel('$\lambda$')
    plt.ylabel('$\sigma_{\lambda}$')
    plt.grid(which='both', axis='both',linewidth = 0.2, linestyle = '--')
    plt.legend(frameon = False,loc=2)
    plt.show()

    
    
    # plot Delta_Lyap*T_0s vs Re and linear fit
    
    
    
    
    plt.errorbar(Reynolds_Numbers,Lyapunov_Errors*T_0s, xerr = Reynolds_Errors,yerr = Error_DeltaL_T, marker ='.', markerfacecolor = 'None',linestyle = 'None',ecolor = '0.6', capsize= 3,elinewidth = 1, capthick = 1, label = 'Data')
    fit3 = plt.plot(values,h(values), color = '#bc5c47',label = 'Linear fit')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$Re$')
    plt.ylabel('$\sigma_{\lambda} \, T_0 \quad (T_0 = \\frac{L}{U})$')
    plt.grid(which='both', axis='both',linewidth = 0.2, linestyle = '--')
    plt.legend(frameon = False,loc=2)
    plt.show()

    
    # plot Delta_Lyap*T_0s vs Re and linear fit
    
    
    
    
    plt.errorbar(Reynolds_Numbers,Lyapunov_Errors*T_0s_b, xerr = Reynolds_Errors,yerr = Error_DeltaL_T_b, marker ='.', markerfacecolor = 'None',linestyle = 'None',ecolor = '0.6', capsize= 3,elinewidth = 1, capthick = 1, label = 'Data ')
    fit3 = plt.plot(values,h_b(values), color = '#bc5c47',label = 'Linear fit')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$Re$')
    plt.ylabel('$\sigma_{\lambda} \, T_0 \quad (T_0 = \\frac{E}{\\varepsilon})$')
    plt.grid(which='both', axis='both',linewidth = 0.2, linestyle = '--')
    plt.legend(frameon = False,loc=2)
    plt.show()


    # plot Delta_Lyap*tau vs Re

    plt.errorbar(Reynolds_Numbers,Lyapunov_Errors*Klmgrv_Microtimes, xerr = Reynolds_Errors, yerr = Error_DeltaL_tau,marker ='.', markerfacecolor = 'None',linestyle = 'None', ecolor = 'g', capsize = 3, elinewidth = 1, capthick = 1)
    plt.xlabel('$Re$')
    plt.ylabel('$\sigma_{\lambda} \, \\tau$')
    plt.xscale('log')
    plt.grid(which='both', axis='both',linewidth = 0.2, linestyle = '--')
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
    
   # plt.plot(Lyapunov_Exponents*T_0s,Error_Lyap_T, marker ='.', markerfacecolor = 'None',linestyle = 'None', label = 'Data ')
   # fit8 = plt.plot(values8,m(values8),color = '#bc5c47' ,label = 'Linear fit')
   # #plt.xscale('log')
   # #plt.yscale('log')
   # plt.xlabel('$\lambda \, T_0$')
   # plt.ylabel('$\sigma_{\lambda T_0}$')
   # plt.grid(which='both', axis='both',linewidth = 0.2, linestyle = '--')
   # plt.legend(frameon = False,loc=2)
   # plt.show()



    
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
    fig5, ax5 = plt.subplots()
    ax5.plot(Reynolds_Numbers, T_0_Errors/T_0s, marker ='.', markerfacecolor = 'None',linestyle = 'None', label = 'Data')
   # fit5 = ax5.plot(values5,j(values5),label = 'Linear fit',color = '#bc5c47')
    plt.xlabel('$T_0$')
    plt.ylabel('$\sigma_{T_0}$')
    plt.grid(which='both', axis='both',linewidth = 0.2, linestyle = '--')
#    plt.title('$\Delta T_0 vs T_0$')
#    plt.text(0.5,0.1,'$\Delta T_0  = m T_0 + b$ \n $ m = %.2f \pm %.2f$ \n $r= %.2f$ \n $b= %.2f$'%(slope5,std_err5,r_value5,intercept5), horizontalalignment = 'center', verticalalignment = 'center',transform = ax.transAxes, bbox = dict(facecolor = 'white'))
   # plt.legend(frameon = False,loc=2)
    plt.show()
    
    
    
    print('T_0 = %.2f +/- %.2f' %(np.mean(T_0s),np.std(T_0s)))
    
    
    
    
    
    # plot Delta_Re vs Re 
    #fig4, ax4 = plt.subplots()
    plt.plot(Reynolds_Numbers, Reynolds_Errors,marker ='.', markerfacecolor = 'None',linestyle = 'None' , label = 'Data')
    fit4 = plt.plot(values,i(values),color = '#bc5c47',label = 'Linear fit')
    plt.xlabel('$Re$')
    plt.ylabel('$\sigma_{Re}$')
    plt.grid(which='both', axis='both',linewidth = 0.2, linestyle = '--')
    plt.legend(frameon = False,loc=2)
    plt.show()

    # plot Delta_L_T/L_T vs Delta_R/R
#    plt.plot(Reynolds_Errors/Reynolds_Numbers,Error_Lyap_T/(Lyapunov_Exponents*T_0s), '.', label= 'Data')
#    fit6 = plt.plot(values6,k(values6),label = 'Linear fit')
#    plt.xlabel('$\\frac{\sigma_{Re}}{Re}$')
#    plt.ylabel('$\\frac{\sigma{\lambda T_0})}{\lambda T_0}$')
#    plt.grid(which = 'both', axis = 'both', linewidth = 0.2, linestyle = '--')
#    plt.legend(shadow = 'true',loc=2)
#    plt.show()


    
    # plot Delta_L/L vs Delta_R/R
    plt.plot(Reynolds_Errors/Reynolds_Numbers,Lyapunov_Errors/Lyapunov_Exponents, marker ='.', markerfacecolor = 'None',linestyle = 'None', label= 'Data')
    fit6 = plt.plot(values6,p(values6),color = '#bc5c47',label = 'Linear fit')
    plt.xlabel('$\\frac{\sigma_{Re}}{Re}$')
    plt.ylabel('$\\frac{\sigma_{\lambda}}{\lambda}$')
    plt.grid(which = 'both', axis = 'both', linewidth = 0.2, linestyle = '--')
    plt.legend(frameon = False,loc=2)
    plt.show()



if __name__ == '__main__':
    main()
