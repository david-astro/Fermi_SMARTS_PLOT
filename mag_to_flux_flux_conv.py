#Reddening corrections and M to flux conversation for the NIR-optical photometric data from the SMART#


import numpy as np
import numpy
import matplotlib.pyplot as plt


#Define variables for the first data
JD, JD_B, B, Berr, JD_V, V, Verr, JD_R, R, Rerr, JD_J, J, Jerr, JD_K, K, Kerr = np.loadtxt('3C454.txt', skiprows=0, unpack=True)

#Start time for each observation subtracted the start time of the observations

#JD_B=[np.subtract(x,min(JD)) for x in JD_B]
#JD_V=[np.subtract(x,min(JD)) for x in JD_V]
#JD_R=[np.subtract(x,min(JD)) for x in JD_R]
#JD_J=[np.subtract(x,min(JD)) for x in JD_J]
#JD_K=[np.subtract(x,min(JD)) for x in JD_K]

#Set the missing observations to zero
MJD=[0.0 if x==999.0 else x-2400000.5 for x in JD]
MJD_B=[0.0 if x==999.0 else x-2400000.5 for x in JD_B]
MJD_V=[0.0 if x==999.0 else x-2400000.5 for x in JD_V]
MJD_R=[0.0 if x==999.0 else x-2400000.5 for x in JD_R]
MJD_J=[0.0 if x==999.0 else x-2400000.5 for x in JD_J]
MJD_K=[0.0 if x==999.0 else x-2400000.5 for x in JD_K]

B=[0.0 if x==999.0 else x for x in B]
Berr=[0.0 if x==999.0 else x for x in Berr]
V=[0.0 if x==999.0 else x for x in V]
Verr=[0.0 if x==999.0 else x for x in Verr]
R=[0.0 if x==999.0 else x for x in R]
Rerr=[0.0 if x==999.0 else x for x in Rerr]
J=[0.0 if x==999.0 else x for x in J]
Jerr=[0.0 if x==999.0 else x for x in Jerr]
K=[0.0 if x==999.0 else x for x in K]
Kerr=[0.0 if x==999.0 else x for x in Kerr]


#Correction parameters of magnitudes of all bands (Cardeli et al. 1989, Schlafy & Finkbeiner 2011)
Ebv=0.0889
Rv=3.1
Av=Ebv*Rv

Alambda_B=1.337*Av
Alambda_V=1.0*Av
Alambda_R=0.751*Av
Alambda_J=0.282*Av
Alambda_K=0.114*Av

#Calculated by "Absorption Law Calculator" (Cardelli et al. 1989, http://www.dougwelch.org/Acurve.html)
#Alambda_B=0.36705
#Alambda_V=0.27819
#Alambda_R=0.23188
#Alambda_J=0.08083
#Alambda_K=0.31515

#Corrected magnitudes
corrB=[0.0 if x==-Alambda_B else x for x in np.subtract(B,Alambda_B)]
corrBmax=[0.0 if x==-Alambda_B else x for x in np.add(corrB,Berr)]
corrBmin=[0.0 if x==-Alambda_B else x for x in np.subtract(corrB,Berr)]
corrV=[0.0 if x==-Alambda_V else x for x in np.subtract(V,Alambda_V)]
corrVmax=[0.0 if x==-Alambda_V else x for x in np.add(corrV,Verr)]
corrVmin=[0.0 if x==-Alambda_V else x for x in np.subtract(corrV,Verr)]
corrR=[0.0 if x==-Alambda_R else x for x in np.subtract(R,Alambda_R)]
corrRmax=[0.0 if x==-Alambda_R else x for x in np.add(corrR,Rerr)]
corrRmin=[0.0 if x==-Alambda_R else x for x in np.subtract(corrR,Rerr)]
corrJ=[0.0 if x==-Alambda_J else x for x in np.subtract(J,Alambda_J)]
corrJmax=[0.0 if x==-Alambda_J else x for x in np.add(corrJ,Jerr)]
corrJmin=[0.0 if x==-Alambda_J else x for x in np.subtract(corrJ,Jerr)]
corrK=[0.0 if x==-Alambda_K else x for x in np.subtract(K,Alambda_K)]
corrKmax=[0.0 if x==-Alambda_K else x for x in np.add(corrK,Kerr)]
corrKmin=[0.0 if x==-Alambda_K else x for x in np.subtract(corrK,Kerr)]


#Flux densities in the unit of mJy (10^-3 ergs cm^-2 sec^-1 Hz^-1) using absolute fluxes from Bessell et al. (1998).

#Effective wavelengths for UBVRIJHKL Cousins-Glass-Johnson system
#                    B         V         R         J         K
#Lmbda_eff (nm)      438       545       641       1220      2190    (Exact values)
#Freq_eff (GHz)      684.46    550.08    467.69    245.73    136.89  (Rounded)
#Energy (10^3 eV)    2.83      2.28      1.93      1.02      0.57    (Rounded)


#Absolute fluxes (corresponding to zero magnitude) for UBVRIJHKL Cousins-Glass-Johnson system
Fzero_B=4.063*10**-20
Fzero_V=3.636*10**-20
Fzero_R=3.064*10**-20
Fzero_J=1.589*10**-20
Fzero_K=0.640*10**-20

Fscale=10**26

FlambdaB=[Fscale*Fzero_B*np.power(10,-np.divide(x,2.5)) if x>1 else x for x in corrB]
FlambdaBmax=[Fscale*Fzero_B*np.power(10,-np.divide(x,2.5)) if x>1 else x for x in corrBmin]
FlambdaBmin=[Fscale*Fzero_B*np.power(10,-np.divide(x,2.5)) if x>1 else x for x in corrBmax]
FlambdaBerr=np.power(np.divide(np.add(np.power(np.subtract(FlambdaBmin,FlambdaB),2),np.power(np.subtract(FlambdaBmax,FlambdaB),2)),2),0.5)

#V band        
FlambdaV=[Fscale*Fzero_V*np.power(10,-np.divide(x,2.5)) if x>1 else x for x in corrV]
FlambdaVmax=[Fscale*Fzero_V*np.power(10,-np.divide(x,2.5)) if x>1 else x for x in corrVmin]
FlambdaVmin=[Fscale*Fzero_V*np.power(10,-np.divide(x,2.5)) if x>1 else x for x in corrVmax]
FlambdaVerr=np.power(np.divide(np.add(np.power(np.subtract(FlambdaVmin,FlambdaV),2),np.power(np.subtract(FlambdaVmax,FlambdaV),2)),2),0.5)

#R band
FlambdaR=[Fscale*Fzero_R*np.power(10,-np.divide(x,2.5)) if x>1 else x for x in corrR]
FlambdaRmax=[Fscale*Fzero_R*np.power(10,-np.divide(x,2.5)) if x>1 else x for x in corrRmin]
FlambdaRmin=[Fscale*Fzero_R*np.power(10,-np.divide(x,2.5)) if x>1 else x for x in corrRmax]
FlambdaRerr=np.power(np.divide(np.add(np.power(np.subtract(FlambdaRmin,FlambdaR),2),np.power(np.subtract(FlambdaRmax,FlambdaR),2)),2),0.5)

#J band
FlambdaJ=[Fscale*Fzero_J*np.power(10,-np.divide(x,2.5)) if x>1 else x for x in corrJ]
FlambdaJmax=[Fscale*Fzero_J*np.power(10,-np.divide(x,2.5)) if x>1 else x for x in corrJmin]
FlambdaJmin=[Fscale*Fzero_J*np.power(10,-np.divide(x,2.5)) if x>1 else x for x in corrJmax]
FlambdaJerr=np.power(np.divide(np.add(np.power(np.subtract(FlambdaJmin,FlambdaJ),2),np.power(np.subtract(FlambdaJmax,FlambdaJ),2)),2),0.5)

#K band
FlambdaK=[Fscale*Fzero_K*np.power(10,-np.divide(x,2.5)) if x>1 else x for x in corrK]
FlambdaKmax=[Fscale*Fzero_K*np.power(10,-np.divide(x,2.5)) if x>1 else x for x in corrKmin]
FlambdaKmin=[Fscale*Fzero_K*np.power(10,-np.divide(x,2.5)) if x>1 else x for x in corrKmax]
FlambdaKerr=np.power(np.divide(np.add(np.power(np.subtract(FlambdaKmin,FlambdaK),2),np.power(np.subtract(FlambdaKmax,FlambdaK),2)),2),0.5)


#Transfer the results
output = '\n'.join('\t'.join(map(str,row)) for row in zip(MJD_B,FlambdaB,FlambdaBerr))
with open('B_band.txt', 'w') as f:
    f.write(output)

output = '\n'.join('\t'.join(map(str,row)) for row in zip(MJD_V,FlambdaV,FlambdaVerr))
with open('V_band.txt', 'w') as f:
    f.write(output)

output = '\n'.join('\t'.join(map(str,row)) for row in zip(MJD_R,FlambdaR,FlambdaRerr))
with open('R_band.txt', 'w') as f:
    f.write(output)

output = '\n'.join('\t'.join(map(str,row)) for row in zip(MJD_J,FlambdaJ,FlambdaJerr))
with open('J_band.txt', 'w') as f:
    f.write(output)

output = '\n'.join('\t'.join(map(str,row)) for row in zip(MJD_K,FlambdaK,FlambdaKerr))
with open('K_band.txt', 'w') as f:
    f.write(output)

#Plot output  
fig, axes = plt.subplots(5,sharex=True)
axes[0].errorbar(MJD_B, FlambdaB, FlambdaVerr, fmt='o', color='blue', ecolor='blue',ms=3,capsize=0)
axes[0].set_title('2008Jun23 - 2015Nov30')
axes[0].set_ylabel('B')
axes[0].axis([min(MJD)-20,max(JD)+20,min(FlambdaB)-1,max(FlambdaB)+1])
axes[1].errorbar(MJD_V, FlambdaV, FlambdaVerr, fmt='o', color='red', ecolor='red',ms=3,capsize=0)
axes[1].set_ylabel('V')
axes[1].axis([min(MJD)-20,max(JD)+20,min(FlambdaV)-4,max(FlambdaV)+4])
axes[2].errorbar(MJD_R, FlambdaR, FlambdaRerr, fmt='o', color='Magenta', ecolor='Magenta',ms=3,capsize=0)
axes[2].set_ylabel('R')
axes[2].axis([min(MJD)-20,max(JD)+20,min(FlambdaR)-4,max(FlambdaR)+4])
axes[3].errorbar(MJD_J, FlambdaJ, FlambdaJerr, fmt='o', color='g', ecolor='g',ms=3,capsize=0)
axes[3].set_ylabel('J')
axes[3].axis([min(MJD)-20,max(JD)+20,min(FlambdaJ)-4,max(FlambdaJ)+4])
axes[4].errorbar(MJD_K, FlambdaK, FlambdaKerr, fmt='o', color='cyan', ecolor='cyan',ms=3,capsize=0)
axes[4].set_ylabel('K')
axes[4].set_xlabel('Time (MJD)')
axes[4].axis([min(MJD)-20,max(JD)+20,min(FlambdaK)-10,max(FlambdaK)+10])
plt.subplots_adjust(left=0.055, bottom=0.06, right=0.99, top=0.94, hspace=0)
plt.show() 

