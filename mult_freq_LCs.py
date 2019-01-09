#Smart and Fermi-LAT data

import numpy as np
import numpy
import matplotlib.pyplot as plt


#Fermi-LAT data

f1 = open('/home/yerevan2/Desktop/3C454.3_86400_t1.txt').read()
f2= open('3c454_gamma.txt', 'w')
replacements = {'"':'', ',':'  ','T':'1', 'F':'0'}
for i in replacements.keys():
    f1 = f1.replace(i, replacements[i])
f2.write(f1)
f2.close


Tmin_MET, Tmax_MET, Flux, Flux_err, Upper_Limit = np.loadtxt('3c454_gamma.txt', skiprows=0, unpack=True)
Tmid_MET=(Tmin_MET+Tmax_MET)/2

Tmin_MJD=Tmin_MET/86400+51910
Tmid_MJD=Tmid_MET/86400+51910
Tmax_MJD=Tmax_MET/86400+51910

FGscale=10**5
Flux=FGscale*Flux
print max(Flux)
Flux_err=FGscale*Flux_err

#SMART data
MJD_B, FlambdaB, FlambdaBerr = np.loadtxt('B_band.txt', skiprows=0, unpack=True)
MJD_V, FlambdaV, FlambdaVerr = np.loadtxt('V_band.txt', skiprows=0, unpack=True)
MJD_R, FlambdaR, FlambdaRerr = np.loadtxt('R_band.txt', skiprows=0, unpack=True)
MJD_J, FlambdaJ, FlambdaJerr = np.loadtxt('J_band.txt', skiprows=0, unpack=True)
MJD_K, FlambdaK, FlambdaKerr = np.loadtxt('K_band.txt', skiprows=0, unpack=True)

fig, axes = plt.subplots(6,sharex=True)
axes[0].errorbar(Tmid_MJD, Flux, yerr=Flux_err, xerr=0.5, fmt='o', color='black', ecolor='black',ms=3,capsize=0)
axes[0].set_title('2008Jun23 - 2015Nov30')
axes[0].set_ylabel('F (10^-6 ph cm^-2 s^-1)')
axes[0].axis([min(Tmid_MJD)-20,max(Tmid_MJD)+20,min(Flux)-1,max(Flux)+1])
axes[1].errorbar(MJD_B, FlambdaB, FlambdaVerr, fmt='o', color='blue', ecolor='blue',ms=3,capsize=0)
axes[1].set_ylabel('B(mJy)')
axes[1].axis([min(Tmid_MJD)-20,max(Tmid_MJD)+20,min(FlambdaB)-1,max(FlambdaB)+1])
axes[2].errorbar(MJD_V, FlambdaV, FlambdaVerr, fmt='o', color='red', ecolor='red',ms=3,capsize=0)
axes[2].set_ylabel('V(mJy)')
axes[2].axis([min(Tmid_MJD)-20,max(Tmid_MJD)+20,min(FlambdaV)-4,max(FlambdaV)+4])
axes[3].errorbar(MJD_R, FlambdaR, FlambdaRerr, fmt='o', color='Magenta', ecolor='Magenta',ms=3,capsize=0)
axes[3].set_ylabel('R(mJy)')
axes[3].axis([min(Tmid_MJD)-20,max(Tmid_MJD)+20,min(FlambdaR)-4,max(FlambdaR)+4])
axes[4].errorbar(MJD_J, FlambdaJ, FlambdaJerr, fmt='o', color='g', ecolor='g',ms=3,capsize=0)
axes[4].set_ylabel('J(mJy)')
axes[4].axis([min(Tmid_MJD)-20,max(Tmid_MJD)+20,min(FlambdaJ)-4,max(FlambdaJ)+4])
axes[5].errorbar(MJD_K, FlambdaK, FlambdaKerr, fmt='o', color='cyan', ecolor='cyan',ms=3,capsize=0)
axes[5].set_ylabel('K (mJy)')
axes[5].set_xlabel('Time (MJD)')
axes[5].axis([min(Tmid_MJD)-20,max(Tmid_MJD)+20,min(FlambdaK)-10,max(FlambdaK)+10])
plt.subplots_adjust(left=0.055, bottom=0.06, right=0.99, top=0.94, hspace=0)
plt.show() 

