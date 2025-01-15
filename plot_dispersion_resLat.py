# -*- coding: utf-8 -*-
"""
Created on Sun Aug 11 05:13:17 2024

@author: zhousu
"""

#import required packages
import numpy as np
import os,numlib,refr_index_full_left
import sys
import matplotlib.pyplot as plt
import matplotlib,space_science

#import WPIT modules
import WPIT.Environment_mod as env
import WPIT.WaveProperties_mod as wave
import WPIT.WPI_mod.EMIC_ion_mod as wpi
 
#input parameters
L_shell=4
density_ratio=15
Byw0=3*10**(-9)  # y-component of the wave magnetic field (3nT)
theta0_deg=0  # wave normal angle
ion_com=np.array([77,20,3])/100
w_wave_norm=0.96
save_dir=r'E:\mywork\my_drafts\nonlinear_affected_by_EMIC_parameters\figs/'
save_fig=True

xlim_kappa=[0,1]
ylim_freq=[0,1]
xlim_Ek=[1e-1,1e6]

alpha_grid=0.36
LineWidth=1.5;LabelSize=20;TickWidth=1.6;TickLen=6
panel_label=['(a)','(b)','(c)','(d)','(e)','(f)']



plt.close('all')
num_row=1
num_column=3
myfig=plt.figure(figsize=(21,7.2))
hf=numlib.axes_create(row=num_row, column=num_column,top=0.08,bottom=0.12,left=0.045,right=0.06,gap_row=0.045,gap_column=0.06)

lat_eq_rad=np.deg2rad(0) #equatorial magnetic latitude (lamda=0 deg)
Beq =env.Bmag_dipole(L_shell,lat_eq_rad) #equatorial magnetic field strength


wce_eq=env.omega_cyclotron(Beq,env.const.qe,env.const.me) #equatorial electron cyclotron frequency
wcHe_eq=env.omega_cyclotron(Beq,env.const.qe,env.const.mHe) #equatorial helium cyclotron frequency
wcH_eq=env.omega_cyclotron(Beq,env.const.qe,env.const.mH) #equatorial hydrogen cyclotron frequency
wcO_eq=env.omega_cyclotron(Beq,env.const.qe,env.const.mO) #equatorial oxygen cyclotron frequency
wave_freq=np.arange(0.0001,1,0.0001)*wcH_eq


wpe_eq=density_ratio*wce_eq #equatorial plasma frequency
ne_eq=(env.const.me*env.const.epsilon0*wpe_eq*wpe_eq)/(env.const.qe*env.const.qe) #calculate equatorial electron density from equatorial plasma frequency
w_proton=np.sqrt(ne_eq*env.const.qi**2/(env.const.mH*env.const.epsilon0))#proton plasma frequency
#env.const.
#Ne=ne_eq

#ne_eq=Ne*1e6
nH_eq=ion_com[0]*ne_eq #77% hydrogen
nHe_eq=ion_com[1]*ne_eq #20% helium
nO_eq=ion_com[2]*ne_eq #3% oxygen
theta0_rad=np.deg2rad(theta0_deg) #convert wave normal angle to rad

aRef_plus=np.ones(len(wave_freq))*np.nan  #refractive index
aRef_minus=aRef_plus.copy()
aRef_left=aRef_plus.copy()
aRef_right=aRef_plus.copy()

#aLat=aRef_plus.copy() 
plr_plus=aRef_plus.copy()
plr_minus=aRef_plus.copy()
aD_stix=aRef_plus.copy()
aKappa_left=aRef_plus.copy()
aKappa_right=aRef_plus.copy()

for iFreq in range(len(wave_freq)):
    aWave_freq=wave_freq[iFreq]
    S0,D0,P0,R0,L0=wave.stix_parameters(aWave_freq, ne_eq, nH_eq, nHe_eq, nO_eq, Beq)
    eta_sq_plus0,eta_sq_minus0,mu0,kappa0,kappaz0,kappax0=wave.refr_index_full(theta0_rad,aWave_freq,S0,P0,R0,L0)
    plr_plus[iFreq]=-(eta_sq_plus0-S0)/D0
    plr_minus[iFreq]=-(eta_sq_minus0-S0)/D0
    ref_left0,ref_right0,eta_sq_plus0,eta_sq_minus0,kappa0,kappaz0,kappax0=refr_index_full_left.refr_index_full_left(theta0_rad,aWave_freq,S0,D0,P0,R0,L0)

    aRef_left[iFreq]=ref_left0
    aRef_right[iFreq]=ref_right0
    aD_stix[iFreq]=D0
    aKappa_left[iFreq]=ref_left0*aWave_freq/env.const.c_light
    aKappa_right[iFreq]=ref_right0*aWave_freq/env.const.c_light
    


#sys.exit()
#find the cutoff frequency of He-band EMIC waves
fp=(wave_freq>=wcO_eq)&(wave_freq<wcHe_eq)
bRef_left=aRef_left[fp]
bWave_freq=wave_freq[fp]
fp=~np.isnan(bRef_left)
He_cutoff_freq=bWave_freq[fp][0]

#find the cutoff frequency of H-band EMIC waves
fp=(wave_freq>=wcHe_eq)&(wave_freq<wcH_eq)
cRef_left=aRef_left[fp]
cWave_freq=wave_freq[fp]
fp=~np.isnan(cRef_left)
H_cutoff_freq=cWave_freq[fp][0]
#find the crossover frequency of He-band EMIC waves
fp=(wave_freq>=wcO_eq)&(wave_freq<wcHe_eq)
bD_stix=aD_stix[fp]
bWave_freq=wave_freq[fp]
bD_stix_cross=np.array([bD_stix[i]*bD_stix[i+1] for i in range(len(bD_stix)-1)])
fp_D=np.where(bD_stix_cross<0)
He_crossover=bWave_freq[fp_D[0][0]+1]
#find the crossover frequency of H-band EMIC waves
fp=(wave_freq>=wcHe_eq)&(wave_freq<wcH_eq)
cD_stix=aD_stix[fp]
cWave_freq=wave_freq[fp]
cD_stix_cross=np.array([cD_stix[i]*cD_stix[i+1] for i in range(len(cD_stix)-1)])
fp_D=np.where(cD_stix_cross<0)
H_crossover=cWave_freq[fp_D[0][0]+1]

wave_freq_left=wave_freq.copy()
wave_freq_right=wave_freq.copy()

fp=(wave_freq==He_crossover)|(wave_freq==H_crossover)
wave_freq_left[fp]=np.nan
wave_freq_right[fp]=np.nan

#plot dispersion relationship of L mode
axes_num=0
hf[axes_num].plot(aKappa_left/(w_proton/env.const.c_light),wave_freq_left/wcHe_eq,color='b')
hf[axes_num].set_xlim(xlim_kappa)
hf[axes_num].set_ylim(ylim_freq)
#plot dispersion relationship of R mode
axes_num=0
hf[axes_num].plot(aKappa_right/(w_proton/env.const.c_light),wave_freq_right/wcHe_eq,color='r')
hf[axes_num].set_xlim(xlim_kappa)
hf[axes_num].set_ylim(ylim_freq)

k=np.arange(0,xlim_kappa[1]*1.1,0.01)*(w_proton/env.const.c_light)

hf[axes_num].plot([xlim_kappa[0],xlim_kappa[1]],np.ones(2)*wcH_eq/wcHe_eq,color='k',ls='-')
hf[axes_num].plot([xlim_kappa[0],xlim_kappa[1]],np.ones(2)*wcHe_eq/wcHe_eq,color='k',ls='-')
hf[axes_num].plot([xlim_kappa[0],xlim_kappa[1]],np.ones(2)*wcO_eq/wcHe_eq,color='k',ls='-')
hf[axes_num].plot([xlim_kappa[0],xlim_kappa[1]],np.ones(2)*wce_eq/wcHe_eq,color='k',ls='-')

hf[axes_num].grid(alpha=alpha_grid)
hf[axes_num].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
hf[axes_num].yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
hf[axes_num].tick_params(which='major',labelsize=LabelSize,width=TickWidth,length=TickLen,tickdir='in',right=True,top=True)
hf[axes_num].tick_params(which='minor',color='k',width=TickWidth*0.5,length=TickLen*0.5,tickdir='in',right=True,top=True)

hf[axes_num].set_xlabel(r'$\rmk/(\omega_{H^+}/c)$',fontsize=LabelSize+1)
hf[axes_num].set_ylabel(r'$\rm\omega/\Omega_{He\_local}$',fontsize=LabelSize+1)

hf[axes_num].text(x=0.91*xlim_kappa[1],y=wcO_eq/wcHe_eq+0.03,s=r'$\rm\Omega_{O\_local}$',ha='center',va='center',fontsize=LabelSize)

hf[axes_num].text(x=0.5,y=1.04,s='Dispersion relationship, '+r'$\rm\psi=$'+str(theta0_deg)+r'$\rm\degree$',ha='center',va='center',fontsize=LabelSize+1,transform=hf[axes_num].transAxes)
hf[axes_num].text(x=0.5,y=0.8,s='L',ha='center',va='center',color='b',fontsize=LabelSize+2,weight='bold',transform=hf[axes_num].transAxes)
hf[axes_num].text(x=0.18,y=0.8,s='R',ha='center',va='center',color='r',fontsize=LabelSize+2,weight='bold',transform=hf[axes_num].transAxes)

hf[axes_num].annotate(text='crossover frequency',color='k',fontsize=LabelSize+1,xy=(0.14,0.4),\
    xycoords='axes fraction',xytext=(0.24,0.35),va='center',textcoords='axes fraction',arrowprops=dict(color='k',headwidth=9,width=1))

#panel 1, plot crossover frequency versus latitude
lat_dip=np.arange(0,90,0.01)
He_co_norm=np.ones(len(lat_dip))*np.nan
Lat_deg=He_co_norm.copy()
aWave_norm=He_co_norm.copy()

for iLat in range(len(lat_dip)):
    aLat_deg=lat_dip[iLat]
    aLat_rad=np.deg2rad(aLat_deg)
    aBfield=env.Bmag_dipole(L_shell,aLat_rad)
    wcHe=env.omega_cyclotron(aBfield,env.const.qe,env.const.mHe) #equatorial helium cyclotron frequency
    wcH=env.omega_cyclotron(aBfield,env.const.qe,env.const.mH)
    aHe_crossover=(He_crossover/wcHe_eq)*wcHe
    He_co_norm[iLat]=aHe_crossover/wcHe
    aWave_norm[iLat]=w_wave_norm*wcHe_eq/wcHe
    Lat_deg[iLat]=aLat_deg


fp=aWave_norm>=(He_crossover/wcHe_eq)

axes_num=1
hf[axes_num].plot(Lat_deg,He_co_norm,color='k',ls='--')
hf[axes_num].plot(Lat_deg[fp],aWave_norm[fp],color='b')

hf[axes_num].set_xlim([0,50])
hf[axes_num].set_ylim([0,1])
hf[axes_num].grid(alpha=alpha_grid)
hf[axes_num].set_xlabel('Latitude '+r'$\rm\lambda\ (deg)$',fontsize=LabelSize+1)
hf[axes_num].set_ylabel(r'$\rm\omega/\Omega_{He\_local}$',fontsize=LabelSize+1)
hf[axes_num].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
hf[axes_num].yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
hf[axes_num].tick_params(which='major',labelsize=LabelSize,width=TickWidth,length=TickLen,tickdir='in',right=True,top=True)
hf[axes_num].tick_params(which='minor',color='k',width=TickWidth*0.5,length=TickLen*0.5,tickdir='in',right=True,top=True)
hf[axes_num].text(x=0.5,y=1.04,s='wave frequency normalized to '+r'$\rm\Omega_{He\_local}$',ha='center',va='center',fontsize=LabelSize+1,transform=hf[axes_num].transAxes)

hf[axes_num].annotate(text=r'$\rm\lambda=$'+str(round(Lat_deg[fp][-1],2))+r'$\rm\degree$',color='k',fontsize=LabelSize+1,xy=(0.53,0.42),\
    xycoords='axes fraction',xytext=(0.65,0.52),va='center',textcoords='axes fraction',arrowprops=dict(color='k',headwidth=9,width=1))

hf[axes_num].annotate(text='crossover frequency',color='k',fontsize=LabelSize+1,xy=(0.14,0.4),\
    xycoords='axes fraction',xytext=(0.24,0.35),va='center',textcoords='axes fraction',arrowprops=dict(color='k',headwidth=9,width=1))

hf[axes_num].text(x=0.36,y=0.87,s=r'$\rm\omega$='+str(w_wave_norm)+r'$\rm\Omega_{He\_eq}$',color='b',ha='center',va='center',fontsize=LabelSize+1,transform=hf[axes_num].transAxes)

#sys.exit()
    
    
    
#panel 2, plot resonance latitude
    
#input parameters
emic_by=np.array([3])*1e-9
Ek_point=[50] # keV
aeq_point=[60]
plot=True
save_fig=True


theta0_deg=1e-3 # wave normal angle
#density_ratio=15
#L_shell=4 #L shell of simulation
m_res=1 #WPI resonance number (0=Landau resonance)
vel_sign=-1
#ion_com=np.array([77,20,3])/100

particle_mass=env.const.mH
particle_charge=env.const.qi



reslat_colorbar_ylim=[0,30]
reslat_levels=np.linspace(reslat_colorbar_ylim[0],reslat_colorbar_ylim[1],1000)
reslat_colorbar_ytk=[0,10,20,30]
reslat_colorbar_ytklb=['0',r'$\rm10^\degree$',r'$\rm20^\degree$',r'$\rm30^\degree$']


s_xlim=[0,90]
s_xtk=np.arange(0,91,15)
s_xtklb=s_xtk
s_ylim=np.log10([1,500])
s_ytk=np.log10([1,10,100])
s_ytklb=['1','10','100']
s_colorbar_ylim=np.log10([1e-2,1e2+0.01])
s_colorbar_ylim=np.log10([0.1,100+0.01])

#s_colorbar_ylim=np.log10([0.1,30])

#s_colorbar_ytk=np.log10([0.01,0.1,1,10,100])
s_colorbar_ytk=[0.01,0.1,1,10,100]
#s_colorbar_ytk=[0.1,1,10,30]
s_colorbar_ytk=[0.1,1,10,100]


s_colorbar_ytklb=['0.1','1','10','100']
s_colorbar_yminortk=np.log10(np.hstack((np.arange(2*1e-2,1e-1,1e-2),np.arange(2*1e-1,1e0,1e-1),np.arange(2*1e0,1e1,1e0),np.arange(2*1e1,1e2,1e1))))


s_levels=(np.linspace(s_colorbar_ylim[0],s_colorbar_ylim[1],1000))

s_levels=np.log10([0.5,1,5,10])
s_label=['0.5','1.0','5','10']

reslat_contour_levels=[10,20]
reslat_contour_label=['10','20']


#s_yminortk=emic_B_colorbar_ytklb=[r'$10^{-3}$',r'$10^{-2}$',r'$10^{-1}$',r'$10^{0}$',r'$10^{1}$',r'$10^{2}$']
s_yminortk=np.log10(np.hstack((np.arange(2*1e0,1e1,1e0),np.arange(2*1e1,1e2,1e1),np.arange(2*1e2,5*1e2+1,1e2))))


# 1. Define simulation parameters
# Here we define all the initial parameters of the simulation in respect with the particle and the wave

### Simulation parameters


lamdaeq=np.deg2rad(0) #equatorial magnetic latitude (lamda=0 deg)
Beq =env.Bmag_dipole(L_shell,lamdaeq) #equatorial magnetic field strength


wce_eq=env.omega_cyclotron(Beq,env.const.qe,env.const.me) #equatorial electron cyclotron frequency
wcHe_eq=env.omega_cyclotron(Beq,env.const.qe,env.const.mHe) #equatorial helium cyclotron frequency
wcH_eq=env.omega_cyclotron(Beq,env.const.qe,env.const.mH) #equatorial hydrogen cyclotron frequency
wcO_eq=env.omega_cyclotron(Beq,env.const.qe,env.const.mO) #equatorial oxygen cyclotron frequency

wpe_eq=density_ratio*wce_eq #equatorial plasma frequency

ne_eq=(env.const.me*env.const.epsilon0*wpe_eq*wpe_eq)/(env.const.qe*env.const.qe) #calculate equatorial electron density from equatorial plasma frequency

nH_eq=ion_com[0]*ne_eq #77% hydrogen
nHe_eq=ion_com[1]*ne_eq #20% helium
nO_eq=ion_com[2]*ne_eq #3% oxygen


theta0_rad=np.deg2rad(theta0_deg) #convert wave normal angle to rad

    
#import refr_index_full_test


aEk0=np.hstack((10**np.linspace(0,1,50),10**np.linspace(1,2,50),10**np.linspace(2,3,50)))

aeq0_deg=np.linspace(0.1,89.99,200)
if not plot:
    
    aEk0=Ek_point.copy()
    aeq0_deg=aeq_point.copy()


#for iorder in range(len(res_order)):

S_out=np.ones((len(emic_by),len(aEk0),len(aeq0_deg)))*np.nan
mirlat_out=S_out.copy()
reslat_out=S_out.copy()

detadt_out=S_out.copy()
aeq_out=S_out.copy()
for iPar in range(len(emic_by)):
    Byw0_sim=emic_by[iPar]  # y-component of the wave magnetic field, unit: nT
    
    
    w_wave=w_wave_norm*wcHe_eq
    lat_dip=np.arange(0,89.9,0.01)
    aRef_plus=np.ones(len(lat_dip))*np.nan  #refractive index
    aRef_minus=aRef_plus.copy()
    aRef_left=aRef_plus.copy()
    aRef_right=aRef_plus.copy()

    aLat=aRef_plus.copy() 
    plr_plus=aRef_plus.copy()
    plr_minus=aRef_plus.copy()
    aD_stix=aRef_plus.copy()

    
    #in the following, find the cut-off frequency for left-hand EMIC waves
    for i in range(len(lat_dip)):        
        aLat_deg=lat_dip[i]
        aLat_rad=np.deg2rad(aLat_deg)
        aBfield=env.Bmag_dipole(L_shell,aLat_rad)
        S0,D0,P0,R0,L0=wave.stix_parameters(w_wave, ne_eq, nH_eq, nHe_eq, nO_eq, aBfield)
        eta_sq_plus0,eta_sq_minus0,mu0,kappa0,kappaz0,kappax0=wave.refr_index_full(theta0_rad,w_wave,S0,P0,R0,L0)
        plr_plus[i]=-(eta_sq_plus0-S0)/D0
        plr_minus[i]=-(eta_sq_minus0-S0)/D0
        ref_left0,ref_right0,eta_sq_plus0,eta_sq_minus0,kappa0,kappaz0,kappax0=refr_index_full_left.refr_index_full_left(theta0_rad,w_wave,S0,D0,P0,R0,L0)

        aRef_left[i]=ref_left0
        aRef_right[i]=ref_right0
        aD_stix[i]=D0

    fp_index=np.where(np.isnan(aRef_left))
    lat_cutoff_deg=lat_dip[fp_index[0][0]-1]
    bD_stix=np.array([aD_stix[i]*aD_stix[i+1] for i in range(len(aD_stix)-1)])
    fp_D=np.where(bD_stix<0)
    lat_crossover_deg=lat_dip[fp_D[0][0]]  #this is the latitude of crossover frequency
    #sys.exit()
    
    for iEkev in range(len(aEk0)):
        Ekev0=aEk0[iEkev]
        #sys.exit()

        for iAeq0 in range(len(aeq0_deg)):
            aeq0_rad=np.deg2rad(aeq0_deg[iAeq0])
            print('EMIC Bwy='+str(round(Byw0_sim*1e9,2))+' nT')
            print('Initial energy: '+str(np.round(Ekev0,2))+' keV, '+'aeq0='+str(round(np.rad2deg(aeq0_rad),1))+' degree')

            aeq_rad=aeq0_rad
            coef=np.array([1,0,0,0,0,3*(np.sin(aeq_rad))**4,-4*(np.sin(aeq_rad))**4])
            result=np.roots(coef)
            for i in range(len(result)):
                if result[i].real>=0 and result[i].imag==0:
                    kappa=result[i].real
                    aMlat_rad=np.arccos(np.sqrt(kappa))
                    aMlat_deg=np.rad2deg(aMlat_rad)
            #print('latitude: '+str(aMlat_deg))
            lamda0_mirr_deg=aMlat_deg #latitude of the mirror point
            #in the following, find the latitude of linear resonance
            #lamda0_end_deg=min(lamda0_mirr_deg-0.0001,lat_cutoff_deg)         
            lamda0_end_deg=min(lamda0_mirr_deg-0.0001,lat_crossover_deg)            

            bLat_deg=np.arange(lamda0_end_deg,0,-0.01)
            if len(bLat_deg)<=1: continue

            bResCon=np.zeros(len(bLat_deg))*np.nan
            cResCon=bResCon.copy()
            for iLamda in range(len(bLat_deg)):
                bLat_rad=np.deg2rad(bLat_deg[iLamda])
                bAlpha=env.aeq2alpha(L_shell,bLat_rad,aeq0_rad)
                bBfield=env.Bmag_dipole(L_shell,bLat_rad)
                bUpar,bUper,bPpar,bPper,bGamma=env.initial_velocity(Ekev0,bAlpha,particle_mass)
    
                #calculate magnetic field strength at ion's initial location using WPIT.Environment_mod.Bmag_dipole routine
                #B0 =env.Bmag_dipole(L_shell,lamda0)
                #calculate electron gyrofrequency at ion's initial location using WPIT.Environment_mod.omega_cyclotron routine
                #bWce=env.omega_cyclotron(bBfield,env.const.qe,env.const.me)
                #calculate helium gyrofrequency at ion's initial location using WPIT.Environment_mod.omega_cyclotron routine
                #bWcHe=env.omega_cyclotron(bBfield,env.const.qe,env.const.mHe)
                #calculate hydrogen gyrofrequency at ion's initial location using WPIT.Environment_mod.omega_cyclotron routine
                bWcyc=env.omega_cyclotron(bBfield,particle_charge,particle_mass)
                #calculate oxygen gyrofrequency at ion's initial location using WPIT.Environment_mod.omega_cyclotron routine
                #bWcO=env.omega_cyclotron(bBfield,env.const.qe,env.const.mO)
                
                #calculate the Stix parameters at ion's initial location using WPIT.WaveProperties_mod.stix_parameters routine
                bS0,bD0,bP0,bR0,bL0=wave.stix_parameters(w_wave, ne_eq, nH_eq, nHe_eq, nO_eq, bBfield)
                #calculate the refractive index and the wavenumber at ion's initial location 
                           
                bRef_left,bRef_right,bEta_sq_plus,bEta_sq_minus,bKappa,bKappaz,bKappax=refr_index_full_left.refr_index_full_left(theta0_rad,w_wave,bS0,bD0,bP0,bR0,bL0)

                
                bResCon[iLamda]=m_res*bWcyc/bGamma+bKappaz*((vel_sign)*bPpar)/(bGamma*particle_mass)-w_wave
            #sys.exit()
            #sys.exit()
            bResCon_sq=[bResCon[i]*bResCon[i+1] for i in range(len(bResCon)-1)]
            cResCon[0:-1]=bResCon_sq
            fp=np.where(cResCon<0)
            # if len(fp[0])<1:
            #     continue
            # else:
            #     sys.exit()
            if len(fp[0])>0:
                lamda0_rad=np.deg2rad(bLat_deg[fp[0][-1]])  #latitude of the linear resonance point
                #2. Find initial proton's local pitch angle
                
                #using WPIT.Environment_mod.aeq2alpha routine
                alpha0=env.aeq2alpha(L_shell,lamda0_rad,aeq0_rad)
                
                #print('\u03B1:',np.rad2deg(alpha0))
                #3. Find initial momentum, velocity and lorentz factor
                
                #using WPIT.Environment_mod.initial_velocity routine
                upar0,uper0,ppar0,pper0,gamma0=env.initial_velocity(Ekev0,alpha0,particle_mass)
                ppar0=ppar0*(-1)
                
                # print('upar0:',upar0,'m/s')
                # print('uper0:',uper0,'m/s')
                # print('ppar0:',ppar0,'Ns')
                # print('pper0:',pper0,'Ns')
                # print('gamma0:',gamma0)
                
                #4. Calculate all the initial parameters
                
                #calculate magnetic field strength at ion's initial location using WPIT.Environment_mod.Bmag_dipole routine
                B0 =env.Bmag_dipole(L_shell,lamda0_rad)
                #calculate electron gyrofrequency at ion's initial location using WPIT.Environment_mod.omega_cyclotron routine
                #wce_0=env.omega_cyclotron(B0,env.const.qe,env.const.me)
                #calculate helium gyrofrequency at ion's initial location using WPIT.Environment_mod.omega_cyclotron routine
                #wcHe_0=env.omega_cyclotron(B0,env.const.qe,env.const.mHe)
                wcyc_0=env.omega_cyclotron(B0,particle_charge,particle_mass)

                #calculate hydrogen gyrofrequency at ion's initial location using WPIT.Environment_mod.omega_cyclotron routine
                #wcH_0=env.omega_cyclotron(B0,env.const.qe,env.const.mH)
                #calculate oxygen gyrofrequency at ion's initial location using WPIT.Environment_mod.omega_cyclotron routine
                #wcO_0=env.omega_cyclotron(B0,env.const.qe,env.const.mO)
                
                #calculate the Stix parameters at ion's initial location using WPIT.WaveProperties_mod.stix_parameters routine
                S0,D0,P0,R0,L0=wave.stix_parameters(w_wave, ne_eq, nH_eq, nHe_eq, nO_eq, B0)
                
                #calculate the refractive index and the wavenumber at ion's initial location 
                # using WPIT.WaveProperties_mod.refr_index_full routine
                eta_sq_plus0,eta_sq_minus0,mu0,kappa0,kappaz0,kappax0=wave.refr_index_full(theta0_rad,w_wave,S0,P0,R0,L0)
                
                ref_left0,ref_right0,eta_sq_plus0,eta_sq_minus0,kappa0,kappaz0,kappax0=refr_index_full_left.refr_index_full_left(theta0_rad,w_wave,S0,D0,P0,R0,L0)


                #calculate resonant velocity and energy at ion's initial location using 
                # WPIT.WaveProperties_mod.resonant_velocity routine
                v_para_res0, v_per_res0, v_tot_res0, E_res0,gamma_res0=wave.resonant_velocity(m_res,w_wave,kappaz0,wcyc_0,alpha0,particle_mass)
                
                #calculate the gradient of hydrogen gyrofrequency at ion's initial location using WPIT.Environment_mod.dwc_ds routine
                dwcds0=env.dwc_ds(wcyc_0,lamda0_rad,L_shell)
                #calculate the gradient of the magnetic field strength at ion's initial location using WPIT.Environment_mod.dB_ds routine
                dBdz0=env.dB_ds(B0,lamda0_rad,L_shell)
                #calculate the wave component amplitudes at ion's initial location 
                # using WPIT.WaveProperties_mod.wave_amplitudes_bell routine
                Bxw0, Byw0, Bzw0, Exw0, Eyw0, Ezw0=wave.wave_amplitudes_bell(ref_left0,P0,D0,S0,Byw0_sim,theta0_rad)
                #calculate wpi parameters at electron's initial location using WPIT.WPI_mod.EMIC_ion_mod.wpi_params routine
                beta0,BwR0,BwL0,EwR0,EwL0,pwR0,pwL0,wR0,wL0=wpi.wpi_params(pper0,kappax0,particle_charge,particle_mass,B0,Exw0,Eyw0,Bxw0,Byw0,gamma0)
                beta0_sim=0
                
                #calulcate electorn plasma frequency
                #wpe_0=env.omega_plasma(ne_eq,env.const.qe,env.const.me)
                #calulcate helium plasma frequency
                #wpHe_0=env.omega_plasma(nHe_eq,env.const.qe,env.const.mHe)
                #calulcate hydrogen plasma frequency
                #wpH_0=env.omega_plasma(nH_eq,env.const.qe,env.const.mH)
                #calulcate oxygen plasma frequency
                #wpO_0=env.omega_plasma(nO_eq,env.const.qe,env.const.mO)
                dwcds=env.dwc_ds(wcyc_0,lamda0_rad,L_shell)
                
                #calculate initial parameters for the investigation of non-linear interactions
                C0_0=wpi.nonlinear_C0(ppar0,kappaz0,m_res,gamma0,particle_charge,particle_mass,wcyc_0,Ezw0)
                C1p_0=wpi.nonlinear_C1p(pper0,ppar0,kappaz0,m_res,particle_charge,particle_mass,gamma0,wR0,EwR0,wcyc_0)
                C1m_0=wpi.nonlinear_C1m(pper0,ppar0,kappaz0,m_res,particle_charge,particle_mass,gamma0,wL0,EwL0,wcyc_0)
                thet_0,wtrsq_0=wpi.nonlinear_theta(C0_0,C1p_0,C1m_0,m_res,beta0)
                dkpar_dt0=0
                #revised by Su Zhou
                lamda0_rad_2=lamda0_rad+0.00001
                #2. Find initial electron's local pitch angle
                
                #using WPIT.Environment_mod.aeq2alpha routine
                alpha0_2=env.aeq2alpha(L_shell,lamda0_rad_2,aeq0_rad)
                
                #3. Find initial momentum, velocity and lorentz factor
                
                #using WPIT.Environment_mod.initial_velocity routine
                upar0_2,uper0_2,ppar0_2,pper0_2,gamma0_2=env.initial_velocity(Ekev0,alpha0_2,particle_mass)
                
                # print('upar0:',upar0,'m/s')
                # print('uper0:',uper0,'m/s')
                # print('ppar0:',ppar0,'Ns')
                # print('pper0:',pper0,'Ns')
                # print('gamma0:',gamma0)
                
                
                #4. Calculate all the initial parameters
                
                #calculate magnetic field strength at ion's initial location using WPIT.Environment_mod.Bmag_dipole routine
                B0_2 =env.Bmag_dipole(L_shell,lamda0_rad_2)
                #calculate electron gyrofrequency at ion's initial location using WPIT.Environment_mod.omega_cyclotron routine
                #wce_0_2=env.omega_cyclotron(B0_2,env.const.qe,env.const.me)
                #calculate helium gyrofrequency at ion's initial location using WPIT.Environment_mod.omega_cyclotron routine
                #wcHe_0_2=env.omega_cyclotron(B0_2,env.const.qe,env.const.mHe)
                #calculate hydrogen gyrofrequency at ion's initial location using WPIT.Environment_mod.omega_cyclotron routine
                #wcH_0_2=env.omega_cyclotron(B0_2,env.const.qe,env.const.mH)
                #calculate oxygen gyrofrequency at ion's initial location using WPIT.Environment_mod.omega_cyclotron routine
                #wcO_0_2=env.omega_cyclotron(B0_2,env.const.qe,env.const.mO)
                
                #calculate the Stix parameters at ion's initial location using WPIT.WaveProperties_mod.stix_parameters routine
                S0_2,D0_2,P0_2,R0_2,L0_2=wave.stix_parameters(w_wave, ne_eq, nH_eq, nHe_eq, nO_eq, B0_2)
                #calculate the refractive index and the wavenumber at ion's initial location 
                # using WPIT.WaveProperties_mod.refr_index_full routine
                
                ref_left0_2,ref_right0_2,eta_sq_plus0_2,eta_sq_minus0_2,kappa0_2,kappaz0_2,kappax0_2=refr_index_full_left.refr_index_full_left(theta0_rad,w_wave,S0_2,D0_2,P0_2,R0_2,L0_2)

                
                dkdl=(kappaz0_2-kappaz0)/(lamda0_rad_2-lamda0_rad)
                dkdz=dkdl*1/(env.const.Re*L_shell*np.cos(lamda0_rad)*np.sqrt(1+3*np.sin(lamda0_rad)**2))
                dkpardt0=(ppar0/(gamma0*particle_mass))*np.cos(theta0_rad)*dkdz
                dkpar_dt0=dkpardt0
    
                H_0=wpi.nonlinear_H(pper0,ppar0,kappaz0,gamma0,m_res,particle_mass,wcyc_0,dkpar_dt0,dwcds,0)
                S_0=wpi.nonlinear_S(H_0,wtrsq_0)
                
                deta_dt0=wpi.detadt(-ppar0,m_res,wcyc_0,gamma0,kappaz0,particle_mass,w_wave)
                
                #outputs
                S_out[iPar,iEkev,iAeq0]=S_0
                
                print('|S|='+str(round(abs(S_0),3)))
                print('-----------------------')
                
                #sys.exit()
                
                reslat_out[iPar,iEkev,iAeq0]=np.rad2deg(lamda0_rad)
                mirlat_out[iPar,iEkev,iAeq0]=lamda0_mirr_deg
                aeq_out[iPar,iEkev,iAeq0]=np.rad2deg(aeq0_rad)
                    
    
        
        
    # Plots of the results
    if not plot:
        continue

    # plot resonance latitude
    axes_num=int(iPar%num_column)+2
    
    reslat_05_plot=hf[axes_num].contour(aeq0_deg,np.log10(aEk0),(reslat_out[iPar,:,:]),levels=reslat_contour_levels,colors='k',linestyles='dashed')
    fmt = {}
    strs = ['1.0']
    for l, s in zip(reslat_05_plot.levels, reslat_contour_label):
        fmt[l] = s
    manual = [(45,2.5),(63,2.5)]

    hf[axes_num].clabel(reslat_05_plot, reslat_05_plot.levels, inline=True, fmt=fmt, fontsize=LabelSize,manual=manual)


    

    reslat_05_plot=hf[axes_num].contourf(aeq0_deg,np.log10(aEk0),(reslat_out[iPar,:,:]),levels=reslat_levels,vmin=reslat_colorbar_ylim[0],vmax=reslat_colorbar_ylim[1],\
                                        cmap='jet',ls='None')
    #hf[axes_num].plot(aeq_point[0],np.log10(Ek_point[0]),marker='^',color='k',markerfacecolor='none',ms=MarkerSize)
    
    hf[axes_num].grid(alpha=alpha_grid)    
    hf[axes_num].set_xlim(s_xlim)
    hf[axes_num].set_xticks(s_xtk,labels=s_xtklb)

    hf[axes_num].set_ylim(s_ylim)
    hf[axes_num].set_yticks(s_ytk,labels=s_ytklb)
    hf[axes_num].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(6))

    hf[axes_num].yaxis.set_minor_locator(matplotlib.ticker.FixedLocator(s_yminortk))
    hf[axes_num].tick_params(which='major',labelsize=LabelSize,width=TickWidth,length=TickLen,tickdir='in',right=True,top=True)
    hf[axes_num].tick_params(which='minor',color='k',width=TickWidth*0.5,length=TickLen*0.5,tickdir='in',right=True,top=True)

    hf[axes_num].set_xlabel(r'$\rm\alpha_{eq}$ (deg)',fontsize=LabelSize+1)
    hf[axes_num].set_ylabel(r'$\rmE_k$ (keV)',fontsize=LabelSize+1)
    #add a colorbar
    pos=hf[axes_num].get_position()
    bar_axes=myfig.add_axes([0.1,0.1,0.1,0.1])
    bar_position=[0.95,pos.y0,0.01,(pos.y1-pos.y0)]
    reslat_Colorbar=plt.colorbar(mappable=reslat_05_plot,cax=bar_axes)
    reslat_Colorbar.ax.set_position(bar_position)
    bar_axes.yaxis.set_ticks_position('right')
    reslat_Colorbar.ax.yaxis.set_ticks(reslat_colorbar_ytk,labels=reslat_colorbar_ytklb,minor=False)
    reslat_Colorbar.ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())

    reslat_Colorbar.ax.tick_params(which='major',labelsize=LabelSize,labelcolor='k',color='k',width=TickWidth,length=TickLen,tickdir='out')
    reslat_Colorbar.ax.tick_params(which='minor',color='k',width=TickWidth*0.5,length=TickLen*0.5)
    #hf[axes_num].text(x=0.5,y=1.04,s='f='+str(aW_ratio)+r'$\rm\Omega_{H^+}$',\
    #                  ha='center',va='center',fontsize=LabelSize,transform=hf[axes_num].transAxes)
    #hf[axes_num].text(x=0.5,y=1.04,s=r'$\rmB_{w}^{y}$='+str(round(Byw0_sim*1e9,1))+' nT',ha='center',va='center',fontsize=LabelSize+2,transform=hf[axes_num].transAxes)
    hf[axes_num].text(x=0.5,y=1.04,s=r'$\rm\omega$='+str(w_wave_norm)+r'$\rm\Omega_{He\_eq}$'+', '+r'$\rm\psi$='+str(int(theta0_deg))+r'$\rm\degree$',\
                  ha='center',va='center',fontsize=LabelSize+1,transform=hf[axes_num].transAxes)
        
    hf[axes_num].text(x=1.19,y=0.5,s=r'$\rm\lambda_{res}\ (deg)$',ha='center',va='center',rotation=90,fontsize=LabelSize,transform=hf[axes_num].transAxes)

    #sys.exit()
    
    
    
    
    
    
    
    
    
    
    
for i in range(len(hf)):
    hf[i].text(x=0.05,y=0.96,s=panel_label[i],color='k',fontsize=LabelSize+3,weight='bold',ha='center',va='center',transform=hf[i].transAxes)

save_name='EMIC_dispersion_resLat'
if save_fig:
    myfig.savefig(save_dir+save_name+'.png')


    





   
    

    

    









