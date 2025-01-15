# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 08:25:06 2024

@author: zhousu
"""

import numpy as np
import matplotlib.pyplot as plt
import numlib,sys,matplotlib
import scipy.io as sio

read_dir=r'E:\mywork\my_drafts\nonlinear_affected_by_EMIC_parameters\data/'
read_name_0='test_particle_diffusion_coef_012.mat'
read_name_1='test_particle_diffusion_coef_345.mat'
read_name_2='test_particle_diffusion_coef_678.mat'
save_dir=r'E:\mywork\my_drafts\nonlinear_affected_by_EMIC_parameters\figs/'
save_fig=True

LineWidth=1;LabelSize=15;TickWidth=1.6;TickLen=6
alpha_grid=0.4
aeq_xlim=[6,90]
adv_aeq_ylim=[-0.5,0.5]
diff_aeq_ylim=[1e-5,1e1]
panel_label=['(a)','(c)','(e)',\
             '(b)','(d)','(f)']


plt.close('all')
num_row=2
num_column=3
myfig=plt.figure(figsize=(19,12))
hf=numlib.axes_create(row=num_row, column=num_column,top=0.05,bottom=0.06,left=0.05,right=0.04,gap_row=0.05,gap_column=0.06)


for iFile in range(3):
    if iFile==0:
        read_name=read_name_0
    elif iFile==1:
        read_name=read_name_1
    else:
        read_name=read_name_2



    data=sio.loadmat(read_dir+read_name)
    
    
    aeq_start=data['aeq_start']
    aeq_end=data['aeq_end']
    ek_start=data['ek_start']
    ek_end=data['ek_end']
    
    
    
    
    emic_intensity=data['emic_intensity'][0]
    Energy=data['Energy'][0]
    eta0_deg=data['eta0_deg'][0]
    freq_norm=data['freq_norm'][0]
    par_number=data['par_number'][0]
    
    pitch_angle=data['pitch_angle'][0]
    transmit_time=data['transmit_time']
    wna=data['wna'][0]
    advection_aeq=np.zeros((len(emic_intensity),len(Energy),len(pitch_angle)))*np.nan
    advection_ek=advection_aeq.copy()
    
    diff_aeq=advection_aeq.copy()
    diff_ek=advection_aeq.copy()
    
    
    color_lb=['r','g','b']
    
    #delta_aeq=data['delta_aeq']
    #delta_Ek=data['delta_Ek']
    
    for iPar in par_number:
        aEmic_nT=emic_intensity[iPar]*1e9
        aFreq_norm=freq_norm[iPar]
        aWna_deg=wna[iPar]
        for iEk in range(len(Energy)):
            for iAeq in range(len(pitch_angle)):
                aeq0_deg=pitch_angle[iAeq]
                aTime=transmit_time[iPar,iEk,iAeq,:]
    
                aeq_start_deg=np.rad2deg(aeq_start[iPar,iEk,iAeq,:])
                aeq_end_deg=np.rad2deg(aeq_end[iPar,iEk,iAeq,:])
                aDaeq_deg=aeq_end_deg-aeq_start_deg
                
                aEk_start=ek_start[iPar,iEk,iAeq,:]
                aEk_end=ek_end[iPar,iEk,iAeq,:]
                aDek=aEk_end-aEk_start
    
                
                
                
                
    
                #aDaeq=np.rad2deg(aeq_start[iPar,iEk,iAeq,:]-aeq_end[iPar,iEk,iAeq,:])
                #aDek=(ek_start[iPar,iEk,iAeq,:]-ek_end[iPar,iEk,iAeq,:])
    
                #aDaeq=np.rad2deg(delta_aeq[iPar,iEk,iAeq,:])
                #aDek=delta_ek[iPar,iEk,iAeq,:]
                bDaeq_dt=aDaeq_deg/aTime
                bDek_dt=aDek/aTime
    
    
                # if aeq0_deg<10:
                    
                #     fp=(abs(aDaeq)<5)
                #     bDaeq_dt=aDaeq[fp]/aTime[fp]
    
                #bAdvection_aeq=np.mean(bDaeq_dt)
                #bAdvection_ek=np.mean(bDek_dt)
    
                #sys.exit()
                #if bAdvection_aeq>0.1:
                #    sys.exit()
                #bAdvection_aeq=np.mean(aDaeq)
    
                
                #bAdvection_ek=np.mean(bDek_dt)
                advection_aeq[iPar,iEk,iAeq]=np.mean(bDaeq_dt)
                advection_ek[iPar,iEk,iAeq]=np.mean(bDek_dt)
                #cal diffusion coef
                #aeqTP=aeq_end_deg
                #aeqTP_mean=np.mean(aeqTP)
                #print(aeqTP_mean)
                cAeq_diff=(aeq_end_deg-np.mean(aeq_end_deg))**2/(2*aTime)
                cEk_diff=(aEk_end-np.mean(aEk_end))**2/(2*aTime)
    
                
                
                diff_aeq[iPar,iEk,iAeq]=np.mean(cAeq_diff)
                diff_ek[iPar,iEk,iAeq]=np.mean(cEk_diff)
    
                #sys.exit()
                
                
                
            axes_num=iFile
            if axes_num==0:
                
                hf[axes_num].plot(pitch_angle,advection_aeq[iPar,iEk,:],color=color_lb[np.mod(iPar,3)],label=r'$\rmB_w^y=$'+str(round(aEmic_nT,1))+' nT')
            elif axes_num==1:
                
                hf[axes_num].plot(pitch_angle,advection_aeq[iPar,iEk,:],color=color_lb[np.mod(iPar,3)],label=r'$\rm\omega=$'+str(aFreq_norm)+r'$\rm\Omega_{He\_eq}$')
            elif axes_num==2:
                hf[axes_num].plot(pitch_angle,advection_aeq[iPar,iEk,:],color=color_lb[np.mod(iPar,3)],label=r'$\rm\psi=$'+str(int(aWna_deg))+r'$\rm^\degree$')

                
            

            hf[axes_num].set_xlim(aeq_xlim)
            hf[axes_num].set_ylim(adv_aeq_ylim)
            hf[axes_num].grid(alpha=alpha_grid)
            hf[axes_num].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
            hf[axes_num].yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
            hf[axes_num].tick_params(which='major',labelsize=LabelSize-1,width=TickWidth,length=TickLen,tickdir='in',right=True,top=True)
            hf[axes_num].tick_params(which='minor',width=TickWidth*0.6,length=TickLen*0.5,tickdir='in',right=True,top=True)
            hf[axes_num].set_xlabel(r'$\rm\alpha_{eq}\ (deg)$',fontsize=LabelSize+3)
            hf[axes_num].set_ylabel(r'$\rm<A_{\alpha_{eq}}>\ (s^{-1})$',fontsize=LabelSize+3)
            
            if iPar==0:
                hf[axes_num].text(x=0.5,y=1.04,s=r'$\rm\omega=$'+str(aFreq_norm)+r'$\rm\Omega_{He^+}$'+', '+r'$\rm\psi=$'+str(int(aWna_deg))+r'$\rm^\degree$',
                              fontsize=LabelSize+3,ha='center',transform=hf[axes_num].transAxes)
            elif iPar==3:
                hf[axes_num].text(x=0.5,y=1.04,s=r'$\rmB_w^y=$'+str(round(aEmic_nT,1))+' nT'+', '+r'$\rm\psi=$'+str(int(aWna_deg))+r'$\rm^\degree$',
                              fontsize=LabelSize+3,ha='center',transform=hf[axes_num].transAxes)
            elif iPar==6:
                hf[axes_num].text(x=0.5,y=1.04,s=r'$\rmB_w^y=$'+str(round(aEmic_nT,1))+' nT, '+r'$\rm\omega=$'+str(aFreq_norm)+r'$\rm\Omega_{He\_eq}$',
                              fontsize=LabelSize+3,ha='center',transform=hf[axes_num].transAxes)
            
            
            axes_num=iFile+3
            if axes_num==3:
                hf[axes_num].plot(pitch_angle,diff_aeq[iPar,iEk,:],color=color_lb[np.mod(iPar,3)],label=r'$\rmB_w^y=$'+str(round(aEmic_nT,1))+' nT')
            elif axes_num==4:
                hf[axes_num].plot(pitch_angle,diff_aeq[iPar,iEk,:],color=color_lb[np.mod(iPar,3)],label=r'$\rm\omega=$'+str(aFreq_norm)+r'$\rm\Omega_{He\_eq}$')
            elif axes_num==5:
                hf[axes_num].plot(pitch_angle,diff_aeq[iPar,iEk,:],color=color_lb[np.mod(iPar,3)],label=r'$\rm\psi=$'+str(int(aWna_deg))+r'$\rm^\degree$')


            hf[axes_num].set_xlim(aeq_xlim)
    
            hf[axes_num].set_ylim(diff_aeq_ylim)
            hf[axes_num].set_yscale('log')
            hf[axes_num].grid(alpha=alpha_grid)
            hf[axes_num].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
            hf[axes_num].tick_params(which='major',labelsize=LabelSize-1,width=TickWidth,length=TickLen,tickdir='in',right=True,top=True)
            hf[axes_num].tick_params(which='minor',width=TickWidth*0.6,length=TickLen*0.5,tickdir='in',right=True,top=True)
            hf[axes_num].set_xlabel(r'$\rm\alpha_{eq}\ (deg)$',fontsize=LabelSize+3)
            hf[axes_num].set_ylabel(r'$\rm<D_{\alpha_{eq}\alpha_{eq}}>\ (s^{-1})$',fontsize=LabelSize+3)
            
            
            #sys.exit()
    axes_num=iFile  
    legend_s = hf[axes_num].legend(handletextpad=0, handlelength=0,loc='lower left',bbox_to_anchor=(0.01,0.01),fontsize=LabelSize)
    for n, text in enumerate( legend_s.texts ):
        text.set_color( color_lb[n] )
    
    axes_num=iFile+3
    legend_s = hf[axes_num].legend(handletextpad=0, handlelength=0,loc='lower left',bbox_to_anchor=(0.01,0.01),fontsize=LabelSize)
    for n, text in enumerate( legend_s.texts ):
        text.set_color( color_lb[n] )
    #sys.exit()
for iP in range(len(hf)):
    hf[iP].text(x=0.02,y=0.94,s=panel_label[iP],color='m',weight='bold',bbox=dict(facecolor='w',edgecolor='None',alpha=0),fontsize=LabelSize+4,transform=hf[iP].transAxes)

 


save_name='diff_coef_50keV'
if save_fig is True:
    myfig.savefig(save_dir+save_name+'.png')
    


