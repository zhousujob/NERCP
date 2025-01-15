# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 12:40:57 2020

@author: zhousu
"""

def axes_create(row,column,top=0.05,bottom=0.05,left=0.05,right=0.05,gap_row=0.05,gap_column=0.05):
    import matplotlib.pyplot as plt
    width=(1-left-right-(column-1)*gap_column)/column;# the width of each panel
    height=(1-top-bottom-(row-1)*gap_row)/row; # the height of each panel
    hf=[]
    for f in range(row*column):
        hf.append(plt.subplot(row,column,f+1));
    del f
    for i in range(row):
        for j in range(column):        
            x=left+(j+1-1)*(width+gap_column) #the x position of the panel
            y=bottom+(row-(i+1))*(height+gap_row)# the y position of the panel
            hf[(i+1-1)*column+j].set_position([x,y,width,height])
    return hf

def num2str(number,length):
    #  'number' is the input number
    # 'length' is the length of the output string
    if len(str(number))>length:
        raise Exception('variable "length" must be longer than "number"')
    aStr=str(number)
    for i in range(length):
        if len(aStr)<length:
            aStr='0'+aStr
    return aStr;

def hhmmss(decimal_hour):
    import math
    if (decimal_hour<0) or (decimal_hour>=24):
        raise Exception('the input time is not in 0-23.9999 range')
    (aDecimal,hh)=math.modf(decimal_hour)
    aMinute=(aDecimal)*60;
    (bDecimal,mm)=math.modf(aMinute)
    ss=round((bDecimal)*60);
    (mm,ss)=(mm+1,0) if ss==60 else (mm,ss)
    (hh,mm)=(hh+1,0) if mm==60 else (hh,mm)
    hh=0 if hh==24 else hh
    hh='0'+str(int(hh)) if hh<10 else str(int(hh))
    mm='0'+str(int(mm)) if mm<10 else str(int(mm))
    ss='0'+str(int(ss)) if ss<10 else str(int(ss))
    return hh,mm,ss  


def background_ssusi(axObj,hemis,sector,mlat_pos,mlat_color,mlat_min=50,mlat_step=10,\
                     font_size=12,line_style='-',line_width=0.5,line_color='black',mlt_dis=4,show_mlt=True,show_mlat=True):
    """
    # axObj: the object of the axes of the panel
    # hemis: hemisphere of the background,'NH' for northern hemisphere,'SH' for southern hemisphere
    # sector: sector of the background,'all_mlt','dayside','nightside','dawnside','duskside'
    # mlat_pos: mlat position
    # mlat_color: mlat color
    # mlat_min: the mlat of the boundary, 45, 50 or 55 degrees MLAT
    # font_size: the fontsize of the mlat and mlt
    # line_style: the style of the line
    # line_width: the width of the line
    # line_color: the color of the line
    # mlt_dis: the radius distance of the MLT text from the boundary
    # show_mlt: whether show the MLT text,true for showing,false for not showing
    # show_mlat: whether show the MLAT text, true for showing, falsie for not showing

    """
    
    
    
    

    import matplotlib.pyplot as plt
    import math
    import numpy as np
    
    if sector=='all_mlt':    
       theta=np.linspace(0,2*math.pi+0.01,10000)
       theta_mlt=np.arange(0,2*math.pi+0.01,math.pi/6)
    elif sector=='dayside':
       theta=np.linspace(0,math.pi+0.01,1e4)
       theta_mlt=np.arange(0,math.pi+0.01,math.pi/6)       
    elif sector=='nightside':
       theta=np.linspace(-math.pi,0+0.01,1e4)
       theta_mlt=np.arange(-math.pi,0+0.01,math.pi/6)       
    elif sector=='dawnside':
       theta=np.linspace(-math.pi/2,math.pi/2+0.01,1e4)
       theta_mlt=np.arange(-math.pi/2,math.pi/2+0.01,math.pi/6)       
    elif sector=='duskside':
       theta=np.linspace(math.pi/2,3*math.pi/2+0.01,1e4)
       theta_mlt=np.arange(math.pi/2,3*math.pi/2+0.01,math.pi/6)       
                    
    if mlat_min==45:
        Lat=np.arange(mlat_min,81,mlat_step)
    elif mlat_min==50:
       Lat=np.arange(mlat_min,81,mlat_step)
    else:
       Lat=np.arange(mlat_min,81,mlat_step)           
    rho=90-Lat
# plot the circles     
    for i in np.arange(len(Lat)):                      
      R_cir=rho[i]*np.ones(len(theta))
      plt.sca(axObj)
      x_cir=R_cir*np.cos(theta);y_cir=R_cir*np.sin(theta)
      plt.plot(x_cir,y_cir,color=line_color,linestyle=line_style,linewidth=line_width)
      xlat=rho[i]*np.cos(mlat_pos);ylat=rho[i]*np.sin(mlat_pos);
      if show_mlat is True:
         if hemis=='NH':
            axObj.text(xlat,ylat,str(Lat[i])+r'$\rm\degree$',color=mlat_color,fontsize=font_size,family='Times New Roman',ha='center',\
                 va='center')
         else:
            axObj.text(xlat,ylat,'-'+str(Lat[i])+r'$\rm\degree$',color=mlat_color,fontsize=font_size,family='Times New Roman',ha='center',\
                 va='center')          
# plot the line for eache MLT         
    for i in theta_mlt:
          x=(90-mlat_min)*np.cos(i)
          y=(90-mlat_min)*np.sin(i)
          axObj.plot([0,x],[0,y],linewidth=line_width,color=line_color,linestyle=line_style)
# give the mlt text string
    if show_mlt is True:
             mlt=(theta_mlt*180/math.pi+90)/15
             fp=mlt>=24
             mlt[fp]=mlt[fp]-24
             for i in np.arange(0,len(mlt)):
                 Rmax=max(rho)
                 aMlt=mlt[i]    
                 x=(mlt_dis+Rmax)*math.cos(theta_mlt[i])
                 y=(mlt_dis+Rmax)*math.sin(theta_mlt[i])
                 bMlt=int(round(aMlt));bMlt=str(bMlt);    
                 if len(bMlt)<2:
                    bMlt='0'+bMlt        
                 axObj.text(x,y,bMlt,color=line_color,fontsize=font_size,ha='center',\
                 va='center')  
    return;

   

    
    
    
    

    
    


def background_superDarn(axObj,boundary=50,width_mlat=20,hemi='north',mlt_text=None,show_mlat=False,show_mlt=False):
    import matplotlib.pyplot as plt
    import math
    import numpy as np
#    plt.close('all')
#    plt.figure(1,figsize=[8,6])
#    ax1=plt.subplot(1,1,1)
#plt.axis('equal')
    plt.sca(axObj)# set axObj to be the current obj
    # draw four lines, up,down, left and right
    max_v=90-boundary
    plt.plot([-max_v,-max_v],[-max_v,max_v],linestyle='-',color='black',linewidth=0.8)
    plt.plot([max_v,max_v],[-max_v,max_v],linestyle='-',color='black',linewidth=0.8)
    plt.plot([-max_v,max_v],[-max_v,-max_v],linestyle='-',color='black',linewidth=0.8)
    plt.plot([-max_v,max_v],[max_v,max_v],linestyle='-',color='black',linewidth=0.8)
    axObj.set_xlim([-max_v,max_v])
    axObj.set_ylim([-max_v,max_v])
    axObj.set_xticks([])
    axObj.set_xticklabels(' ')
    axObj.set_yticks([])
    axObj.set_yticklabels([])
    axObj.axis('scaled')
    # draw 24 meridian lines at 24 MLT
    gap_mlt=1
    mlt=np.arange(0,24,gap_mlt)
    maxLat=80;co_maxLat=90-maxLat;# the maximum latitude of plotting meridian lines
    minLat=20;co_minLat=90-minLat;# the minimum latitude of plotting meridian lines
    for i in mlt:  
        theta=(i*15+270)*math.pi/180
        x1=co_maxLat*np.cos(theta)    
        y1=co_maxLat*np.sin(theta)
        x2=co_minLat*np.cos(theta)
        y2=co_minLat*np.sin(theta)
        plt.plot([x1,x2],[y1,y2],color='grey',linestyle=':',linewidth=0.5)
        if len(list(set([i]).intersection(set(mlt_text))))>0:# showing the MLT            
            x=(max_v+3)*np.cos(theta)
            y=(max_v+3)*np.sin(theta)
            if i<10:
                s_mlt='0'+str(i)
            else:
                s_mlt=str(i)
            axObj.text(x,y,s=s_mlt,ha='center',va='center',fontname='Times New Roman')
        
       
# plot the circles
    for i in np.arange(80,boundary-1,-width_mlat):
        # draw circles
        radius=90-i
        cir=plt.Circle([0,0],radius,color='black',linestyle=':',linewidth=0.5,fill=False)
        axObj.add_patch(cir)
        # show the mlat in the figure
        if show_mlat is True:
            pos_mlt=8.5  # the magnetic local time where latitudes are shown
            theta_show=(pos_mlt*15+270)*math.pi/180
            x=(90-i)*np.cos(theta_show)
            y=(90-i)*np.sin(theta_show)
            if hemi=='north':
               latL=str(i)
            else:
               latL='-'+str(i)
            axObj.text(x,y,s=latL+'$^o$',ha='center',va='center',fontname='Times New Roman')
 
def rbsp_gse2mgse(var_input,wgse):
    import numpy as np
    
#MagData=(aEpoch,bMag)
    if np.shape(var_input[1])!=np.shape(wgse):
        raise ValueError('the dimensions of var_input[1] and wgse are not consistent ')
    aEpoch=var_input[0]
    aVar=var_input[1]
    zgse = [0,0,1]

    aWgse=np.zeros((len(aEpoch),3))*np.nan
    datx = np.zeros(len(aEpoch))*np.nan
    daty = np.zeros(len(aEpoch))*np.nan
    datz = np.zeros(len(aEpoch))*np.nan
    
    for j in range(len(aEpoch)):
        aWgse[j,:] = wgse[j,:]/np.sqrt(wgse[j,0]**2 + wgse[j,1]**2 + wgse[j,2]**2)
        Ymgse = -1*np.cross(aWgse[j,:],zgse)
        aYmgse = Ymgse/(np.sqrt(Ymgse[0]**2 + Ymgse[1]**2 + Ymgse[2]**2))
        Zmgse = np.cross(aWgse[j,:],aYmgse)
        aZmgse = Zmgse/(np.sqrt(Zmgse[0]**2 + Zmgse[1]**2 + Zmgse[2]**2))
        Xmgse = np.cross(aYmgse,aZmgse)
        aXmgse = Xmgse/(np.sqrt(Xmgse[0]**2 + Xmgse[1]**2 + Xmgse[2]**2))
        #Project data along MGSE axes
        datx[j] = sum(aVar[j,:]*aXmgse)
        daty[j] = sum(aVar[j,:]*aYmgse)
        datz[j] = sum(aVar[j,:]*aZmgse)
    VarMgse=np.transpose(np.array((datx,daty,datz)))
    return VarMgse
    
    