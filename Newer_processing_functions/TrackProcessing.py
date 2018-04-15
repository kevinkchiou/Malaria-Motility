# -*- coding: utf-8 -*-
"""
Spyder Editor

created by Kevin K Chiou 4/15/2018

Automated processing For use with C. Hopp's data files from automated tracking
"""

import pandas as pd
import numpy as np

def cleanData(df):
    nt = (len(df.columns)-1) // 3
    df_out = pd.DataFrame()
    df_out['time'] = df[df.columns[0]]/1000.0
    for i in range(nt):
        idx=[1+3*i,2+3*i,3+3*i]
        vals=df[df.columns[idx]]
        metric = [vals.diff().abs().sum().values,idx]
        exclude = metric[1][np.argmin(metric[0])]
        idx.remove(exclude)
        df_out[''.join([str(i),'x'])] = df[df.columns[idx[0]]]
        df_out[''.join([str(i),'y'])] = df[df.columns[idx[1]]]
    return df_out

def timeWinLinInterp(df,dt=3):
    df_out=pd.DataFrame()
    t=df.time.values
    pos=df.drop('time',axis=1).values
    tf = np.arange(np.min(t),np.max(t),dt)
    posf = np.zeros((len(tf),pos.shape[1]))
    
    for i,v in enumerate(tf):
        if i==0:
            posf[0,:] = pos[0,:]
        else:
            idx=len(t[t<=v])-1
            dt1,dt2 = v - t[idx], t[idx+1]-v
            posf[i,:] = pos[idx,:]*dt2/(dt1+dt2) + pos[idx+1,:]*dt1/(dt1+dt2)
    
    df_out['time']=tf
    for i,v in enumerate(df.columns):
        if v != 'time':
            df_out[v] = posf[:,i-1]
    df_vel = df_out.drop('time',axis=1).diff()/dt
    df_vt = pd.DataFrame()
    for i in range(df_vel.shape[1]//2):
        df_vt[str(i)] = np.sqrt(df_vel.values[:,2*i]**2 + df_vel.values[:,2*i+1]**2)
    
    return df_out,df_vel,df_vt

if __name__ == "__main__":
    import os
    
    for (_,_,filearrays) in os.walk('./'):
        for file in filearrays:
            if ((file.find('.xls')>-1) & (file.find('_processed.xls')==-1)):
                df0 = pd.read_excel(file,header=1)
                df1 = cleanData(df0)
                dfp,dfv,dfvt = timeWinLinInterp(df1)

                sname = file.split('.')
                sname[-1]='_processed.xls'
                writer = pd.ExcelWriter(''.join(sname))
                dfp.to_excel(writer,'Positions')
                dfv.to_excel(writer,'Velocities')
                dfvt.to_excel(writer,'Speed')
                writer.save()
