# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 12:11:05 2022

@author: mason

This is a file that contains various functions that I've written for classes and research'
"""
import numpy as np
import math
import scipy
import matplotlib.pyplot as plt
from my_funcs import *
import scipy.stats as stats
from datetime import datetime, timedelta
from scipy.stats import chi2

def my_mean(data):
    ''' This is a function to find the mean of a data set '''
    N=0
    added = 0
    for i in range(len(data)):
        if data[i] != 'nan':
            added += data[i]
            N += 1
    avg = added/N
    return avg
    
def my_var(data):
    #Function to find the standard deviation of a data set
    mean = my_mean(data)
    ss = np.zeros_like(data)
    N = 0
    for i in range(len(data)):
        if data[i] != 'nan':
           ss[i] = (data[i])**2 #find the distance away from the mean
           N += 1
        else:
               ss[i] = float('nan')
    var = my_mean(ss) - mean**2
    
    return var

def my_std_filt(data,n,nor):
    k = 0
    
    while k < nor:
        mean = my_mean(data)
        var = my_var(data)
        std = math.sqrt(var)
        for i in range(len(data)):
            if np.logical_or(data[i] > (n*std+mean), data[i] < (-n*std+mean)) == True:
                data[i] = 'nan'
        k +=1
        return data

def my_gauss(size_in,n):
    # function to create a Guassian distribution
    gauss = np.random.normal(loc=0, scale=1, size=(size_in,n))
    return gauss

def my_laplace(size_in,n):
    # function to create a Laplacian distribution
    x = np.random.rand(size_in,n)-0.5

    laplace = np.sign(x)*np.log(1-2*np.absolute(x))
    return laplace

def my_funky(size_in,n):
    # function to create the funky distribution
	x = 2*np.random.rand(size_in,n) - 1
	funky = np.sign(x)*np.sqrt(np.absolute(x))
	
	return funky

def hist_gauss(data):
    std = np.nanstd(data)   #the std dev ignoring NaN
    mean = np.nanmean(data)  #mean of data ignoring NaN
    
    x = np.arange(mean-6*std,mean+6*std,.01)   #array for the pdf equation to run through
    normal_pdf = 1/np.sqrt(2*3.1415926)/std*np.exp(-(x-mean)**2/(2*std**2)) #equation for normal pdf
    plt.hist(data, density=True, bins=100)  #histogram of data
    plt.plot(x,normal_pdf)
    
    

def  lagcorr(xin,yin,nlag):

    ''' 
    this function takes two vectors of equal length, xin and yin, and
    computes the correlation between xin and yin lagged by nlag,
    where nlag is amount of the lag in indices. 
    '''

    if nlag > 0:
        xlag = xin[0:-nlag]
        ylag =  yin[nlag:]
    
    if nlag < 0:
        xlag = xin[-nlag:]
        ylag = yin[0:nlag]
        
    if nlag == 0:
        xlag = xin
        ylag = yin
            
        
    nan_ind = np.logical_and(np.isfinite(xlag), np.isfinite(ylag))
    xlag = xlag[nan_ind]
    ylag = ylag[nan_ind]
    
    #correlating the lagged vectors
    corr = np.corrcoef(xlag,ylag)
    return corr 

def my_reg(ts):
    N = len(ts)
    x1 = np.ones(N+1)             #constant
    x2 = np.arange(0,N+1)/N     #linear
    M = 2
    X = [x1,x2]
    X = np.array(X)
    D = np.zeros([M,M])
    z = np.zeros(M)
    A = np.zeros(M)
    Y = np.zeros(N)
    for i in np.arange(0,M):
        for j in np.arange(0,M):
            ng = 0
            sxx = 0
            zm = 0
            for k in np.arange(0,N):
                if np.isfinite(ts[k]) == 1:
                    ng += 1
                    sxx += X[i,k]*X[j,k]
                    zm += X[j,k]*ts[k]
            D[i,j] = 1/ng*sxx
            z[j] = 1/ng*zm
    D_inv = np.linalg.inv(D)
    for i in np.arange(0,M):
        for j in np.arange(0,M):
            A[i] += D_inv[i,j]*z[j]
    for i in np.arange(0,N):
        for j in np.arange(0,M):
            Y[i] += A[j]*X[j,i]
    return Y,A

def my_tidalreg(ts,t,tide):
    N = len(ts)
    tide_int = np.interp(t,tide[:,0],tide[:,1])
    N_tide = len(tide_int)
    x1 = np.ones(N_tide)             #constant
    x2 = tide_int   #linear
    M = 2
    X = [x1,x2]
    X = np.array(X)
    #X = np.transpose(np.array(X))
    D = np.zeros([M,M])
    z = np.zeros(M)
    A = np.zeros(M)
    Y = np.zeros(N)
    for i in np.arange(0,M):
        for j in np.arange(0,M):
            ng = 0
            sxx = 0
            zm = 0
            for k in np.arange(0,N):
                if np.isfinite(ts[k]) == 1:
                    ng += 1
                    sxx += X[i,k]*X[j,k]
                    zm += X[j,k]*ts[k]
            D[i,j] = 1/ng*sxx
            z[j] = 1/ng*zm
    D_inv = np.linalg.inv(D)
    for i in np.arange(0,M):
        for j in np.arange(0,M):
            A[i] += D_inv[i,j]*z[j]
    for i in np.arange(0,N):
        for j in np.arange(0,M):
            Y[i] += A[j]*X[j,i]
    return Y,A

def my_tidal_harmonic(ts,t):
    N = len(ts)
    tides_cos = np.cos(2*3.1415/(12.4*60*60)*t)
    tides_sin = np.sin(2*3.1415/(12.4*60*60)*t)
    N_tide = len(tides_sin)
    x1 = np.ones(N_tide)             #constant
    x2 = tides_cos   #linear
    x3 = tides_sin
    M = 3
    X = [x1,x2,x3]
    X = np.array(X)
    #X = np.transpose(np.array(X))
    D = np.zeros([M,M])
    z = np.zeros(M)
    A = np.zeros(M)
    Y = np.zeros(N)
    for i in np.arange(0,M):
        for j in np.arange(0,M):
            ng = 0
            sxx = 0
            zm = 0
            for k in np.arange(0,N):
                if np.isfinite(ts[k]) == 1:
                    ng += 1
                    sxx += X[i,k]*X[j,k]
                    zm += X[j,k]*ts[k]
            D[i,j] = 1/ng*sxx
            z[j] = 1/ng*zm
    D_inv = np.linalg.inv(D)
    for i in np.arange(0,M):
        for j in np.arange(0,M):
            A[i] += D_inv[i,j]*z[j]
    for i in np.arange(0,N):
        for j in np.arange(0,M):
            Y[i] += A[j]*X[j,i]
    return Y,A

def datetime2matlabdn(dt):
   ord = dt.toordinal()
   mdn = dt + timedelta(days = 366)
   frac = (dt-datetime(dt.year,dt.month,dt.day,0,0,0)).seconds / (24.0 * 60.0 * 60.0)
   return mdn.toordinal() + frac

def datenum2datetime(datenum):
    """
    Convert Matlab datenum into Python datetime.

    :param datenum: Date in datenum format
    :return:        Datetime object corresponding to datenum.
    """
    days = datenum % 1
    hours = days % 1 * 24
    minutes = hours % 1 * 60
    seconds = minutes % 1 * 60
    milliseconds = seconds % 1 * 1000
    return datetime.fromordinal(int(datenum)) \
           + timedelta(days=int(days)) \
           + timedelta(hours=int(hours)) \
           + timedelta(minutes=int(minutes)) \
           + timedelta(seconds=int(seconds)) \
           + timedelta(milliseconds=round(milliseconds)) \
           - timedelta(days=366)
def my_long_lag(x,y,minlag,maxlag,dt):
    K = int(abs((maxlag-minlag)/dt));
    corr = np.zeros(2*K)
    xlag = np.zeros(2*K)
    N = len(x)
    M = maxlag - minlag
    
    for k in np.arange(K):
        ng = 0
        sx = 0
        sxx = 0
        sy = 0
        syy = 0
        sxy = 0
        for i in np.arange(N-(k+minlag)):
            x1 = x[i]
            y1 = y[i+(k+minlag)]
            if np.logical_and(np.isfinite(x1)==1,np.isfinite(y1)==1):
                ng += 1
                sx += x1
                sy += y1
                sxx += x1**2
                syy += y1**2
                sxy += x1*y1
        corr[k+K] = sxy/ng - (sx/ng*sy/ng)/(((sxx/ng - (sx/ng)**2)**(1/2))*((syy/ng - (sy/ng)**2)**(1/2)))
        xlag[k+K] = dt*k
        corr_high = corr
    for k in np.arange(K):
        ng = 0
        sx = 0
        sxx = 0
        sy = 0
        syy = 0
        sxy = 0
        for i in np.arange(N-(k+minlag)):
            x1 = x[i]
            y1 = y[i+(k+minlag)]
            if np.logical_and(np.isfinite(x1)==1,np.isfinite(y1)==1):
                ng += 1
                sx += x1
                sy += y1
                sxx += x1**2
                syy += y1**2
                sxy += x1*y1
        corr[k+K] = sxy/ng - (sx/ng*sy/ng)/(((sxx/ng - (sx/ng)**2)**(1/2))*((syy/ng - (sy/ng)**2)**(1/2)))
        xlag[-k+K] = -dt*k
        corr_low = corr
        
    SA = (1/2/M)*((sum(corr_low**2))+sum(corr_high)**2)
    N_star = 1/SA
    return N_star,corr_low,corr_high
        

def f_spectra(ts, dt, nbands, nens, alpha):
    N = len(ts)
    DOF = nbands*nens*2
    # Ensemble average
    K = N // nens
    jj = np.arange(0, K // 2 + 1)
    fk = jj / (K * dt)
    df = (fk > 0) & (fk < 1 / (2 * dt))
    fk = fk[df]

    So = np.zeros(2*(1+len(fk)))
    Sj = np.zeros((2*(1+len(fk)), nens))

    for k in range(nens):
        tsk = ts[k * K:(k + 1) * K]
        tsk = tsk - np.mean(tsk)
        Gj = (1 / K) * np.fft.fft(tsk)
        G_pro = np.abs(Gj) ** 2
        Sj[:, k] = 2 * K * dt * G_pro
        So = So + Sj[:, k]

    Sk = So / nens

    # Band average
    B = len(Sk) // nbands
    S = np.zeros(B)
    f = np.zeros(B)

    for b in range(B):
        S[b] = np.mean(Sk[b * nbands:(b + 1) * nbands])
        f[b] = np.mean(fk[b * nbands:(b + 1) * nbands])

    crit = chi2.ppf(alpha,DOF)
    CI = [DOF / chi2.ppf(1-alpha/2,DOF), DOF / chi2.ppf(alpha/2,DOF)]
    return S, f, DOF, CI

def my_chi2(a,df):
    
    lower = scipy.stats.chi2.ppf(a,df)
    upper = scipy.stats.chi2.ppf(1-a/2,df)
    conf = np.array([lower,upper])
    
    return conf
    
    
def my_cross_spectrum(ts1,ts2,dt,Nb,Ne):
        
        N = len(ts1)
        
        K = np.floor(N/Ne)
        j = np.arange(K/2)
        k = np.arange(K)
        fe = j/K/dt
        S1 = np.zeros_like(k,'complex')
        S2 = np.zeros_like(k,'complex')

        Co = np.zeros_like(k,'complex')
        Qo = np.zeros_like(k,'complex')
        Se1 = np.zeros_like(k,'complex')
        Se2 = np.zeros_like(k,'complex')
        for i in np.arange(Ne):
            ens_ts1 = ts1[int(i*K):int((i+1)*K)]
            ens_ts1 = ens_ts1 - np.mean(ens_ts1)
            G1 = np.fft.fft(ens_ts1)/K
            A1 = np.real(G1)
            B1 = np.imag(G1)
            Se1= K*dt*(A1**2 + B1**2)
            S1 = S1 + Se1
            
            
            ens_ts2 = ts2[int(i*K):int((i+1)*K)]
            ens_ts2 = ens_ts2 - np.mean(ens_ts2)
            G2 = np.fft.fft(ens_ts2)/K
            A2 = np.real(G2)
            B2 = np.imag(G2)
            Se2 = K*dt*(A2**2 + B2**2)
            S2 = S2 + Se2
                
            C = K*dt*(A1*A2 + B1*B2)
            Co = Co + C
            Q = K*dt*(B1*A2 - A1*B2)
            Qo = Qo + Q
        
        S_ens1 = S1/Ne
        S_ens2 = S2/Ne
        Co_ens = Co/Ne
        Qo_ens = Qo/Ne
        Cf = np.zeros_like(k,'complex')
        Qf = np.zeros_like(k,'complex')
        Sb1 = np.zeros_like(k,'complex')
        Sb2 = np.zeros_like(k,'complex')
        f = np.zeros_like(k,'complex')
        B = int(np.floor(len(S_ens1)/Nb))
        
        for i in np.arange(B):
            Sb1[int(i)] = np.mean(S_ens1[int((i-1)*Nb):int(i*B)])
            Sb2[int(i)] = np.mean(S_ens2[int((i-1)*Nb):int(i*B)])
            Cf[int(i)] = np.mean(Co_ens[int((i-1)*Nb):int(i*B)])
            Qf[i] = np.mean(Qo_ens[int((i-1)*Nb):int(i*B)])
            f[i] = np.mean(fe[int((i-1)*Nb):int(i*B)])
        
        S12 = Cf - 1j*Qf
        Coh2 = (Cf**2 + Qf**2)/S1/S2
        phi = np.arctan(-Qf/C)*180/3.1415
        
        return S12, Coh2, phi, f