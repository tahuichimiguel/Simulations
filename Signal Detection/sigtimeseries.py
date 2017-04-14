import numpy as np
import matplotlib.pyplot as plt

# Detect Locations of Sudden Changes in Signal
def delta_slope(series,tail_dist=10,lead_dist=10, advance = 50,threshold=1.1 ,back=0):
    """
    Identify locations of sudden change using a 3-Node Stencil
    Very High Sensitivity 
    Does Very Poorly in Presence of Noise
    """
    out = []   
    temp=0

    tail=0
    center = tail_dist
    lead_dist = tail_dist+lead_dist
    while(center<(len(series)-lead_dist)):
        tail = center-tail_dist
        lead = center+lead_dist
        
        temp = np.abs(series[lead]-series[center])/np.abs(series[center]-series[tail])
        if( temp>threshold ):
            out.append(center-back)
            center+=advance
            continue
        center+=1  
    return out
    
    
def delta_outlier(series, window=50, advance = 100, Z=6,back=1):
    """
    Identify locations of sudden change accounting for Gaussian Noise Distribution 
    Moderate Sensitivity
    Robust to Noise
    """
    
    if advance<window:
        print('Advance Must Be Greater Than or Equal to Window')
        return
    
    out = []
    
    S = np.std(series[0:window])
    X_bar= np.mean(series[0:window])
    i=window
    while(i<(len(series))):
        if(np.abs(series[i]-X_bar)>Z*S):
            out.append(i-back)
            i+=advance
            S = np.std(series[i-window:i])
            X_bar= np.mean(series[i-window:i])
            continue
        i+=1        
    return out

# Calculate Time to Reach Specified Fraction of Final Value
def delta_rate(series,start_indxs,end_indxs,base_dist=50, span_dist=50,threshold=0.9):
    """
    Counts the number of samples it takes for the change in signal to reach
    the fraction of the maximum change specified by threshold
    """
    
    if len(start_indxs)!=len(end_indxs):
        print('Start Indices and End Indices Must Have Same Length')
        return
        
    out=[]
    baseline=0
    span=0
    delta100=0
    for i in range(len(start_indxs)):
        baseline = np.mean(series[(start_indxs[i]-base_dist):start_indxs[i]])
        span = np.mean(series[(end_indxs[i]-span_dist):end_indxs[i]])
        delta100=np.abs(span-baseline)
        
        count=0
        indx=start_indxs[i]        
        while(np.abs(series[indx]-baseline)<threshold*delta100):
            count+=1
            indx+=1    
        out.append(count)
        
    return out
    
# Detect Sudden Change is Signal Noise (Assumming Gaussian Noise)
def delta_noise_norm(series,window=30,advance=30, F_crit=3):
    out = []
        
    indx=1
    ref_var=np.var(series[0:window])
    last_test=0
    test_var=np.var(series[window:window+window])
    while(indx<len(series)-window):
        if test_var>ref_var:
            F=test_var/ref_var
        else:
            F=ref_var/test_var
    
        if F>F_crit:
            
            if ref_var!=last_test:
                out.append((indx,test_var,ref_var))
                
            last_test=test_var    
            ref_var=np.var(series[indx:indx+window])
            test_var=np.var(series[indx+window:indx+window+window])
            indx=indx+advance
            continue
        
        indx+=1
        test_var=np.var(series[indx:indx+window])

    return out
    
def spans(series,start_indxs,end_indxs,base_dist=50, span_dist=50):
    out=[]
    
    for i in range(len(start_indxs)):
         baseline=np.mean(series[start_indxs[i]-base_dist:start_indxs[i]])
         span=np.mean(series[end_indxs[i]-span_dist:end_indxs[i]])
         out.append(span-baseline)
     
    return out
    
    
