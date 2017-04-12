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
    
    
def delta_outlier(series, window=50, advance = 50, Z=3,back=1):
    """
    Identify locations of sudden change accounting for Gaussian Noise 
    Moderate Sensitivity
    Very Robust to Noise
    """
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

# Calculate Rate of Sudden Change
def delta_rate(series)
    
# Detect Sudden Change is Signal Noise
def delta_noise(series):
    
    
# Test Sequence (Clean Signal)
testA = np.zeros(300)
testA = np.append(testA,np.ones(300)*400)
testA = np.append(testA,np.zeros(300))
testA = np.append(testA,np.ones(300)*400)
testA = np.append(testA,np.zeros(300))
testA = np.append(testA,np.ones(300)*400)
testA = np.append(testA,np.zeros(300))
testA = np.append(testA,np.ones(300)*400)
testA = np.append(testA,np.zeros(300))
testA = np.append(testA,np.ones(300)*400)

#Noise Signal
noise = np.random.rand(len(testA))*100

#Noisy Signal
noisy_testA =testA+noise

plt.plot(testA)
print(delta_slope(testA,tail_dist=2,lead_dist=2, advance = 10,threshold=1.2 ,back=1))
#print(delta_outlier(test,advance=30,Z=3))
