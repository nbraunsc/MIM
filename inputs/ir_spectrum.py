import itertools
import numpy as np
import scipy as sp
from scipy import stats
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys

np.set_printoptions(suppress=True)

folder = sys.argv[1]
#freq = np.array([1968.0385289, 3942.44084817, 4223.63201615])
#y = np.array([3.74769303, 62.00750759, 47.17832453])

path = os.getcwd() + "/" + folder
x = np.load(path + 'freq.npy')
ydata = np.load(path + 'intensity.npy')

freq = np.array(x[68:])
y = np.array(ydata[68:])
print(freq)
print(freq.shape)
print(type(freq))
y_log = np.log(y)
print(y_log)
print(y_log.shape, type(y_log))

def gaussian(freq, y_log):
    xvalues = []
    yvalues = []
    for i in range(0, len(freq)):
        x_data = np.arange(freq[i]-100, freq[i]+100, 10)
        mean = freq[i] #center of peak (freq)
        variance = 100/(np.sqrt((2*np.pi)*y_log[i]))
        #variance = 1/(np.sqrt(2*np.pi)*y_log[i]))
        y_data = stats.norm.pdf(x_data, mean, variance) #generated yvalues
        scale_value = y_log[i]/np.max(y_data)
        scaledy = y_data*scale_value 
        xvalues.append(x_data)
        yvalues.append(scaledy)
        #yvalues.append(y_data)
    return xvalues, yvalues

#def lorentizan(x, x0, gamma):
#    return 1/(np.pi*gamma) * (gamma**2/((x-x0)**2 + gamma**2))

def lorentizan(x, x0, gamma, scale):
    return scale * ((gamma/2)/((x-x0)**2 + (gamma/2)**2))

def mult_lorentizan(freq, y, peak_width):
    xvalues = []
    yvalues = []
    for i in range(0, freq.shape[0]):
        x_data = np.arange(freq[i]-200, freq[i]+200, 10)
        print(len(x_data))
        x0 = freq[i]
        scale = 1/(np.pi*y[i])
        tempy = []
        for j in x_data:
            tempy.append(lorentizan(j, x0, peak_width, scale))
        scale_value = y[i]/np.max(np.array(tempy))
        scaledy = np.array(tempy)*scale_value 
        xvalues.append(list(x_data))
        yvalues.append(list(scaledy))
        #yvalues.append(tempy)
    return xvalues, yvalues


xvalues, yvalues = mult_lorentizan(freq, y, 50)
List_flatx = list(itertools.chain(*xvalues))
List_flaty = list(itertools.chain(*yvalues))

#xvalues, yvalues = gaussian(freq, y_log)

finalx = np.array(List_flatx).flatten()
print(finalx, finalx.shape)
print(type(finalx))
finaly = np.array(List_flaty).flatten()
print(finaly, finaly.shape)
print(type(finaly))

final = np.zeros((finalx.shape[0], 4))
final[:,0] = finalx #freq in cm-1
final[:,1] = finaly #intensities in km/mol   
final[:,2] = -np.log10(finaly/100)    #%abs
final[:,3] = finaly/100 #%transmittance
pd.DataFrame(final).to_csv("xy.csv")


plt.plot(finalx, finaly)
#plt.xlim([np.min(finalx)-5, np.max(finalx)+5])
plt.xlim([5000, 300])
#plt.ylim(max(y_data), min(y_data))
plt.show()


