import sys
import os
import matplotlib.pyplot as plt
import numpy as np

one = sys.argv[1]
two = sys.argv[2]
three = sys.argv[3]
four = sys.argv[4]
folders = [one, two, three, four]
gridpoints = np.arange(900, 3200, 1)

def get_xy(folder):
    path = os.getcwd() + "/" + folder
    vib = np.load(path + 'freq.npy')
    freqs = vib[6:] 
    intense = np.load(path + 'intensity.npy')
    return freqs, intense

def ys(alpha, freqs, intense):
    y = np.zeros((len(gridpoints)))
    for x in range(0, len(gridpoints)):
        for f in range(0, len(freqs)):
            value = -alpha*(gridpoints[x]-freqs[f])**2
            y[x] += intense[f]*np.exp(value)
    return y

#print(gridpoints)
#print(y)

#y1 = ys(10)
#y2 = ys(0.02)

y_values = []
for i in range(0, len(folders)):
    freq, intensity = get_xy(folders[i])
    yi = ys(0.02, freq, intensity)
    y_values.append(yi)
    print("Done with folder:", folders[i])
#plt.plot(gridpoints, y1, '-', color='grey')
#plt.plot(gridpoints, y2, '-', color='blue')
#plt.xlim(900, 1600)
#plt.ylim(0, 50)
#plt.show()
#exit()

fig, (a1, a2, a3, a4) = plt.subplots(nrows=4, sharex=True, sharey=True)# subplot_kw=dict(frameon=False)) # frameon=False removes frames
a1.plot(gridpoints, y_values[0], '-', color = 'grey', label='103 Diamond', linewidth=1)
#a1.fill_between(x1, y1, color = 'grey', alpha=0.8)
a1.legend()
a2.plot(gridpoints, y_values[1], '-', color = 'green', label='163 Diamond', linewidth=1)
#a2.fill_between(x2, y2, color = 'green', alpha=0.8)
a2.legend()
a3.plot(gridpoints, y_values[2], '-', color = 'purple', label='240 Diamond', linewidth=1)
#a3.fill_between(x3, y3, color = 'purple', alpha=0.8)
a3.legend()
a4.plot(gridpoints, y_values[3], '-', color = 'blue', label='593 Diamond', linewidth=1)
#a4.fill_between(x4, y4, color = 'blue', alpha=0.8)
a4.legend()

plt.subplots_adjust(hspace=.0)
fig.supylabel('Intensity')
fig.supxlabel('Frequency ($cm^{-1}$)')

plt.xlim(900, 1600)
plt.ylim(0, 60)
plt.show()
