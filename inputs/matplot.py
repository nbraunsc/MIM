import matplotlib.pyplot as plt
import numpy as np

x1 = np.load('x1.npy')
y1 = np.load('y1.npy')

x2 = np.load('x2.npy')
y2 = np.load('y2.npy')

x3 = np.load('x3.npy')
y3 = np.load('y3.npy')

x4 = np.load('x4.npy')
y4 = np.load('y4.npy')

fig, (a1, a2, a3, a4) = plt.subplots(nrows=4, sharex=True, sharey=True)# subplot_kw=dict(frameon=False)) # frameon=False removes frames


a1.plot(x1, y1, '-', color = 'grey', label='103 Diamond', linewidth=1)
a1.fill_between(x1, y1, color = 'grey', alpha=0.8)
a1.legend()
a2.plot(x2, y2, '-', color = 'green', label='163 Diamond', linewidth=1)
a2.fill_between(x2, y2, color = 'green', alpha=0.8)
a2.legend()
a3.plot(x3, y3, '-', color = 'purple', label='240 Diamond', linewidth=1)
a3.fill_between(x3, y3, color = 'purple', alpha=0.8)
a3.legend()
a4.plot(x4, y4, '-', color = 'blue', label='593 Diamond', linewidth=1)
a4.fill_between(x4, y4, color = 'blue', alpha=0.8)
a4.legend()
#plt.xlabel('Frequency (cm^-1)')
#plt.ylabel('Intensity')
plt.xlim([1600, 1000])
plt.ylim([0, 50])
plt.subplots_adjust(hspace=.0)
#fig.text(0.5, 0.04, 'common X', ha='center')
#fig.text(0.04, 0.5, 'Intensity', va='center', rotation='vertical')
fig.supylabel('Intensity')
fig.supxlabel('Frequency ($cm^{-1}$)')
plt.show()
