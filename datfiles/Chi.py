import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.interpolate import spline

t12,dm21,nsi,chi2=np.loadtxt("chisq.dat",dtype=float,delimiter="\t",usecols=(0,2,4,6),comments="#",unpack=True)
#t12,dm21,nsi,chi2=np.loadtxt("chisq.dat",dtype=float,delimiter="\t",usecols=(0,2,5,6),comments="#",unpack=True)

N=100


nsi_axismin=-2
nsi_axismax=2

histogram1D=np.zeros((N))
ax=np.zeros(N)
ay=np.zeros(N)




for i in range(0,N):
		histogram1D[i]=4000

chi2min=2e10



nsimax = +2.0
nsimin = -2.0

dnsi=(nsimax-nsimin)/N

for i in range(0,N):
	ax[i]=nsimin+dnsi*i



for i in range(0,len(chi2)-1):

	line = int((nsi[i]-nsimin)/dnsi)

	if chi2[i]<histogram1D[line]:
		histogram1D[line]=chi2[i]

	if chi2[i]<chi2min:
		chi2min=chi2[i]
		bfnsi=nsi[i]

ax1=[]
ay1=[]
for i in range(0,N):

		if histogram1D[i]-chi2min < 100.0:
			ax1.append(ax[i])
			ay1.append(histogram1D[i]-chi2min)			



fig = plt.figure()
axes=fig.add_subplot(1,1,1, axisbg="1.0")



#axes.scatter(ax,histogram1D-chi2min,c='red',s=10)
axes.plot(ax1,ay1)

#xnew = np.linspace(nsi_axismin, nsi_axismax, 300)
#power_smooth = spline(ax1, ay1, xnew)
#plt.plot(xnew,power_smooth)


axes.plot([-10,10],[1,1])
axes.plot([-10,10],[2.71,2.71])


print(chi2min)

axes.set_xlim([nsi_axismin,nsi_axismax])
axes.set_ylim([0,3])
axes.set_xlabel('Re$\{[\\tilde{T}]_{e\\mu}\}$')
axes.set_ylabel('$\\Delta\\chi^2$')

plt.title('KamLAND')
plt.legend(('NSI','68% C.L.','90% C.L.'),loc='upper right')
plt.savefig('ReSeu.pdf') 
plt.show()




