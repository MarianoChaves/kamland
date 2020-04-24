import numpy as np
import matplotlib.pyplot as plt
import math

t12,dm21,nsi,chi2=np.loadtxt("chisq.dat",dtype=float,delimiter="\t",usecols=(0,2,4,6),comments="#",unpack=True)

N=1000



histogram2D=np.zeros((N,N))
ax=np.zeros(N)
ay=np.zeros(N)




for i in range(0,N-1):
	for j in range(0,N-1):
		histogram2D[i][j]=300.0

chi2min=2e10


t12max = math.log10(10.0)
t12min = math.log10(0.1)

nsimax = +10.0
nsimin = -10.0

dt12=(t12max-t12min)/N
dnsi=(nsimax-nsimin)/N

for i in range(0,N-1):
	ax[i]=math.pow(10,t12min+dt12*i)
	ay[i]=nsimin+dnsi*i



for i in range(0,len(chi2)-1):

	line = int((t12[i]-t12min)/dt12)
	collumm = int((nsi[i]-nsimin)/dnsi)

	if chi2[i]<histogram2D[line][collumm]:
		histogram2D[line][collumm]=chi2[i]

	if chi2[i]<chi2min:
		chi2min=chi2[i]
		bft12=pow(10,t12[i])
		bfnsi=nsi[i]



fig = plt.figure()
axes=fig.add_subplot(1,1,1, axisbg="1.0")



ax1=[]
ax2=[]
ax3=[]
ay1=[]
ay2=[]
ay3=[]

for i in range(0,N-1):
	for j in range(0,N-1):
		if histogram2D[i][j]<chi2min+11.83:
			ax1.append(ax[i])
			ay1.append(ay[j])
		if histogram2D[i][j]<chi2min+6.18:
			ax2.append(ax[i])
			ay2.append(ay[j])
		if histogram2D[i][j]<chi2min+2.30:
			ax3.append(ax[i])
			ay3.append(ay[j])




axes.scatter(ax1, ay1,c='red',edgecolors='none',s=20)
axes.scatter(ax2, ay2,c='blue',edgecolors='none',s=20)
axes.scatter(ax3, ay3,c='green',edgecolors='none',s=20)



axes.scatter(bft12,bfnsi, alpha=1,c="black",edgecolors='none',s=30)

axes.set_xscale('log')
axes.set_xlim([0.1,10])
axes.set_ylim([-2,2])
axes.set_ylabel('Im$\{[\\tilde{S}]_{e\\mu}\}$')
axes.set_xlabel('$\\tan^2\\theta_{12}$')

plt.title('KamLAND')
plt.legend(('1$\\sigma$','2$\\sigma$','3$\\sigma$'),loc='upper right')
plt.show()




