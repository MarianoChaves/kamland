import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.cm as cm
import math

t12,dm21,nsi,chi2=np.loadtxt("chisq.dat",dtype=float,delimiter="\t",usecols=(0,2,4,6),comments="#",unpack=True)

N=65




histogram2D=np.zeros((N,N))
ax=np.zeros(N)
ay=np.zeros(N)




for i in range(0,N):
	for j in range(0,N):
		histogram2D[i][j]=4000

chi2min=2e10


t12max = math.log10(10.0)
t12min = math.log10(0.1)

nsimax = +9.0
nsimin = -9.0

dt12=(t12max-t12min)/N
dnsi=(nsimax-nsimin)/N

for i in range(0,N):
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
ay1=[]
az1=[]


for i in range(0,N):
	for j in range(0,N):
		if histogram2D[i][j]<chi2min+20.0:
			ax1.append(ax[i])
			ay1.append(ay[j])
			az1.append(histogram2D[i][j]-chi2min)




#triang = tri.Triangulation(ax1, ay1)
#refiner = tri.UniformTriRefiner(triang)
#tri_refi, z_test_refi = refiner.refine_field(az1, subdiv=3)
#axes.set_aspect('equal')
#axes.tricontourf(tri_refi,z_test_refi,levels=[0.0,2.30,6.18,11.83],cmap='Reds')


#axes.tricontour(ax1,ay1,az1,levels=[2.30,6.18,11.83],linestyles=['solid','dashed','dashdot'])
axes.tricontourf(ax1,ay1,az1,levels=[0.0,2.30,6.18,11.83],cmap='Reds')
#axes.tripcolor(ax1,ay1,az1)

axes.scatter(bft12, bfnsi,c='black',edgecolors='none',s=20)

proxy = [plt.Rectangle((0,0),1,1,fc = pc.get_facecolor()[0]) 
    for pc in axes.collections]

plt.legend(proxy, ['1$\\sigma$','2$\\sigma$','3$\\sigma$'],loc='upper right')

y_axis_NAME='Re$\{[\\tilde{S}]_{e\\tau}\}$'
x_axis_NAME='$\\tan^2\\theta_{12}$'

axes.set_xscale('log')
axes.set_xlim([0.1,10])
axes.set_ylim([-9,9])
axes.set_ylabel(y_axis_NAME)
axes.set_xlabel(x_axis_NAME)

plt.title('KamLAND')
#plt.legend(('1$\\sigma$','2$\\sigma$','3$\\sigma$'),loc='upper right')
plt.savefig('th12_ReSeu.pdf') 
plt.show()




