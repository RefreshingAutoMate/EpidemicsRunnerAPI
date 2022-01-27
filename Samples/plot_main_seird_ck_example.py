import matplotlib.pyplot as plt
import numpy as np
import os
file_directory = os.path.dirname(__file__)
print(file_directory)
alldata=np.loadtxt(file_directory+"/aseird_ck_output.txt",delimiter=",")

t=alldata[:,0]
u1=alldata[:,1]
u2=alldata[:,2]
print(u1[15])
fig, (ax1,ax2)=plt.subplots(2)
ax1.plot(t,alldata[:,1],label="Area1 Susceptibles")
ax1.plot(t,alldata[:,2],label="Area1 Exposed")
ax1.plot(t,alldata[:,3],label="Area1 Infected")
ax1.plot(t,alldata[:,4],label="Area1 Recovered")
ax1.plot(t,alldata[:,5],label="Area1 Deceased")
#ax1.plot(t,alldata[:,6],label="Area1 C")
#ax1.plot(t,alldata[:,7],label="Area1 K")
#ax1.plot(t,u2,label="Area1 Exposed")

ax2.plot(t,alldata[:,8],label="Area2 Susceptibles")
ax2.plot(t,alldata[:,9],label="Area2 Exposed")
ax2.plot(t,alldata[:,10],label="Area2 Infected")
ax2.plot(t,alldata[:,11],label="Area2 Recovered")
ax2.plot(t,alldata[:,12],label="Area2 Deceased")
ax2.plot(t,alldata[:,13],label="Area2% C")
ax2.plot(t,alldata[:,14],label="Traveling Exposed from Area2")

ax1.set_xlabel("time")
ax2.set_xlabel("time")
ax1.legend()
ax2.legend()
plt.show()