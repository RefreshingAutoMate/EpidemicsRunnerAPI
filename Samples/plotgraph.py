import matplotlib.pyplot as plt
import numpy as np
import os
file_directory = os.path.dirname(__file__)
print(file_directory)
alldata=np.loadtxt(file_directory+"/linear_implicit_output.txt",delimiter="\t")

t=alldata[:,0]
u1=alldata[:,1]
u2=alldata[:,2]

fig, (ax1,ax2)=plt.subplots(2)
ax1.plot(t,u1,label="u0 approx")
ax1.plot(t,u2,label="u1 approx")

hs=np.zeros(shape=(len(t),1))
hs[0]=0.01
for i in range(len(hs)-1):
	hs[i+1]=t[i+1]-t[i]


ax2.plot(t,hs,label="stepsize")

ax1.set_xlabel("t")
ax2.set_xlabel("h")
ax1.legend()
ax2.legend()
plt.show()