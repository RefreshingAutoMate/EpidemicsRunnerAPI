import covpred
import matplotlib.pyplot as plt
import numpy as np

class InitialConditions:
	def __init__(self, u0, alpha,name):
		self.u0=u0
		self.alpha=alpha
		self.name=name
		self.kappa=0.5
		self.theta=0.1
		self.pp=1
		self.qq=1
		self.omega=2e-5
		self.gamma=0.1
		self.epsilon=1e-1
		self.phi=8e-3
		
alpha=[]
alpha_limits=[20]
kappa=[]
theta=[]
qq=[]
pp=[]
sigma=[]
omega=[]
gamma=[]
epsilon=[]
phi=[]
u0=[] #initial conditions of all regions
t_start=0
t_end=17

region_initials=[]

region_initials.append(InitialConditions([760000,0,0,0,0,0,0],[2.08e-6,0],"SK_Frankfurt"))
region_initials.append(InitialConditions([160000,1500,0,0,0,0,0],[0,1e-6],"SK_Offenbach"))

for i in range(0,len(region_initials)):
	alpha+=region_initials[i].alpha
	kappa.append(region_initials[i].kappa)
	theta.append(region_initials[i].theta)
	qq.append(region_initials[i].qq)
	pp.append(region_initials[i].pp)
	phi.append(region_initials[i].phi)
	epsilon.append(region_initials[i].epsilon)	
	gamma.append(region_initials[i].gamma)		
	u0+=region_initials[i].u0
	
	#Nobody from SK_Frankfurt travels to SK_Offenbach
	if i==0:
		omega.append(0)
		for j in range(len(region_initials)):
			sigma.append(0)
			
	#But people from SK_Offenbach are traveling to SK_Frankfurt		
	else:
		omega.append(1e-4)
		for j in range(len(region_initials)):
			sigma.append(1/(len(region_initials)-1))

model=covpred.SEIRD_CK(len(region_initials),alpha, alpha_limits,kappa,theta,qq,pp, sigma, omega, gamma,epsilon, phi)
model.change_minimum_stepsize(0.00001)
model.change_linear_implicit_maximum_stepsize(0.0001)

timepoints, sim_data=model.run_linear_implicit(t_start,u0,t_end)

sim_data=np.array(sim_data).reshape(len(timepoints),int(len(sim_data)/len(timepoints))) #reshape data so that result is not one long contiguous array anymore


fig, (ax1,ax2)=plt.subplots(2)

ax1.plot(timepoints,sim_data[:,0],label="Susceptibles")
ax1.plot(timepoints,sim_data[:,1],label="Exposed")
ax1.plot(timepoints,sim_data[:,2],label="Infected")
ax1.plot(timepoints,sim_data[:,3],label="Recovered")
ax1.plot(timepoints,sim_data[:,4],label="Deceased")
ax1.plot(timepoints,sim_data[:,6],label="Traveling Susceptibles from " +region_initials[0].name )
ax1.plot(timepoints,sim_data[:,7],label="Traveling Exposed from "+region_initials[0].name )

ax2.plot(timepoints,sim_data[:,5],label="Susceptibles")
ax2.plot(timepoints,sim_data[:,6],label="Exposed")
ax2.plot(timepoints,sim_data[:,7],label="Infected")
ax2.plot(timepoints,sim_data[:,8],label="Recovered")
ax2.plot(timepoints,sim_data[:,9],label="Deceased")
ax2.plot(timepoints,sim_data[:,10],label="Traveling Susceptibles from " +region_initials[1].name)
ax2.plot(timepoints,sim_data[:,11],label="Traveling Exposed from "+region_initials[1].name )

ax1.set_title(region_initials[0].name)
ax2.set_title(region_initials[1].name)
ax1.set_xlabel("time")
ax2.set_xlabel("time")
ax1.set_ylabel("People")
ax2.set_ylabel("People")
ax1.legend()
ax2.legend()
plt.show()

