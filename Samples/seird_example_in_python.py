import covpred
import matplotlib.pyplot as plt
import numpy as np


alpha=4.95954479066937e-07
kappa=0.356035567977659
theta=4.14932000304998e-07
qq=8
pp=5


u0=[753056,2714,0,0,72]
t_start=0
t_end=43

alphas=[4.95954479066937e-07,4.95954479066937e-8]
alpha_limits=[30]

model=covpred.SEIRD_VARA(alphas, alpha_limits,kappa,theta,qq,pp)

'''
model=covpred.SEIRD(alpha,kappa,theta,qq,pp)

timepoints, sim_data=model.run_linear_implicit(t_start,u0,t_end)
'''
model.change_minimum_stepsize(0.01)
model.change_linear_implicit_maximum_stepsize(0.1)
timepoints, sim_data=model.run_linear_implicit(t_start,u0,t_end)

sim_data=np.array(sim_data).reshape(len(timepoints),int(len(sim_data)/len(timepoints))) #reshape data so that result is not one long contiguous array anymore

plt.plot(timepoints,sim_data[:,0],label="Susceptibles")
plt.plot(timepoints,sim_data[:,1],label="Exposed")
plt.plot(timepoints,sim_data[:,2],label="Infected")
plt.plot(timepoints,sim_data[:,3],label="Recovered")
plt.plot(timepoints,sim_data[:,4],label="Deceased")
plt.legend()
plt.xlabel('time')
plt.ylabel('People')
plt.show()
cumulated_exposed=model.get_cumulated_exposed()
cumulated_infected=model.get_cumulated_infected()
plt.plot(timepoints,cumulated_exposed)

plt.show()
