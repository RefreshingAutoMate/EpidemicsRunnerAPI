#This example calibrates the alpha and qq parameters with respect to the dataset defined in subset.lua

import covpred 

def evaluate_model(t0,tend,v):
	index=v.contains_position("alpha")
	
	alphas=[v.get_param(index).get_value_as_double(),4.95954479066937e-7]
	alpha_limits=[63]
	kappa=0.356035567977659
	theta=4.14932000304998e-07
	index=v.contains_position("qq")
	qq=v.get_param(index).get_value_as_double()
	pp=5	
	u0=[753056,2714,0,0,72] #initial conditions
	t_start=0 #simulation start time
	t_end=62 #simulation end time
	
	model=covpred.SEIRD_VARA(alphas, alpha_limits,kappa,theta,qq,pp)
	model.change_minimum_stepsize(0.01)
	model.change_linear_implicit_maximum_stepsize(0.1)
	timepoints, sim_data=model.run_linear_implicit(t_start,u0,t_end) #running the SEIRD model
	
	#We are outputting cumulative data
	cumulated_exposed=model.get_cumulated_exposed()
	cumulated_infected=model.get_cumulated_infected()	
	timepoints_filtered=[]
	datapoints_filtered=[]
	for i in range(len(timepoints)):
		timepoints_filtered.append(covpred.EFloat64(timepoints[i]))
		datapoints_filtered.append(covpred.EFloat64(cumulated_infected[i]))

	return list([list(timepoints_filtered),list(datapoints_filtered)])

options=covpred.PSOOptions()
options.set_max_iterations(10)
options.set_n_groups(5)
options.set_n_particles(30)

enum=covpred.ConfigOutput.File

estimated_parameters=covpred.EVar64Manager()

evaluator=covpred.EpidemicsEvaluation("./","subset_target.lua",evaluate_model)
optimizer=covpred.ParticleSwarmOptimizerEpidemics(options,evaluator)

optimizer.run(estimated_parameters,list(["alpha","qq"]),list([covpred.EFloat64(0),covpred.EFloat64(1e-6),covpred.EFloat64(1),covpred.EFloat64(10)]))
print("Estimated parameters:")
index=estimated_parameters.contains_position("alpha")
print("alpha:" +str(estimated_parameters.get_param(index).get_value_as_double()))
index=estimated_parameters.contains_position("qq")
print("qq:" +str(estimated_parameters.get_param(index).get_value_as_double()))
