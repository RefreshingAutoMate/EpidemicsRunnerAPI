#include "../models/seird_ck.h"
#include "../models/writer.h"
#include<vector>
#include <iostream>
#include<string>

/*
This function evaluates the SEIRD model and writes the output to file.
*/

class InitialConditions{
	
	public:
	InitialConditions(std::vector<double> _u0, std::vector<double> _alpha, std::string _name):u0(_u0),alpha(_alpha),name(_name){
	}
	
	std::vector<double> u0;
	double kappa=0.5;
	double theta=0.1;
	double pp=1;
	double qq=1;
	double omega=2e-5;
	double gamma=0.1;
	double epsilon=1e-1;
	double phi=8e-3;
	std::vector<double> alpha;
	
	std::string name;
	
};

int main(int argc, char *argv[]){

	std::vector<double> alpha;
	std::vector<double> alpha_limits={20};
	
	std::vector<double> kappa;
	std::vector<double> theta;
	std::vector<double> qq;
	std::vector<double> pp;
	std::vector<double> sigma;
	std::vector<double> omega;
	std::vector<double> gamma;
	std::vector<double> epsilon;
	std::vector<double> phi;
	
	std::vector<InitialConditions> region_initials;
	
	region_initials.push_back(InitialConditions({760000,0,0,0,0,0,0},{2.08e-6,0},"SK_Frankfurt"));
	region_initials.push_back(InitialConditions({160000,1500,0,0,0,0,0},{0, 1e-6},"SK_Offenbach"));
	
	constexpr int regions=2;

	std::vector<double> u0;
	for (int i=0;i<regions;i++){
		for (int j=0;j<region_initials[i].alpha.size();j++){
			alpha.push_back(region_initials[i].alpha[j]);
		}
		if (i==0){
			omega.push_back(0);
			for (int j=0;j<regions;j++){
				sigma.push_back(0.0); 
			}
		}
		else{
			omega.push_back(1e-4);
			for (int j=0;j<regions;j++){
				sigma.push_back(1); 
			}
		}
		kappa.push_back(region_initials[i].kappa);
		theta.push_back(region_initials[i].theta);
		qq.push_back(region_initials[i].qq);
		pp.push_back(region_initials[i].pp);
		//omega.push_back(region_initials[i].omega);
		phi.push_back(region_initials[i].phi);
		gamma.push_back(region_initials[i].gamma);
		epsilon.push_back(region_initials[i].epsilon);
		u0.push_back(region_initials[i].u0[0]);
		u0.push_back(region_initials[i].u0[1]);		
		u0.push_back(region_initials[i].u0[2]);	
		u0.push_back(region_initials[i].u0[3]);	
		u0.push_back(region_initials[i].u0[4]);	
		u0.push_back(region_initials[i].u0[5]);	
		u0.push_back(region_initials[i].u0[6]);	
	}
	
	double t_start=0;
	double t_end=14;
	
	ug::epi::SEIRD_CK<std::vector<double>,regions> seird_model(alpha, alpha_limits,kappa,theta,qq,pp, sigma, omega, gamma,epsilon, phi);
	
	
	seird_model.change_minimum_stepsize(0.00001);
	seird_model.change_linear_implicit_maximum_stepsize(0.0001);
	seird_model.change_linear_implicit_tol(0.0000001);

	auto result=seird_model.run_linear_implicit(t_start,u0,t_end);
	auto timepoints=std::get<0>(result);
	auto data=std::get<1>(result);	
	std::cout<<"Done\n";
	
	std::vector<std::string> names;
	for (int i=0;i<region_initials.size();i++){
		names.push_back(region_initials[i].name+"_Susceptibles");
		names.push_back(region_initials[i].name+"_Exposed");
		names.push_back(region_initials[i].name+"_Infected");
		names.push_back(region_initials[i].name+"_Recovered");
		names.push_back(region_initials[i].name+"_Deceased");
		names.push_back(region_initials[i].name+"_C");
		names.push_back(region_initials[i].name+"_K");
	}
	
	
	ug::epi::write_data(argv[0], "seird_ck_output.txt", timepoints, data,names,"#","",",");
	
	std::cout<<"Ende\n";
}
