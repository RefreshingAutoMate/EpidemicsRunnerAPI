#include <iostream>

#include "../models/utility.h"

		template<class T>
		class Problem{
			using F=typename T::value_type;
			
			private:
			F ht=0.1;

			void calc_values(F t, std::array<F,2>& u, std::vector<F>& res){
				std::array<F,2> k1=system(u);
				std::array<F,2> v;
				v[0]=u[0]+ht*0.5*k1[0];
				v[1]=u[1]+ht*0.5*k1[1];
				std::array<F,5> k2= system(v);

				v[0]=u[0]+ht*0.5*k2[0];
				v[1]=u[1]+ht*0.5*k2[1];
				std::array<F,5> k3= system(v);

				v[0]=u[0]+ht*k3[0];
				v[1]=u[1]+ht*k3[1];
				std::array<F,5> k4= system(v);

				u[0]=u[0]+(1.0/6.0)*ht*(k1[0]+2*k2[0]+2*k3[0]+k4[0]);
				u[1]=u[1]+(1.0/6.0)*ht*(k1[1]+2*k2[1]+2*k3[1]+k4[1]);
				res.push_back(u[0]);
				res.push_back(u[1]);	
			}	
			
			
			public:
			const std::array<std::string,2> names={"X1","X2"}; /**< Names of the various SEIRD classes */
			SEIRD();
			
			/*! Sets the step size for the ordinary differential equation solvers.
			@param[in] _ht Step size for the ODE solvers		
			*/	
			void change_minimum_stepsize(F _ht) {
				ht = _ht;
			}
			
			/*! Solves the SEIRD model with an explicit solver of fourth order. Because the solver is explicit, instabilities in the solution
			 * profile can occur. It is recommended to use the linear implicit solver. Results are stored in a vector of vectors. Results are not written to file as of yet.
			@param[in] t0 Simulation start time		
			@param[in] u0 Initial conditions for the five S-E-I-R-D classes
			@param[in] tend end time
			@param[out] std::tuple<std::vector<F>,std::vector<F>> Vector of vectors containing the simulation output			
			*/				
			std::tuple<std::vector<F>,std::vector<F>> run(F t0, const T u0, F tend){
				std::vector<F> res;
				std::vector<F> ts;
				
				res.push_back(u0[0]);
				res.push_back(u0[1]);

				ts.push_back(t0);
				std::array<F,5> u={u0[0],u0[1]};
				
				F t=t0+ht;
				while(t<=tend){
					calc_values(t,u,res);
					ts.push_back(t);
					t+=ht;
				}
				if (t!=tend){
					calc_values(tend,u,res);
					ts.push_back(tend);
				}
				return std::make_tuple(ts,res);			
			}
			
			//!< Returns the system matrix of the ordinary differential equations system determined by the SEIRD model evaluated at time t.
			std::array<F, 2> system(std::array<F, 2>& u, F t=0) {
				std::array<F, 2> res;

				res[0]=1+u[0]*u[0]*u[1]-4*u[0];
				res[1]=3*u[0]-u[0]*u[0]*u[1];
				return res;
			}

			//!< Returns the Jacobi matrix of the ordinary differential equations system determined by the SEIRD model evaluated at time t.
			std::array<F,4> jacobian(const std::array<F,2>& u,F t=0) {
				std::array<F, 4> res;

				res[0]=2*u[1]*u[0]-4;
				res[1]=u[0]*u[0];
				res[2]=3-2*u[1]*u[0];
				res[3]=-u[0]*u[0];

				return res;
			}
			
			/*! Solves the SEIRD model with a linear implicit solver of second order. Results are stored in a vector of vectors. Results are not written to file as of yet.
			@param[in] t0 Simulation start time		
			@param[in] u0 Initial conditions for the five S-E-I-R-D classes
			@param[in] tend end time
			@param[out] std::tuple<std::vector<F>,std::vector<F>> vector of vectors containing the simulation output			
			*/	
			std::tuple<std::vector<F>,std::vector<F>> run_linear_implicit(F t0, const T& u0, F tend) {
				std::array<F, 2> u = { u0[0],u0[1]};
				ug::epi::utility::LinearImplicitSolver23<std::array<F,2>,std::array<F,4>,Problem,F> solver(this,2);
				solver.change_minimum_stepsize(ht);
				solver.change_maximum_stepsize(1);
				solver.change_tol(0.0001);
				auto result=solver.run(t0, u, tend);
				return std::make_tuple(result.first,result.second);
			}
		
		};


int main(){
	
	
	Problem<std::vector<double>> model;
	std::vector<double> u0={1.01,3};
	double t0=0;
	double tend=20;
	double c1=0.8;
	double c2=0;
	model.change_minimum_stepsize(0.001);
	auto result=model.run_linear_implicit(t0,u0,tend);
	std::vector<double> timepoints_linear_implicit = std::get<0>(result);
	std::vector<double> data_linear_implicit = std::get<1>(result);
	ug::epi::write_data("", "linear_implicit_output.txt", timepoints_linear_implicit, data_linear_implicit, model.names, "#");
	
}