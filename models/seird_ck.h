#pragma once
#include <vector>
#include <array>
#include <tuple>
#include <string>
#include "utility.h"

namespace ug{
	namespace epi{
	
		template<class T>
		class SEIRD_CK{
			using F=typename T::value_type;
			
			private:
			std::vector<F> alpha;
			std::vector<F> alpha_limits;
			int N;
			int system_dim;
			std::vector<F> kappa;
			std::vector<F> theta;
			std::vector<F> qq;
			std::vector<F> pp;
			std::vector<F> sigma;
			std::vector<F> omega;
			std::vector<F> gamma;
			std::vector<F> epsilon;
			std::vector<F> phi;
			
			F ht=0.01;
			F ht_max=0.125;
			std::vector<F> cumulated_infected_of_last_run=std::vector<F>(N,F(-1));
			std::vector<F> cumulated_exposed_of_last_run=std::vector<F>(N,F(-1));
			std::vector<std::vector<F>> cumulated_infected_of_last_run_container;
			std::vector<std::vector<F>> cumulated_exposed_of_last_run_container;	
			
			double tol=0.01;
			
			F eval_alpha(F t,int index){
				bool interpolation=false;
				int alphaindex=0;
				int limitindex=-1;
				
				F eps=F(0.5);
				F offset_alpha=index*(alpha.size()/N);
				
				for (int i=0;i<alpha_limits.size();i++){
					if (t>(alpha_limits[i]-eps)){
						alphaindex=i+1;
						limitindex=i;
						if (t>=alpha_limits[i]){
							interpolation=false;
						}
						else{
							interpolation=true;
						}
					}
					
				}
	
				if (interpolation==false){
					return alpha[offset_alpha+alphaindex];
				}else{
					
					F x0= alpha_limits[limitindex]-eps;
					F x1=alpha_limits[limitindex];
					F y0=alpha[offset_alpha+alphaindex-1];
					F y1=alpha[offset_alpha+alphaindex];
					F R;
					F C;
					F c;
					F d=std::sqrt((double)((y0-y1)*(y0-y1)));
					if (y1>y0){
						R=x1;
						C=x0;
						c=y0;
					}
					else{
						R=x0;
						C=x1;
						c=y1;
					}
					return y0+d*std::exp(double(-((t-R)*(t-R))/((t-C)*(t-C))));
				}
				
			}
			
			public:
			const std::array<std::string,7> names={"Susceptibles","Exposed", "Infected", "Recovered", "Deaths", "Traveling Susceptibles", "Traveling Exposed"};
			SEIRD_CK(int _N, std::vector<F> _alpha,std::vector<F>& _alpha_limits ,std::vector<F>& _kappa, std::vector<F>& _theta,std::vector<F>& _qq, std::vector<F>& _pp, std::vector<F>& _sigma, std::vector<F>& _omega, std::vector<F> _gamma,std::vector<F>& _epsilon, std::vector<F>& _phi):N(_N),system_dim(_N*7),alpha(_alpha),alpha_limits(_alpha_limits),kappa(_kappa),theta(_theta),qq(_qq),pp(_pp), sigma(_sigma), omega(_omega),phi(_phi), epsilon(_epsilon), gamma(_gamma){
				
			}	
			

			void change_minimum_stepsize(F _ht) {
				ht = _ht;
			}
			
			void change_linear_implicit_tol(double _tol){
				tol=_tol;
			}
			void change_linear_implicit_maximum_stepsize(double _ht_max){
				ht_max=_ht_max;
			}			
			
			void update_metainfo(std::vector<F>& u, F t, F h){	
			
				for (int i=0;i<N;i++){
					F _alpha=eval_alpha(t,i);
					F G = u[i*system_dim+0]; // Gesunde (Susceptibles)
					F A = u[i*system_dim+1]; // Angesteckte (Exposed)
					F K = u[i*system_dim+2]; // Kranke (Infected)							
					cumulated_exposed_of_last_run[i]+=_alpha * G * A*h;
					cumulated_infected_of_last_run[i]+=(kappa[i] / qq[i]) * A*h;
					cumulated_exposed_of_last_run_container[i].push_back(cumulated_exposed_of_last_run[i]);
					cumulated_infected_of_last_run_container[i].push_back(cumulated_infected_of_last_run[i]);								
				}
		
			}
			
			std::vector<F> system(std::vector<F>& u, F t) {
				
				std::vector<F> res(7*N);
				for (int i=0;i<N;i++){
					
					F S=u[i*7];
					F E=u[i*7+1];
					F I=u[i*7+2];
					F R=u[i*7+3];
					F D=u[i*7+4];
					F C=u[i*7+5];
					F K=u[i*7+6];
					F sum=0;
					for (int j=0;j<N;j++){
						if (j!=i){
							sum+=sigma[j*N+i]*u[j*7+5]*E;
						}
					}
					for (int j=0;j<N;j++){
						if (j!=i){
							sum+=sigma[j*N+i]*S*u[j*7+6];
						}
					}
					res[i*7]=-eval_alpha(t,i)*(S*E+sum)-omega[i]*S+gamma[i]*C; // S
					sum=0;
					for (int j=0;j<N;j++){
						if (j!=i){
							sum+=sigma[j*N+i]*u[j*7+5]*E;
						}
					}					
					for (int j=0;j<N;j++){
						if (j!=i){
							sum+=sigma[j*N+i]*S*u[j*7+6];
						}
					}
					
					
					res[i*7+1]=eval_alpha(t,i)*(S*E+sum)-phi[i]*E+epsilon[i]*K-(1/qq[i])*E; //E
					
					res[i*7+2]=(kappa[i]/qq[i])*E-(1/pp[i])*I; //I
					
					res[i*7+3]=((1-kappa[i])/qq[i])*E+((1-theta[i])/qq[i])*I;    //R

					res[i*7+4]=(theta[i]/qq[i])*I; //D
					
					sum=0;
					for (int j=0;j<N;j++){
						if (j!=i){
							sum+=sigma[i*N+j]*C*u[j*7+1];
						}
					}
					
					res[i*7+5]=omega[i]*S-gamma[i]*C-sum; //C
					
					sum=0;
					for (int j=0;j<N;j++){
						if (j!=i){
							sum+=alpha[j]*sigma[i*7+j]*C*u[j*7+1];
						}
					}
					
					res[i*7+6]=phi[i]*E-epsilon[i]*K+sum; //K
				}
				
				return res; 
			}
		
			std::vector<F> jacobian(const std::vector<F>& u, F t) {
				std::vector<F> res(N*7*N*7);
			
				for (int i = 0; i < res.size(); i++) {
					res[i] = 0;
				}
				
				int stride=system_dim;
				int iter=0;			
				for (int i=0;i<N;i++){
					F S=u[i*7];
					F E= u[i*7+1];
					F I=u[i*7+2];
					F R=u[i*7+3];
					F D=u[i*7+4];
					F C=u[i*7+5];
					F K=u[i*7+6];
					F sum=0;
					
					//S
					for (int j=0;j<N;j++){
						if (j!=i){
							sum+=sigma[j*N+i]*u[j*7+6];
						}
					}
					res[i*7*stride+i*7]=-eval_alpha(t,i)*(E+sum)-omega[i];
					res[i*7*stride+i*7+1]=-eval_alpha(t,i)*S;
								
					sum=0;
					for (int j=0;j<N;j++){
						if (j!=i){
							sum+=eval_alpha(t,i)*sigma[j*system_dim+i]*u[j*7+5];
						}
					}
					
					//res[i*7*stride+i*7+1]-=sum;
					//res[i*7*stride+i*7+5]=gamma[i];
	
					for (int j=0;j<N;j++){
						if (j!=i){
							res[i*7*stride+j*7+5]=-eval_alpha(t,i)*sigma[j*system_dim+i]*E;
							res[i*7*stride+j*7+6]=-eval_alpha(t,i)*S*sigma[j*system_dim+i];
						}
					}
					
					//E
					sum=0;
					for (int j=0;j<N;j++){
						if (j!=i){
							sum+=sigma[j*N+i]*u[j*7+6];
						}
					}
					
					res[(i*7+1)*stride+i*7]=eval_alpha(t,i)*(E+sum);
					sum=0;
					for (int j=0;j<N;j++){
						if (j!=i){
							sum=sigma[j*system_dim+i]*u[j*7+5];
						}
					}
					res[(i*7+1)*stride+i*7+1]=eval_alpha(t,i)*(S+sum)-phi[i]-(1/(qq[i]));
					
					res[(i*7+1)*stride+i*7+6]=epsilon[i];
					
					for (int j=0;j<N;j++){
						if (j!=i){
							res[(i*7+1)*stride+j*7+5]=eval_alpha(t,i)*sigma[j*system_dim+i]*E;
							res[(i*7+1)*stride+j*7+6]=eval_alpha(t,i)*sigma[j*system_dim+i]*S;
						}
					}			
					
					//I
					res[(i*7+2)*stride+i*7+1]=kappa[i]/qq[i];
					res[(i*7+2)*stride+i*7+2]=-1/pp[i];
					
					//R
					res[(i*7+3)*stride+i*7+1]=(1-kappa[i])/qq[i];
					res[(i*7+3)*stride+i*7+2]=(1-theta[i])/qq[i];
					
					//D
					res[(i*7+4)*stride+i*7+3]=theta[i]/qq[i];
					
					//C
				
					res[(i*7+5)*stride+i*7]=omega[i];
					sum=0;
					for (int j=0;j<N;j++){
						if (j!=i){
							sum=eval_alpha(t,j)*sigma[i*system_dim+j]*u[j*7+1];
						}
					}
					res[(i*7+5)*stride+i*7+5]=-gamma[i]-sum;
					
					for (int j=0;j<N;j++){
						if (j!=i){
							res[(i*7+5)*stride+j*7+1]=eval_alpha(t,i)*sigma[i*system_dim+j]*C;
						}
					}	
					
					//K				
					res[(i*7+6)*stride+i*7+1]=phi[i];
					
					sum=0;
					for (int j=0;j<N;j++){
						if (j!=i){
							sum=eval_alpha(t,j)*sigma[i*system_dim+j]*u[j*7+1];
						}
					}	
					
					for (int j=0;j<N;j++){
						if (j!=i){
							res[(i*7+5)*stride+j*7+5]=eval_alpha(t,j)*sigma[i*system_dim+j]*C;
						}
					}							
					
					res[(i*7+6)*stride+i*7+5]=sum;
					
					res[(i*7+6)*stride+i*7+6]=-epsilon[i];
					
					}
					
				return res;
			}
			std::tuple<std::vector<F>,std::vector<F>> run_linear_implicit(F t0, const T& u0, F tend, std::vector<std::vector<F>>* cumulated_exposed=nullptr, std::vector<std::vector<F>>* cumulated_infected=nullptr) {
				std::vector<F> u(7*N);
				for (int i=0;i<system_dim;i++){
					u[i]=u0[i];
				}
				
				for (int i=0;i<N;i++){
					cumulated_exposed_of_last_run[i]=u0[i*7+1];
					cumulated_infected_of_last_run[i]=u0[i*7+2];
					cumulated_exposed_of_last_run_container.push_back(std::vector<F>(1,u0[i*7+1]));
					cumulated_infected_of_last_run_container.push_back(std::vector<F>(1,u0[i*7+2]));	
				
				}
				utility::LinearImplicitSolver23<std::vector<F>,std::vector<F>,SEIRD_CK,F> solver(this,system_dim);
				solver.change_minimum_stepsize(ht);
				solver.change_maximum_stepsize(ht_max);
				solver.change_tol(tol);
				auto result=solver.run(t0, u, tend);
		
				if (cumulated_exposed != nullptr){
					*cumulated_exposed=cumulated_exposed_of_last_run_container;
					if (cumulated_infected != nullptr){
						*cumulated_infected=cumulated_infected_of_last_run_container;
					
					}
				}
				return std::make_tuple(result.first,result.second);
			}
			
			
			std::vector<std::vector<F>> get_cumulated_exposed() const{
				return cumulated_exposed_of_last_run_container;
			}
			std::vector<std::vector<F>> get_cumulated_infected() const{
				return cumulated_infected_of_last_run_container;
			}	
					
		};

	
	
	}
}

