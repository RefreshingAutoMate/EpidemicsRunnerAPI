#pragma once
#include "../../ConstrainedOptimization/core/parameters.h"
#include "../../ConstrainedOptimization/core/transformation.h"
#include "writer.h"
#include <utility>
#include <vector>

namespace ug {
	namespace epi {
		namespace utility {
			
			template<class T>
			class has_update_metainfo{
				typedef char yes[1];
				typedef char no[2];
				
				template<class C> 
				static yes& test(decltype(&C::update_metainfo));
				
				template<class C> 
				static no& test(...);
				
				public:
				static const bool value=sizeof(test<T>(0)) == sizeof(char);
				
			};
			
			
			
			template<class T1,class T2, class C,class F>
			class LinearImplicitSolver23 {

				C* model;
				int dim;
				F h_min=0.25;
				bool StoreToFile=false; //if true, results are stored to file
				OutputWriter<T1,F>* ow;
				std::string filepath="";
				std::string filename="output";	

				std::vector<F> create_identity_matrix() {
					std::vector<F>  result(dim*dim,F(0.0));
					for (int i = 0; i < dim; i++) {
						result[i + i * dim] = F(1.0);
					}
					return result;
				}

				std::vector<F> sumAB(F alpha, std::vector<F>& A, F beta, std::vector<F>& B) {
					std::vector<F>  result(dim * dim, F(0.0));

					for (int i = 0; i < dim;i++) {
						for (int j = 0; j < dim; j++) {
							result[i * dim + j] = alpha*A[i * dim + j] + beta*B[i * dim + j];
						}
					}
					return result;
				}
			bool class_has_metainfo=false;
			
			//double tol=0.0000000025;
			double tol=0.000025;				
			int p=2;
			
			double h_max=1;
			
			public:
				void set_store_to_file(bool _store_to_file, std::string _filepath, std::string _filename, OutputWriter<T1,F>* _ow){
					StoreToFile=_store_to_file;
					filepath=_filepath;
					filename=_filename;
					ow=_ow;
				}

				LinearImplicitSolver23(C* _model, int _dim) : model(_model), dim(_dim) {
					class_has_metainfo=has_update_metainfo<C>::value;
					//std::cout<<"Class has metainfo:"<<class_has_metainfo<<"\n";
				}
				
				void change_maximum_stepsize(F _h_max) {
					h_max=_h_max;
				}
				
				void change_minimum_stepsize(F _h_min) {
					h_min = _h_min;
				}
			
				std::pair<std::vector<F>,std::vector<F>> run_nonadaptive(F t0,  T1 u,  F tend) {
					if (u.size() != dim) {
						std::cerr << "Error: Input dimension of u0 in LinearImplicitSolver23 is different than the dimension previously given in the constructor\n";
						std::cerr << "u.size() = " << u.size() << std::endl;
						std::cerr << "dim = " << dim << std::endl;
					}
					if (StoreToFile){
						(*ow).write_to_file(filepath, filename+std::to_string(0)+".txt",t0,u);
					}
					std::vector<F> datapoints;

					for (int i = 0; i < dim; i++) {
						datapoints.push_back(u[i]);
					}
					std::vector<F> timepoints = {t0};
					std::vector<F> I;
					std::vector<F> J(dim * dim,F(0.0));
					std::vector<F> Qt(dim * dim);
					std::vector<F> R(dim * dim);
					std::vector<F> M(dim*dim); //M=(I-ahJ)
					std::vector<F> q1(dim);
					std::vector<F> k1(dim);
					std::vector<F> k2(dim);
					std::vector<F> k3(dim);
					
					F h = h_min;
					F t = t0;
					//std::copy(u0.begin(), u0.end(), u.begin());
					F a = 1 / (2 + 1.41);
			
					std::vector<F> u_copy(dim);
					

					int iter = 1;

					F h_prev=h;
					
					
					while (t<=tend+h) {		
						if constexpr (has_update_metainfo<C>::value){
							model->update_metainfo(u,t,h);
						}		

						std::copy(u.begin(), u.end(), u_copy.begin());
						//k1
						T1 temp = model->system(u,t); //it only has dim entries but otherwise type errors in this old version of dgemm
						std::vector<F> fy(dim);
						std::copy(temp.begin(),temp.end(), fy.begin());
						
						I = create_identity_matrix();
						T2 temp2=model->jacobian(u,t);
						std::copy(temp2.begin(), temp2.end(), J.begin());
						M=sumAB(F(1), I, F(-a * h), J);
						co::dc::qr<typename std::vector<F>::iterator,F> (M.begin(), dim, dim, Qt.begin(), R.begin());
						//co::dc::qr<typename std::vector<F>::iterator> (M.begin(), dim, dim, Qt.begin(), R.begin());
						co::mul::dgemm_nn(dim, 1, dim, F(1.0), Qt.begin(), 1, dim, fy.begin(), 1, 1, F(0.0), q1.begin(), 1, 1);
						co::dc::backwards_substitution<F>(R.begin(), k1.begin(), 1, q1.begin(), dim);
						
						//k2
						for (int i = 0; i < dim; i++) {
							u[i] = u_copy[i] + F(0.5) * h * k1[i];
							//std::cout<<I[i]<<"\t";
						}
					
						temp = model->system(u,t+F(0.5)*h); //it only has dim entries but otherwise type errors in this old version of dgemm
						temp2=model->jacobian(u,t+F(0.5)*h);
						std::copy(temp.begin(), temp.end(), fy.begin());
						co::mul::dgemm_nn(dim, 1, dim, -a*h, J.begin(), 1, dim, k1.begin(), 1, 1, F(1.0), fy.begin(), 1, 1);
						co::mul::dgemm_nn(dim, 1, dim, F(1.0), Qt.begin(), 1, dim, fy.begin(), 1, 1, F(0.0), q1.begin(), 1, 1);
						co::dc::backwards_substitution<F>(R.begin(), k2.begin(), 1, q1.begin(), dim);
						t+=h;
						for (int i = 0; i < dim; i++) {
									u[i] = u_copy[i] + h*k2[i];
						}	
						if (StoreToFile==false){
							for (int i = 0; i < dim; i++) {
								datapoints.push_back(u[i]);
							}
							timepoints.push_back(t);
						}
						else{
							(*ow).write_to_file(filepath, filename+std::to_string(iter)+".txt",t,u);
						}	
						iter++;
					}
					return std::make_pair(timepoints, datapoints);
				}
				std::pair<std::vector<F>,std::vector<F>> run(F t0,  T1 u,  F tend) {
					if (u.size() != dim) {
						std::cerr << "Error: Input dimension of u0 in LinearImplicitSolver23 is different than the dimension previously given in the constructor\n";
						std::cerr << "u.size() = " << u.size() << std::endl;
						std::cerr << "dim = " << dim << std::endl;
							return std::make_pair<std::vector<F>,std::vector<F>>({},{});
					}
					else if(h_min>=h_max){
						std::cerr<<"The solver stepsize hmin is greater or equal hmax\n";
						return std::make_pair<std::vector<F>,std::vector<F>>({},{});
					}
					if (StoreToFile){
						(*ow).write_to_file(filepath, filename+std::to_string(0)+".txt",t0,u);
					}
					std::vector<F> datapoints;

					for (int i = 0; i < dim; i++) {
						datapoints.push_back(u[i]);
					}
					std::vector<F> timepoints = {t0};
					std::vector<F> I;
					std::vector<F> J(dim * dim,F(0.0));
					std::vector<F> Qt(dim * dim);
					std::vector<F> R(dim * dim);
					std::vector<F> M(dim*dim); //M=(I-ahJ)
					std::vector<F> q1(dim);
					std::vector<F> k1(dim);
					std::vector<F> k2(dim);
					std::vector<F> k3(dim);


					F t = t0;
					//std::copy(u0.begin(), u0.end(), u.begin());
					F a = 1 / (2 + 1.41);
					F d31 = -(4 + 1.41) / (2 + 1.41);
					F d32 = (6+1.41)/ (2+1.41);
					std::vector<F> u_copy(dim);
					std::copy(u.begin(), u.end(), u_copy.begin());

					int iter = 1;
					F h =0.5*(h_max-h_min);

					bool stop_iterations=false;
					bool min_stepsize_reached=false;					
					double u_prev_norm1=0;
					for (int i=0;i<dim;i++){
						u_prev_norm1+=u[i]*u[i];
					}
					
					while (t<=tend) {
						stop_iterations=false;
						double err=999999999;
						double toln=tol*u_prev_norm1;
						int iter2=0;
						
						while (!stop_iterations){
							iter2++;
												
							//k1
							T1 temp = model->system(u,t); //it only has dim entries but otherwise type errors in this old version of dgemm
							std::vector<F> fy(dim);
							std::copy(temp.begin(),temp.end(), fy.begin());
							
							I = create_identity_matrix();
							T2 temp2=model->jacobian(u,t);
							std::copy(temp2.begin(), temp2.end(), J.begin());
							M=sumAB(F(1), I, F(-a * h), J);
							co::dc::qr<typename std::vector<F>::iterator,F> (M.begin(), dim, dim, Qt.begin(), R.begin());
							//co::dc::qr<typename std::vector<F>::iterator> (M.begin(), dim, dim, Qt.begin(), R.begin());
							co::mul::dgemm_nn(dim, 1, dim, F(1.0), Qt.begin(), 1, dim, fy.begin(), 1, 1, F(0.0), q1.begin(), 1, 1);
							co::dc::backwards_substitution<F>(R.begin(), k1.begin(), 1, q1.begin(), dim);
							
							//k2
							for (int i = 0; i < dim; i++) {
								u[i] = u_copy[i] + F(0.5) * h * k1[i];
								//std::cout<<I[i]<<"\t";
							}
						
							temp = model->system(u,t+F(0.5)*h); //it only has dim entries but otherwise type errors in this old version of dgemm
							temp2=model->jacobian(u,t+F(0.5)*h);
							std::copy(temp.begin(), temp.end(), fy.begin());
							co::mul::dgemm_nn(dim, 1, dim, -a*h, J.begin(), 1, dim, k1.begin(), 1, 1, F(1.0), fy.begin(), 1, 1);
							co::mul::dgemm_nn(dim, 1, dim, F(1.0), Qt.begin(), 1, dim, fy.begin(), 1, 1, F(0.0), q1.begin(), 1, 1);
							co::dc::backwards_substitution<F>(R.begin(), k2.begin(), 1, q1.begin(), dim);

							//k3 (useful for adaptive stepsizes)
							
							
							//std::cout<<"toln"<<toln<<"\n";
							err=0;
							for (int i = 0; i < dim; i++) {
			
									u[i] = u_copy[i] + h*k2[i];
									F diff = h*k2[i] -((h / 6) * (k1[i] + 4 * k2[i] + k3[i]));
									err+=diff*diff;
									
									//std::cout<<u[i]<<" vs." <<u_copy[i]<<"\n";
							}	
						//err=1e-15
							if(err==0){
							err+=1e-30;
								
							}
					
							/*
							if (iter2>20){
							std::cout<<"err: "<<err<<"\n";	
							std::cout<<"h: "<<h<<"\n";
							std::cout<<"hnext: "<<h*std::pow(toln/err,1.0/(p+1))<<"\n";
							std::cout<<"toln: "<<toln<<"\n";							
							std::cout<<"t: "<<t<<"\n";
							std::cout<<"iter2: "<<iter2<<"\n";
							std::cout<<((toln/err))<<"\n";
							std::cout<<((toln/err)>=1)<<"\n";
							std::cin.get();
							}*/
						//	std::cout<<h<<"\n";
			
							if ((toln/err)>=1 || min_stepsize_reached){
								if constexpr (has_update_metainfo<C>::value){
									model->update_metainfo(u,t,h);
								}	
								
								std::copy(u.begin(), u.end(), u_copy.begin());
							//		std::cout<<"err: "<<err<<"\n";
						//	td::cout<<"h: "<<h<<"\n";
						//	td::cin.get();
							//	std::cout<<"drin\n";
							//	std::cout<<"t: "<<t<<"\n";
								//calculating yi+1	
								t+=h;
								if (StoreToFile==false){
									for (int i = 0; i < dim; i++) {
										datapoints.push_back(u[i]);
									}
									timepoints.push_back(t);
								}
								else{
									(*ow).write_to_file(filepath, filename+std::to_string(iter)+".txt",t,u);
								}	
								min_stepsize_reached=false;
								stop_iterations=true;
								
								u_prev_norm1=0;
								for (int i=0;i<dim;i++){
									u_prev_norm1+=u[i]*u[i];
								}
								
								//h_prev=h;
							}
							else{
								std::copy(u_copy.begin(), u_copy.end(), u.begin());
							}
							
													
							h=h*std::pow(toln/err,1.0/(p+1));
							
							
							h=(h>=h_max)?h_max:h;
							
							
							//std::cout<<"err: "<<err<<"\n";
							//std::cout<<"h: "<<h<<"\n";
							//std::cout<<"t: "<<t<<"\n";
							
							
							
							if (h<h_min || iter2>30){
								//std::cout<<"hmin: "<<h<<"   t: "<<t<<"\n";
								h=h_min;
								min_stepsize_reached=true;
								//std::cout<<"min reached: "<<h<<"\n";
							}
							
							if (t+h>=tend){
								h=tend-t;
								
							}
						}
						if (t>=tend){
							break;
						}								
					//	std::cout<<"h: "<<h<<" iter "<<t<<"\n";			
						iter++;
						//t+=h;
					}

					return std::make_pair(timepoints, datapoints);
				}
				
			void change_tol(double _tol){
				tol=_tol;
			}
			
			};
		}
	}

}


