// varying concentration of chemoattractant, reaction diffusion equation for chemoattractant

#include "Aboria.h"
#include <random>
#include <map>
#include <Eigen/Core>
#include <algorithm> // std::unique_copy
// writing on a text file
#include <iostream>
#include <fstream>

#include <cmath> // for power

#include <math.h>

// just to check if changes work

using namespace std;
using namespace Aboria;
using namespace Eigen; // objects VectorXf, MatrixXf

int main() {

	
	int lattice_size = 100;

	// create matrices for the current value and the updated one
	MatrixXf chemo(lattice_size,lattice_size), chemo_new(lattice_size,lattice_size); 	

	// initially gradient is uniform internally
	for (int i =1;i<lattice_size-1;i++){
		for (int j=1;j<lattice_size-1;j++){
			chemo(i,j)=1;
			chemo_new(i,j)=1;// this is for later updates, for now only interior, so set it to 0
		}
	}

	// there is no chemoattractant on the boundary, to avoid all the cells gathering there, chemoattractant there is consumed at a slower rate than anywthere else

	for (int i =0;i<lattice_size;i++){	
		chemo(i,0)=0;
		chemo_new(i,0)=0;// this is for later updates, for now only interior, so set it to 0
	}

	for (int i =0;i<lattice_size;i++){	
		chemo(0,i)=0;
		chemo_new(0,i)=0;// this is for later updates, for now only interior, so set it to 0
	}

	for (int i =0;i<lattice_size;i++){	
		chemo(lattice_size -1,i)=0;
		chemo_new(lattice_size -1,i)=0;// this is for later updates, for now only interior, so set it to 0
	}

	for (int i =0;i<lattice_size;i++){	
		chemo(i,lattice_size -1)=0;
		chemo_new(i,lattice_size -1)=0;// this is for later updates, for now only interior, so set it to 0
	}



	// save data to plot chemoattractant concentration in MATLAB
	/*ofstream output("varying_chemo.txt");
	  
	for (int i =0;i<lattice_size;i++){
		for(int j = 0;j<lattice_size;j++){
        	output << chemo(i,j) << " "; 
		}
		output << "\n" << std::endl; 
    	}*/


	/*
	 * update changes in chemoattractant concentration based on reaction-diffusion equation	
	 */

	// parameters
	double D = 0.1; // to 10^5 \nu m^2/h diffusion coefficient
	double t = 0; // initialise time, redundant
	double dt =1; // time step, redundant
	int t_final = 60; // final time	
	int dx = 1; // space step in x direction
	int dy = 1; // space step in y direction
	double kai = 0.0001; // to 1 /h production rate of chemoattractant



	// initially assume that the domain length is fixed
	
	int domain_len = 100;
	double domain_len_der = 0; // for now this is useless

	int mesh_per_domain = lattice_size/domain_len; // the ratio between total mesh points and domain length

	// parameters for internalisation

	double R = 7.5; // \nu m cell radius
	int lam = 100; // to 1000 /h chemoattractant internalisation
	
	 

	// initialise internalisation matrix
	MatrixXf intern(lattice_size,lattice_size);

	/*
	 * initial cells	
	 */

	const size_t N = 10;
	//ABORIA_VARIABLE(velocity,vdouble2,"velocity")
	typedef Particles<std::tuple<>, 2> particle_type;
	//typedef Particles<std::tuple<>,2,std::vector,bucket_search_serial> particle_type;
	typedef particle_type::position position;
	particle_type particles(N);
	std::default_random_engine gen;
	std::uniform_real_distribution<double> uniform(1,domain_len);  
	for (int i=0; i<N; ++i) {
	    get<position>(particles)[i] = vdouble2(1,uniform(gen)); // x=2, uniformly in y
	
	}

	// save particles before they move
	vtkWriteGrid("before",0,particles.get_grid(true));


	while (t < t_final){
		

		// internalisation
		for (int i =0;i<lattice_size;i++){
			for(int j = 0;j<lattice_size;j++){
				//go thorugh all the cells
				for (int k =0; k<particles.size();k++){
					vdouble2 x;
					x = get<position>(particles[k]);
					intern(i,j) = intern(i,j) + exp(- (domain_len*domain_len*(i-x[0]*mesh_per_domain)*(i-x[0]*mesh_per_domain)+(j-x[1]*mesh_per_domain)*(j-x[1]*mesh_per_domain))/(2*R*R));
					//cout << "xcoord" << x[0] << endl;
					//cout << "ycoord" << x[1] << endl;
				}			
			}
    		}		

		
		// internal part of the chemoattractant, finite difference for reaction diffusion equation
		for (int i=1;i<lattice_size-1;i++){
			for (int j=1;j<lattice_size-1;j++){
				chemo_new(i,j) = dt * (D*((1/(domain_len*domain_len))* (chemo(i+1,j)-2*chemo(i,j)+chemo(i-1,j))/(dx*dx) + (chemo(i,j+1)- 2* chemo(i,j)+chemo(i,j-1))/(dy*dy)  ) - (chemo(i,j)*lam / (2*M_PI*R*R)) * intern(i,j) + kai*chemo(i,j)*(1-chemo(i,j)) - domain_len_der/domain_len *chemo(i,j) ) + chemo(i,j);
			//cout << "print the internalisation term " << intern(i,j) << endl;
			//cout << "new chemo " << chemo_new(i,j) << endl;
			//cout << "chemo " << chemo(i,j) << endl;
			}
		}
		

		// change only the interior
		/*for (int i=0;i<lattice_size;i++){
			for (int j=0;j<lattice_size;j++){
				chemo(i,j) = chemo_new(i,j);
			}
		}*/
		chemo = chemo_new;
		
		/*
		 * cells move 1 step up the cell-induced gradient
		 */


		for (int i=0; i < particles.size(); i++){


			vdouble2 x;
			x = get<position>(particles[i]);
			//cout << "print coord of chemo x1 "<< round(x)[0]*mesh_per_domain+1 << endl;
			//cout << "print coord of chemo y1 "<< (round(x)[1])*mesh_per_domain +1 << endl;
			//cout << "print coord of chemo x "<< round(x)[0]*mesh_per_domain-1 << endl;
			//cout << "print coord of chemo y "<< (round(x)[1])*mesh_per_domain -1 << endl;
			
			
			// approximate coordinates

			int xplus1 = round(x)[0]*mesh_per_domain+1;
			int xminus1 = round(x)[0]*mesh_per_domain-1;
			int xeq = round(x)[0]*mesh_per_domain;
			int yplus1 = round(x)[1]*mesh_per_domain+1;
			int yminus1 = round(x)[1]*mesh_per_domain-1;
			int yeq = round(x)[1]*mesh_per_domain;

			if (xplus1<0){xplus1=0;}
			else if(xplus1>=lattice_size){xplus1=lattice_size -1;} 
			if (xminus1<0){xminus1=0;}
			else if(xminus1>=lattice_size){xminus1=lattice_size -1;}
			if (xeq<0){xeq=0;}
			else if(xeq>=lattice_size){xeq=lattice_size -1;}
			if (yplus1<0){yplus1=0;}
			else if(yplus1>=lattice_size){yplus1=lattice_size -1;}
			if (yminus1<0){yminus1=0;}
			else if(yminus1>=lattice_size){yminus1=lattice_size -1;}
			if (yeq<0){yeq=0;}
			else if(yeq>=lattice_size){yeq=lattice_size -1;}
			


			//double change_x = (chemo(xplus1,yeq)-chemo(xminus1,yeq))/(chemo(xeq,yeq)); relative

			//double change_y = (chemo(xeq,yplus1)-chemo(xeq,yminus1))/(chemo(xeq,yeq)); relative


			double change_x = (chemo(xplus1,yeq)-chemo(xminus1,yeq));

			double change_y = (chemo(xeq,yplus1)-chemo(xeq,yminus1)); 


			cout << "step size in x direction " << change_x << endl;

			cout << "step size in y direction " << change_y << endl;
 					
			get<position>(particles)[i] += vdouble2(change_x,change_y);

			// do not move if maximum is negative

	
		}//end of < particle_size
		//particles.update_positions();
		
	t = t +dt;

	}// end of for loop with t<t_final

	
// save particles after they move
	vtkWriteGrid("after",0,particles.get_grid(true));

// save data for final chemoattractant concentration
ofstream output("final_chemo.txt");
	  
	for (int i =0;i<lattice_size;i++){
		for(int j = 0;j<lattice_size;j++){
        	output << chemo(i,j) << " "; 
		//cout << chemo(i,j) << " ";
		}
		//cout << "\n " << endl;
		output << "\n" << std::endl; 
    	}	

}
