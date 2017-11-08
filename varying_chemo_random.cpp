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
	int cell_size = 2; // cell size relative to mesh

	// create matrices for the current value and the updated one
	MatrixXf chemo(lattice_size,lattice_size), chemo_new(lattice_size,lattice_size); 	

	// initially gradient is uniform
	for (int i =0;i<lattice_size;i++){
		for (int j=0;j<lattice_size;j++){
			chemo(i,j)=1;
			chemo_new(i,j)=0;// this is for later updates, for now only interior, so set it to 0
		}
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
	double dt = 1; // time step, redundant
	int t_final = 50; // final time	
	int dx = 1; // space step in x direction
	int dy = 1; // space step in y direction
	double kai = 0.0001; // to 1 /h production rate of chemoattractant
	int step = 1;


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

	const size_t N = 5;
	//ABORIA_VARIABLE(velocity,vdouble2,"velocity")
	typedef Particles<std::tuple<>, 2> particle_type;
	//typedef Particles<std::tuple<>,2,std::vector,bucket_search_serial> particle_type;
	typedef particle_type::position position;
	particle_type particles(N);
	std::default_random_engine gen;
	std::uniform_real_distribution<double> uniform(0,lattice_size);  
	for (int i=0; i<N; ++i) {
	    get<position>(particles)[i] = vdouble2(4,uniform(gen)); // x=2, uniformly in y
	
	}

	// save particles before they move
	vtkWriteGrid("before",0,particles.get_grid(true));

	// Update positions based on the gradient

	int N_steps = 1; // number of times the cells move up the gradient

	// choose a set of random number between 0 and 2pi
		std::default_random_engine gen1;
		std::uniform_real_distribution<double> uniformpi(0,2*M_PI);
		VectorXf random_angle(N_steps*particles.size());  

		for (int i = 0; i<N_steps*particles.size();i++){
			random_angle(i) = uniformpi(gen1);		
			cout << "angle to move " << random_angle(i) << endl;
		}

	int rand_num_count =0;


	for (int t=0;t<t_final;t++){
		

		// internalisation
		for (int i =0;i<lattice_size;i++){
			for(int j = 0;j<lattice_size;j++){
				//go thorugh all the cells
				for (int k =0; k<particles.size();k++){
					vdouble2 x;
					x = get<position>(particles[k]);
					intern(i,j) = intern(i,j) + exp(- (domain_len*domain_len*(i-x[0]*mesh_per_domain)*(i-x[0]*mesh_per_domain)+(j-x[1]*mesh_per_domain)*(j-x[1]*mesh_per_domain))/(2*R*R));
				}			
			}
    		}		

		
		// internal part of the chemoattractant, finite difference for reaction diffusion equation
		for (int i=1;i<lattice_size-1;i++){
			for (int j=1;j<lattice_size-1;j++){
				chemo_new(i,j) = dt * (D*((1/(domain_len*domain_len))* (chemo(i+1,j)-2*chemo(i,j)+chemo(i-1,j))/(dx*dx) + (chemo(i,j+1)- 2* chemo(i,j)+chemo(i,j-1))/(dy*dy)  ) - (chemo(i,j)*lam / (2*M_PI*R*R)) * intern(i,j) + kai*chemo(i,j)*(1-chemo(i,j)) - domain_len_der/domain_len *chemo(i,j) ) + chemo(i,j);
			cout << "print the internalisation term " << intern(i,j) << endl;
			cout << "new chemo " << chemo_new(i,j) << endl;
			cout << "chemo " << chemo(i,j) << endl;
			}
		}
		

		// change only the interior
		/*for (int i=1;i<lattice_size-1;i++){
			for (int j=1;j<lattice_size-1;j++){
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
		
				//cout << "x coord " << x[0] << endl;	


			cout << " code stops here"<< sin(random_angle(rand_num_count) << endl;
			// check if chemoattractant concentration is greater at the random direction we selected
			//sign(cos(random_angle(rand_num_count)))
			//sign(sin(random_angle(rand_num_count)))
		
			// make sure that I choose the right direction when taking into account cell radius
			int sign_x, sign_y;
			if(sin(random_angle(rand_num_count))<0){
				sign_x=-1;
			}else{sign_x=1;}

			if(cos(random_angle(rand_num_count))<0){
				sign_y=-1;
			}else{sign_y=1;}
				

			// check if the gradient in the other position is larger, if yes, move to that position, x changes by sin and y to cos, because of the the chemo is defined. 

			cout << "xcoord  "<< round(x)[0] << endl;
			cout << "xcoord new "<< round(x)[1] << endl;
			cout << "ycoord "<< round(x[0]+sin(random_angle(rand_num_count))+sign_x*cell_size) << endl;
			cout << "ycoord new " << round(x[1]+ cos(random_angle(rand_num_count))+sign_y*cell_size) << endl;

			cout << "chemo at current position "<< chemo(round(x)[0],round(x)[1]) << endl;

			cout << "chemo at other position "<< chemo(round(x[0]+sin(random_angle(rand_num_count))+sign_x*cell_size),round(x[1]+ cos(random_angle(rand_num_count))+sign_y*cell_size)) << endl;

			if (chemo(round(x)[0],round(x)[1]) < chemo(round(x[0]+sin(random_angle(rand_num_count))+sign_x*cell_size),round(x[1]+ cos(random_angle(rand_num_count))+sign_y*cell_size))){
				cout << "x coord " << x[0] << endl;
				cout << "x up " << sin(random_angle(rand_num_count))+sign_x*cell_size << endl;
				cout << "y coord " << x[1] << endl;
				cout << "y up " << cos(random_angle(rand_num_count))+sign_y*cell_size << endl;	
				cout << "chemo coonc in current site " << chemo(round(x)[0],round(x)[1])<< endl;	
				cout << "chemo coonc in other site " << chemo(round(x[0] +sin(random_angle(rand_num_count)) +sign_x*cell_size),round(x[1]+cos(random_angle(rand_num_count))+sign_y*cell_size))<< endl;	
				//get<position>(particles)[i] += vdouble2(sin(random_angle(rand_num_count))+sign_x*cell_size, cos(random_angle(rand_num_count))+sign_y*cell_size);
				get<position>(particles)[i] += vdouble2(sin(random_angle(rand_num_count)), cos(random_angle(rand_num_count)));
			}

				rand_num_count += 1; // update random number count
			
		}
		


	}// end of for loop with t<t_final

	
// save particles after they move
vtkWriteGrid("after",0,particles.get_grid(true));

// save data for final chemoattractant concentration
ofstream output("final_chemo.txt");
	  
	for (int i =0;i<lattice_size;i++){
		for(int j = 0;j<lattice_size;j++){
        	output << chemo(i,j) << " "; 
		}
		output << "\n" << std::endl; 
    	}	

}
