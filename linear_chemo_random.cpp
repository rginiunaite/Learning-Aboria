// cells move towords high concentration of fixed gradient

#include "Aboria.h"
#include <random>
#include <map>
#include <Eigen/Core>
#include <algorithm> // std::unique_copy
// writing on a text file
#include <iostream>
#include <fstream>


#include <math.h>

// just to check if changes work

using namespace std;
using namespace Aboria;
using namespace Eigen; // objects VectorXf, MatrixXf

int main() {



	const int length = 100; // length of the chemoattractant vector

	int cell_size = 2; // cell size relative to mesh

	MatrixXf chemo(length, length);	

	// generate gradient which linearly increases
	for (int i = 0;i<length;i++){
		for (int j = 0; j< length; j++){
			chemo(i,j) = i*i;
		}			
	}



	
	// save data to plot chemoattractant concentration in MATLAB
	ofstream output("matrix.txt");
	  
	for (int i =0;i<length;i++){
		for(int j = 0;j<length;j++){
        	output << chemo(i,j) << " "; 
		}
		output << "\n" << std::endl; 
    	}	

	

	    /*
	     * 2D domain with a few randomly placed particles
	     */

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
	std::uniform_real_distribution<double> uniform(1,length);  
	for (int i=0; i<N; ++i) {
	    get<position>(particles)[i] = vdouble2(4,uniform(gen)); // x=2, uniformly in y
	
	}

	// save particles before they move
	vtkWriteGrid("before",0,particles.get_grid(true));

	
	

	// Update positions based on the gradient

	int N_steps = 4; // number of times the cells move up the gradient
	

	// choose a set of random number between 0 and 2pi
		std::default_random_engine gen1;
		std::uniform_real_distribution<double> uniformpi(0,2*M_PI);
		VectorXf random_angle(N_steps*particles.size());  

		for (int i = 0; i<N_steps*particles.size();i++){
			random_angle(i) = uniformpi(gen1);		
			cout << "angle to move " << random_angle(i) << endl;
		}

	int rand_num_count =0;

for (int j=0; j<N_steps;j++){
	for (int i=0; i < particles.size(); i++){
		

		vdouble2 x;
		x = get<position>(particles[i]);
		
			//cout << "x coord " << x[0] << endl;	



		// check if chemoattractant concentration is greater at the random direction we selected
		//sign(cos(random_angle(rand_num_count)))
		//sign(sin(random_angle(rand_num_count)))
		
		// make sure that I choose the right direction
		int sign_x, sign_y;
		if(cos(random_angle(rand_num_count))<0){
			sign_x=-1;
		}else{sign_x=1;}

		if(sin(random_angle(rand_num_count))<0){
			sign_y=-1;
		}else{sign_y=1;}
			

		if (chemo(round(x)[0],round(x)[1]) < chemo(round(x[0] + cos(random_angle(rand_num_count))+sign_x*cell_size),round(x[1]+sin(random_angle(rand_num_count))+sign_y*cell_size))){
			cout << "x coord " << x[0] << endl;
			cout << "x up " << cos(random_angle(rand_num_count))+sign_x*cell_size << endl;
			cout << "y coord " << x[1] << endl;
			cout << "y up " << sin(random_angle(rand_num_count))+sign_y*cell_size << endl;	
			cout << "chemo coonc in current site " << chemo(round(x)[0],round(x)[1])<< endl;	
			cout << "chemo coonc in other site " << chemo(round(x[0] + cos(random_angle(rand_num_count))+sign_x*cell_size),round(x[1]+sin(random_angle(rand_num_count))+sign_y*cell_size))<< endl;	
			get<position>(particles)[i] += vdouble2(cos(random_angle(rand_num_count))+sign_x*cell_size,sin(random_angle(rand_num_count))+sign_y*cell_size);
		}

			rand_num_count += 1; // update random number count
			
	}
	//particles.update_positions();

}

	// save particles after they move
	vtkWriteGrid("fixed",0,particles.get_grid(true)); 

}
