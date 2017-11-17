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



	const int length_x = 200; // length of the chemoattractant vector
	const int length_y = 50;


	
	const double diameter = 1; // diameter in which there have to be no cells



	int cell_size = 2; // cell size relative to mesh

	MatrixXf chemo(length_x, length_y);	

	// generate gradient which linearly increases
	for (int i = 0;i<length_x;i++){
		for (int j = 0; j< length_y; j++){
			chemo(i,j) = i; // this time it increases as down the rows, but this is consistent with the movement defined afterwards, we initial set cells on x=1, varying y
		}			
	}


	  
	/*for (int i =0;i<length;i++){
		for(int j = 0;j<length;j++){
        	output << chemo(i,j) << ", "; 
		}
		output << "\n" << std::endl; 
    	}*/

	// three columns for x, y, z
	
	// form a matrix which would store x,y,z

	/* MatrixXf chemo_3col(length_x*length_y,4);

	// x, y coord, 1st and 2nd columns respectively
	int k = 0;
		// it has to be 3D for paraview
	while (k<length_x*length_y){
		for (int i = 0;i<length_x;i++){
			for (int j = 0; j< length_y; j++){
				chemo_3col(k,0) = i; 
				chemo_3col(k,1) = j;
				chemo_3col(k,2) = 0;
				k += 1;
			}			
		}
	}

	// u column
	for (int i=0;i<length_x*length_y;i++){
		chemo_3col(i,3) = chemo(chemo_3col(i,0),chemo_3col(i,1));
	}

	// save data to plot chemoattractant concentration in MATLAB
	ofstream output("matrix_3col.csv");
	
		output << "x, y, z, u" << "\n" << endl;


	for (int i=0;i<length_x*length_y;i++){
		for(int j=0;j<4;j++){
			output << chemo_3col(i,j) << ", ";
		}
		output << "\n" << endl;
	}	*/

	

	    /*
	     * 2D domain with a few randomly placed particles
	     */

	/*
	 * initial cells	
	 */

	const size_t N = 50;
	//ABORIA_VARIABLE(velocity,vdouble2,"velocity")
	typedef Particles<std::tuple<>, 2> particle_type;
	//typedef Particles<std::tuple<>,2,std::vector,bucket_search_serial> particle_type;
	typedef particle_type::position position;
	particle_type particles(N);
	std::default_random_engine gen;
	std::uniform_real_distribution<double> uniform(1,length_y);  
	for (int i=0; i<N; ++i) {
	    get<position>(particles)[i] = vdouble2(4,uniform(gen)); // x=2, uniformly in y
	
	}

	// save particles before they move
	vtkWriteGrid("before",0,particles.get_grid(true));

	/*
         * initialise neighbour search with 2d cuboid domain,
         * periodic in x and y
         */

        particles.init_neighbour_search(vdouble2(0,0),vdouble2(length_x,length_y),vbool2(false,false));

        
	

	// Update positions based on the gradient

	int N_steps = 300; // number of times the cells move up the gradient
	

	// choose a set of random number between 0 and 2pi
		std::default_random_engine gen1;
		std::uniform_real_distribution<double> uniformpi(0,M_PI); // can only move forward
		VectorXf random_angle(N_steps*particles.size());  

		for (int i = 0; i<N_steps*particles.size();i++){
			random_angle(i) = uniformpi(gen1);		
			//cout << "angle to move " << random_angle(i) << endl;
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
		
		// make sure that I choose the right direction when taking into account cell radius
		int sign_x, sign_y;
		if(sin(random_angle(rand_num_count))<0){
			sign_x=-1;
		}else{sign_x=1;}

		if(cos(random_angle(rand_num_count))<0){
			sign_y=-1;
		}else{sign_y=1;}
			
		//cout << "sin of the angle " << sin(random_angle(rand_num_count)) << endl;
		//cout << "cos of the angle " << cos(random_angle(rand_num_count)) << endl;

		
		// check if the x coordinates are not out of the domain, if they are, ignore that step
		
		//cout << "chemo coord x before the test " << round(x[0]+sin(random_angle(rand_num_count))+sign_x*cell_size) << endl;
		//cout << "chemo coord y before the test " << round(x[0]+cos(random_angle(rand_num_count))+sign_y*cell_size) << endl;

		if (round(x[0]+sin(random_angle(rand_num_count))+sign_x*cell_size)>-1 && round(x[0]+sin(random_angle(rand_num_count))+sign_x*cell_size)<length_x && round(x[1]+ cos(random_angle(rand_num_count))+sign_y*cell_size) >-1 && round(x[1]+ cos(random_angle(rand_num_count))+sign_y*cell_size)<length_y ){

		//cout << "sin of the angle (inside the loop) " << sin(random_angle(rand_num_count)) << endl;
		//cout << "cos of the angle (inside the loop) " << cos(random_angle(rand_num_count)) << endl;

		// check if the gradient in the other position is larger, if yes, move to that position, x changes by sin and y to cos, because of the the chemo is defined. 

			if (chemo(round(x)[0],round(x)[1]) < chemo(round(x[0]+sin(random_angle(rand_num_count))+sign_x*cell_size),round(x[1]+ cos(random_angle(rand_num_count))+sign_y*cell_size))){
				//cout << "x coord " << x[0] << endl;
				//cout << "x up " << sin(random_angle(rand_num_count))+sign_x*cell_size << endl;
				//cout << "y coord " << x[1] << endl;
				//cout << "y up " << cos(random_angle(rand_num_count))+sign_y*cell_size << endl;	
				//cout << "chemo coonc in current site " << chemo(round(x)[0],round(x)[1])<< endl;	
				//cout << "chemo coonc in other site " << chemo(round(x[0] +sin(random_angle(rand_num_count)) +sign_x*cell_size),round(x[1]+cos(random_angle(rand_num_count))+sign_y*cell_size))<< endl;	
				//get<position>(particles)[i] += vdouble2(sin(random_angle(rand_num_count))+sign_x*cell_size, cos(random_angle(rand_num_count))+sign_y*cell_size);

	// check if there are no cells around the position where I want to move

                /*
                 * loop over all neighbouring particles within a euclidean distance
                 * of size "diameter"
                 */
		x += vdouble2(sin(random_angle(rand_num_count)),cos(random_angle(rand_num_count)));

		//cout << "Position "<< x << endl;
		int count_position = 0;
		bool free_position = true; // check if the neighbouring position is free

		// if this loop is entered, it means taht there is another cell where I want to move 
                for (auto tpl: euclidean_search(particles.get_query(),x,diameter)) {

			count_position += 1; // just to check if this works

                    free_position = true;
                    //break;
                }
		//cout << "print position " << count_position << endl;
		if (free_position == true){
			get<position>(particles)[i] += vdouble2(sin(random_angle(rand_num_count)), cos(random_angle(rand_num_count)));
		}

			}
		}
			rand_num_count += 1; // update random number count			
	}
	//particles.update_positions();
		// save particles after they move
			    /*
		     * on every i/o step write particle container to a vtk
		     * unstructured grid file
		     */
		    std::cout << "." << std::flush;
	#ifdef HAVE_VTK
	vtkWriteGrid("particles",j,particles.get_grid(true));
    
	#endif

	//for (int i =0;i<5;++i){cout << "koks rezas " << i <<endl;}
	
}


}
