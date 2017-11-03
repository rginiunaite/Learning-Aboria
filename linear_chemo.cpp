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
			chemo(i,j) = i; // down the rows increases
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

	const size_t N = 5;
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

	int N_steps = 30; // number of times the cells move up the gradient
	int step = 1;

for (int j=0; j<N_steps;j++){
	for (int i=0; i < particles.size(); i++){


			vdouble2 x;
			x = get<position>(particles[i]);
		
			//cout << "x coord " << x[0] << endl;

			//cout << "y coord " << x[1] << endl;	
		// choose internal positions
		if (round(x)[0]>1 && round(x)[0] < length && round(x)[1]>1 && round(x)[1] < length ){
		
		// check which of the four directions we obtain the biggest gradient

					VectorXf directions(4);

					double up = chemo(round(x)[0]+cell_size,round(x)[1])- chemo(round(x)[0],round(x)[1]);		

			cout << "up " << up << endl;
					double down = chemo(round(x)[0]-cell_size,round(x)[1])- chemo(round(x)[0],round(x)[1]);
			cout << "down " << down << endl;
					double right = chemo(round(x)[0],round(x)[1]+cell_size)- chemo(round(x)[0],round(x)[1]);
			cout << "right " << right << endl;

					double left = chemo(round(x)[0],round(x)[1]-cell_size)- chemo(round(x)[0],round(x)[1]);

			cout << "left " << left << endl;

					directions(0) = up;
					directions(1) = down;
					directions(2) = right;
					directions(3) = left;

		
					// find maximum direction, if maximum is not unique, it chooses the one with the lowest coefficient, I will have to change that then it chooses randomly
					Eigen::VectorXf::Index max_index;
					double max_dir_index = directions.maxCoeff(&max_index);// do not forget from 0 to 3
					//cout << "max index " << max_index << endl;
					double change;

			// do not move if maximum is negative

			if (directions(max_index)<=0)
			{
				get<position>(particles)[i] += vdouble2(0,0);
			}
			// if not, then move up the gradient
			else{
			//cout << "old position " << get<position>(particles)[i] << endl;
					if (max_index == 0)// up
					{
						get<position>(particles)[i] += vdouble2(step,0);
					}
					if (max_index == 1)// down
					{
						get<position>(particles)[i] += vdouble2(-step,0);
					}
					if (max_index == 2)// right
					{
						get<position>(particles)[i] += vdouble2(0,step);
					}
					if (max_index == 3)// left
					{
						get<position>(particles)[i] += vdouble2(0,-step);
					}
			//cout << "new position " << get<position>(particles)[i] << endl;
			}
		}
	
	}
	//particles.update_positions();

}

	// save particles after they move
	vtkWriteGrid("fixed",0,particles.get_grid(true)); 

}
