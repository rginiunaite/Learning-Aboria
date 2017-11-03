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
			chemo(i,j) = i; // this time it increases as down the rows, but this is consistent with the movement defined afterwards, we initial set cells on x=1, varying y
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
	

	// choose a set of random number between 0 and 2pi
		std::default_random_engine gen1;
		std::uniform_real_distribution<double> uniformpi(0,2*M_PI);
		VectorXf random_angle(N_steps*particles.size());  

		for (int i = 0; i<N_steps*particles.size();i++){
			random_angle(i) = uniformpi(gen1);		
			cout << "angle to move " << random_angle(i) << endl;
		}

	int rand_num_count = 0;

for (int j=0; j<N_steps;j++){
	for (int i=0; i < particles.size(); i++){
		

		vdouble2 x;
		x = get<position>(particles[i]);
		
			//cout << "x coord " << x[0] << endl;		

			// approximate coordinates

			int xplus1 = round(x)[0]+1;
			int xminus1 = round(x)[0]-1;
			int xeq = round(x)[0];
			int yplus1 = round(x)[1]+1;
			int yminus1 = round(x)[1]-1;
			int yeq = round(x)[1];

			//cout << "chemo x coord plus before " << xplus1 << endl;
			//cout << "chemo x coord minus before " << xminus1 << endl;
			//cout << "chemo y coord plus before " << yplus1 << endl;
			//cout << "chemo y coord minus before " << yminus1 << endl;

			if (xplus1<0){xplus1=0;}
			else if(xplus1>=length){xplus1=length -1;} 
			if (xminus1<0){xminus1=0;}
			else if(xminus1>=length){xminus1=length -1;}
			if (xeq<0){xeq=0;}
			else if(xeq>=length){xeq=length -1;}
			if (yplus1<0){yplus1=0;}
			else if(yplus1>=length){yplus1=length -1;}
			if (yminus1<0){yminus1=0;}
			else if(yminus1>=length){yminus1=length -1;}
			if (yeq<0){yeq=0;}
			else if(yeq>=length){yeq=length -1;}
			


			//double change_x = (chemo(xplus1,yeq)-chemo(xminus1,yeq))/(chemo(xeq,yeq)); relative

			//double change_y = (chemo(xeq,yplus1)-chemo(xeq,yminus1))/(chemo(xeq,yeq)); relative
			//cout << "chemo x coord plus  " << xplus1 << endl;
			//cout << "chemo x coord minus  " << xminus1 << endl;
			//cout << "chemo y coord plus " << yplus1 << endl;
			//cout << "chemo y coord minus  " << yminus1 << endl;


			double change_x = (chemo(xplus1,yeq)-chemo(xminus1,yeq));

			double change_y = (chemo(xeq,yplus1)-chemo(xeq,yminus1)); 


			//cout << "step size in x direction " << change_x << endl;

			//cout << "step size in y direction " << change_y << endl;

			//cout << "position x coord before " << get<position>(particles)[i][0] << endl;

			//cout << "position y coord before " << get<position>(particles)[i][1] << endl;
 					
			//get<position>(particles)[i] += vdouble2(change_x,change_y);

			get<position>(particles)[i] += vdouble2(change_x,change_y);

			//cout << "position x coord " << get<position>(particles)[i][0] << endl;

			//cout << "position y coord " << get<position>(particles)[i][1] << endl;
	
		}//end of < particle_size
				
}//end of <N_steps


	// save particles after they move
	vtkWriteGrid("fixed",0,particles.get_grid(true)); 

}
