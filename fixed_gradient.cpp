#include "Aboria.h"
#include <random>
#include <map>
#include <Histogram.h>
#include <Eigen/Core>
#include <algorithm>    // std::max

// just to check if changes work


using namespace Aboria;
using Eigen::MatrixXd;

int main() {


	// Gaussian gradient

	std::default_random_engine generator;
        std::normal_distribution<double> dist(0.0,1.0); // for x and y variables same mean and variance

	
	const int length = 100;

	double chemo_x[length], chemo_y[length];

	for (int i =0;i<length;i++){
		for (int j=0;j<length;j++){
			chemo_x[i] = dist(generator);
			chemo_y[i] = dist(generator);
			//std::cout << "Gaussian distribution x " << chemo_x[i] << std::endl;
			//std::cout << "Gaussian distribution y " << chemo_y[i] << std::endl;
		}
	}



	// Create a histogram of values gradient values 
    int lattice_size =25;
   
    int histogram[lattice_size][lattice_size]={};
    for (int i =0;i<length;i++){
	for(int j = 0;j<length;j++){
		if ((int)(chemo_x[i]*10) +10 < 0 || (int)(chemo_y[j]*10)+10 < 0){
			int a=1; // nothing
		}else
        	histogram [(int)(chemo_x[i]*10)+10] [(int)(chemo_y[j]*10)+10] ++;
		//std::cout << "x coordinate times 10 " << (int)(chemo_x[i]*10) << std::endl;
		//std::cout << "y coordinate times 10 " << (int)(chemo_y[j]*10) << std::endl;
	}
    }	
 
	


  for (int i = 0; i < lattice_size; i++ ) {
      for (int j = 0; j < lattice_size; j++ ) {        
            cout << histogram[i][j] << " ";
      }
      cout << endl;
    }
		
	MatrixXd m(2,2);

	m(0,0) = 3;
	m(1,1) = 5;

	int myint[] = {1,2};


	std::cout << "Matrix created using eigen " << std::endl; // nesamone
	std::cout << std::max_element (myint,myint +2) << std::endl;


	


	    /*
	     * 2D domain with a few randomly placed particles
	     */

	const size_t N = 100;
	//ABORIA_VARIABLE(velocity,vdouble2,"velocity")
	typedef Particles<std::tuple<>, 2> particle_type;
	//typedef Particles<std::tuple<>,2,std::vector,bucket_search_serial> particle_type;
	typedef particle_type::position position;
	particle_type particles(N);
	std::default_random_engine gen;
	std::uniform_real_distribution<double> uniform(-1,1);
	for (int i=0; i<N; ++i) {
	    get<position>(particles)[i] = vdouble2(uniform(gen),uniform(gen));
	}
	 


	// saves the file
	vtkWriteGrid("fixed",0,particles.get_grid(true));

}
