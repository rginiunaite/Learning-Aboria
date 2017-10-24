#include "Aboria.h"
#include <random>
#include <map>
#include <Eigen/Core>
#include <algorithm> // std::unique_copy

#include <math.h>

// just to check if changes work

using namespace std;
using namespace Aboria;
using namespace Eigen; // objects VectorXf, MatrixXf

int main() {


	// Gaussian gradient

	std::default_random_engine generator;
        std::normal_distribution<double> dist(0.0,1.0); // for x and y variables same mean and variance

	const int length = 100000; // length of the chemoattractant vector

	VectorXf chemo_x(length), chemo_y(length);	

	// generate Gaussian gradient for x and y variables
	for (int i =0;i<length;i++){

		chemo_x(i) = dist(generator);
		chemo_y(i) = dist(generator);
			
	}

	// find minimum value, so that we would we would not have negative entries of the histogram matrix which we will in afterwards

	double min_x = chemo_x.minCoeff();
	double min_y = chemo_y.minCoeff();

	// subtract it from each of the elements (since it is negative )
	for (int i =0;i<length;i++){

		chemo_x(i) -= min_x;
		chemo_y(i) -= min_y;			
	}

	// find max elements
	double max_x = chemo_x.maxCoeff();
	double max_y = chemo_y.maxCoeff();

	// find biggerr maximum coefficient, because we want to have a square histogram

	double bigger_max;

	if (max_x<max_y){
		bigger_max = max_y;
	}
	else{
		bigger_max = max_x;
	}

	// size of the histogram

	int lattice_size = (int)(bigger_max*10); // check if it is different, if I use round
	std::cout << "size of histogram: " << lattice_size << std::endl;
	
	// matrix for histogram of specified size
	
	MatrixXf histogram(lattice_size+1,lattice_size+1); // +1 to make sure that the size is big enough
	std::cout << "length " << length << std::endl;
   
	// create a histogram for the Gaussian x,y
	for (int i =0;i<length;i++){		
    		histogram ((int)(chemo_x[i]*10),(int)(chemo_y[i]*10)) ++;
		//std::cout << "i " << i << std::endl;
		//std::cout << "j " << j << std::endl;
		//std::cout << "x coordinate times 10 " << (int)(chemo_x[i]*10)+1 << std::endl;
		//std::cout << "y coordinate times 10 " << (int)(chemo_y[j]*10)+1 << std::endl;
    	}	


	std::cout << "Matrix created using eigen " << std::endl; 
	
	for (int i =0;i<lattice_size;i++){
		for(int j = 0;j<lattice_size;j++){
        	std::cout << histogram(i,j) << " "; 
		}
		std::cout << "\n" << std::endl; 
    	}	
	

	    /*
	     * 2D domain with a few randomly placed particles
	     */

	// step size;

	double step = 10;

	const size_t N = 100;
	//ABORIA_VARIABLE(velocity,vdouble2,"velocity")
	typedef Particles<std::tuple<>, 2> particle_type;
	//typedef Particles<std::tuple<>,2,std::vector,bucket_search_serial> particle_type;
	typedef particle_type::position position;
	particle_type particles(N);
	std::default_random_engine gen;
	std::uniform_real_distribution<double> uniform(0,lattice_size); // uniformly distributed around the histogram
	for (int i=0; i<N; ++i) {
	    get<position>(particles)[i] = vdouble2(uniform(gen),uniform(gen));
	}

	cout << "step " << step << endl;
	// Update positions based on the gradient

	for (auto& x: get<position>(particles)) {


		// choose internal positions
		if (round(x)[0]>1 && round(x)[0] < lattice_size && round(x)[1]>1 && round(x)[1] < lattice_size ){
		
			// check which of the four directions we obtain the biggest gradient

		VectorXf directions(4);

		double up = histogram(round(x)[0],round(x)[1]+1)- histogram(round(x)[0],round(x)[1]);		

cout << "up " << up << endl;
		double down = histogram(round(x)[0],round(x)[1]-1)- histogram(round(x)[0],round(x)[1]);
cout << "down " << down << endl;
		double right = histogram(round(x)[0]+1,round(x)[1])- histogram(round(x)[0],round(x)[1]);
cout << "right " << right << endl;

		double left = histogram(round(x)[0]-1,round(x)[1])- histogram(round(x)[0],round(x)[1]);

cout << "left " << left << endl;

		directions(0) = up;
		directions(1) = down;
		directions(2) = right;
		directions(3) = left;

	// save particles before they move
	vtkWriteGrid("before_fixed",0,particles.get_grid(true));
		
		// find maximum direction, if maximum is not unique, it chooses the one with the lowest coefficient, I will have to change that then it chooses randomly
		Eigen::VectorXf::Index max_index;
		double max_dir_index = directions.maxCoeff(&max_index);// do not forget from 0 to 3
		std::cout << "max index " << max_index << std::endl;
		double change;

		if (max_index == 0)// up
		{
			x += vdouble2(0,step);
		}
		if (max_index == 1)// down
		{
			x += vdouble2(0,-step);
		}
		if (max_index == 2)// right
		{
			x += vdouble2(step,0);
		}
		if (max_index == 3)// left
		{
			x += vdouble2(-step,0);
		}
	    }
	}
	particles.update_positions();
	 


	// save particles after they movev
	vtkWriteGrid("fixed",0,particles.get_grid(true));

}
