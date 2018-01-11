// cells move towords high concentration of fixed gradient

#include "Aboria.h"
#include <random>
#include <map>
#include <Eigen/Core>
#include <algorithm> // std::unique_copy
#include <iostream>// writing on a text file
#include <fstream>
#include <math.h>
#include <assert.h>

using namespace std;
using namespace Aboria;
using namespace Eigen; // objects VectorXf, MatrixXf

int main() {


	// model parameters

	int length_x = 24;//240;//40;//4; // length of the chemoattractant vector, for fixed domain
	int length_x_change = 24;//240;
	//int new_length_x_change = length_x_change;
	double initial_domain_length = 24;
	double domain_length = 24; //this variable is for the actual domain length

	double old_length = 24;
	const int length_y = 12;//120;//20;//4;
	const double diameter = (2*7.5)/10;//2 // diameter in which there have to be no cells, equivalent to size of the cell
	double cell_radius = (7.5)/10;//0.5; // cell size relative to mesh
	int N_steps = 100; // number of times the cells move up the gradient
	const size_t N = 4; // number of cells
	double l_filo = 27.5/10;//2; // sensing radius
	double diff_conc = 0.15; // how much concentration has to be bigger, so that the cell moves
	int freq_growth = 1; // domain grows linear every freq_growth time step

	// domain growth parameters

	/*double L_0 = 20;//300;
	double a = 0.08;
	//double t_s = -16;
 	int L_inf = 110;//1100;//100;//16;//870;*/
	double L_0 = initial_domain_length;//300;
	double a = 0.08/10;
	//double t_s = -16;
 	int L_inf = 50;//1100;//100;//16;//870;


	// parameters for the dynamics of chemoattractant concentration


	double D = 1/10; // to 10^5 \nu m^2/h diffusion coefficient
	double t = 0; // initialise time, redundant
	double dt = 0.1/10; // time step	
	int dx = 1; // space step in x direction
	int dy = 1; // space step in y direction
	double kai = 0.0001/10; // to 1 /h production rate of chemoattractant
	

	// parameters for internalisation

	double R = 7.5/10; // \nu m cell radius
	int lam = 100/10;//(100)/10; // to 1000 /h chemoattractant internalisation


	MatrixXf chemo(length_x, length_y), chemo_new(length_x,length_y);	

	// initialise internalisation matrix
	MatrixXf intern(length_x,length_y);

	// generate gradient which linearly increases
	for (int i = 0;i<length_x;i++){
		for (int j = 0; j< length_y; j++){
			chemo(i,j) = 1; // uniform concentration initially
			chemo_new(i,j) = 1; // this is for later updates, not sure if I will need it
			//chemo_change_len(i,j) = 1;
		}			
	}

	// domain hasn't grown there yet
	/*for(int i = length_x;i<L_inf;i++){
		for (int j = 0; j< length_y; j++){
			chemo_change_len(i,j) = 0;
		}
	}*/


	// three columns for x, y, z
	
	// form a matrix which would store x,y,z

	 MatrixXf chemo_3col(length_x*length_y,4), chemo_3col_ind(length_x*length_y,2); // need for because that is how paraview accepts data, third dimension is just zeros

	// x, y coord, 1st and 2nd columns respectively
	int k = 0;
		// it has to be 3D for paraview
	while (k<length_x*length_y){
		for (int i = 0;i<length_x;i++){
			for (int j = 0; j< length_y; j++){
				chemo_3col_ind(k,0) = i; 
				chemo_3col_ind(k,1) = j;
				chemo_3col(k,2) = 0;
				k += 1;
			}			
		}
	}


	// y and x (initially) column  
	for (int i=0;i<length_x*length_y;i++){
		chemo_3col(i,1) = chemo_3col_ind(i,1);
		chemo_3col(i,0) = chemo_3col_ind(i,0);
	}	
	
	
	// u column
	for (int i=0;i<length_x*length_y;i++){
		chemo_3col(i,3) = chemo(chemo_3col_ind(i,0),chemo_3col_ind(i,1));
	}

	// save data to plot chemoattractant concentration in MATLAB
	ofstream output("matrix_3col.csv");
	
		output << "x, y, z, u" << "\n" << endl;


	for (int i=0;i<length_x*length_y;i++){
		for(int j=0;j<4;j++){
			output << chemo_3col(i,j) << ", ";
		}
		output << "\n" << endl;
	}
	

	    /*
	     * 2D domain with a few randomly placed particles
	     */

	/*
	 * initial cells of fixed radius 	
	 */

	//ABORIA_VARIABLE(velocity,vdouble2,"velocity")
 	ABORIA_VARIABLE(radius,double,"radius")
	typedef Particles<std::tuple<radius>, 2> particle_type;
	//typedef Particles<std::tuple<>,2,std::vector,bucket_search_serial> particle_type;
	typedef particle_type::position position;
	particle_type particles;
	std::default_random_engine gen;
	std::uniform_real_distribution<double> uniform(2,length_y-1);

	/*
         * initialise neighbour search with 2d cuboid domain,
         * periodic in x and y
         */

        particles.init_neighbour_search(vdouble2(0,0), vdouble2(length_x,length_y), vbool2(false,false));	
  
	for (int i=0; i<N; ++i) {
		bool free_position = false;
		particle_type::value_type p;
		get<radius>(p) = cell_radius;
		while(free_position == false){
	    		get<position>(p) = vdouble2(cell_radius,uniform(gen)); // x=2, uniformly in y
			free_position = true;
		        /*
		         * loop over all neighbouring particles within "diameter=2*radius" distance
		         */
		        for (auto tpl: euclidean_search(particles.get_query(),get<position>(p),diameter)) {
		            /*
		             * tpl variable is a tuple containing:
		             *  (0) -> neighbouring particle value_type
		             *  (1) -> relative position of neighbouring particle
		             *         from query point
		             */
		            const vdouble2& dx = std::get<1>(tpl);
		            const particle_type::value_type& j = std::get<0>(tpl);
		            if (dx.norm() <  diameter) {
		                free_position = false;
		                break;
		            }
			}
		}
		particles.push_back(p);
	}
		
	// save particles before they move
	vtkWriteGrid("before",0,particles.get_grid(true));


	// choose a set of random number between 0 and pi, to avoid more rejections when it goes backwords (it would always be rejected)
	std::default_random_engine gen1;
	std::uniform_real_distribution<double> uniformpi(0,2*M_PI); // can only move forward
	VectorXf random_angle(N_steps*particles.size()*particles.size()*N_steps);  



	// make sure I have enough because more particles will appear
	for (int i = 0; i<N_steps*particles.size()*particles.size()*N_steps;i++){
		random_angle(i) = uniformpi(gen1);		
		//cout << "angle to move " << random_angle(i) << endl;
	}

	int rand_num_count = 0;

	
	int domain_len_der = 0; // for now assume linear growth
	double update = 0;



		// choose diff_domain random numbers from 2 to length_x, to account for positions where extra entry will be added
	//cout << "stops, this is diff_domain " << diff_domain << endl;
	
		/* nesamone
		std::default_random_engine gen2;
		std::uniform_real_distribution<double> uniformlen(2,length_x-1); // random position in x
		VectorXf random_num_long(1000);  //vector to store lost of random numbers random 
		cout << "unhappy here" << endl;
		for (int i = 0; i<1000;i++){
			random_num_long(i) = int(uniformlen(gen2));
		cout << "unhappy here, LT" << endl;		

		}

		int rand_count =0; // this will loop through long random vector
	*/

	for (int t = 0; t < N_steps; t++){
	// insert new cells at the start of the domain at insertion time (have to think about this insertion time)

		 if (t % 5 == 0 ){
			bool free_position = false;
			particle_type::value_type p;
			get<radius>(p) = cell_radius;
			
		    		get<position>(p) = vdouble2(cell_radius,uniform(gen)); // x=2, uniformly in y
				free_position = true;
				/*
				 * loop over all neighbouring particles within "dem_diameter" distance
				 */
				for (auto tpl: euclidean_search(particles.get_query(),get<position>(p),diameter)) {
				    /*
				     * tpl variable is a tuple containing:
				     *  (0) -> neighbouring particle value_type
				     *  (1) -> relative position of neighbouring particle
				     *         from query point
				     */
				    const vdouble2& dx = std::get<1>(tpl);
				    const particle_type::value_type& j = std::get<0>(tpl);
				    if (dx.norm() <  diameter) {
				        free_position = false;
				        break;
				    }
				}
			if (free_position == true){
			particles.push_back(p);}
		} 


	/////////////////////////////////////	
		// grow domain 

		
		domain_len_der = 0;

		if (t % freq_growth == 0){
			//cout << "are you never in " << endl;

			//new_length_x = int((L_inf*exp(a*t))/ (L_inf/L_0 + exp(a*t) - 1) );

			//domain_length = domain_length + 1.0;
			domain_length = ((L_inf*exp(a*t))/ (L_inf/L_0 + exp(a*t) - 1) );
			//length_x_change = int( length_x_change+1); // change in the domain length

			
			
			//update = (double(diff_domain)/double(domain_length));
			//cout << "update "<< update << endl;
			//cout << "diff domain " << diff_domain << endl;
			domain_len_der = ((a*L_inf*exp(a*t))/ (L_inf/L_0 + exp(a*t) - 1) - (a*L_inf*exp(2*a*t))/(L_inf/L_0 + exp(a*t) - 1) );
			cout << "derivative " << domain_len_der;
			//domain_len_der = 1;
		}
			//cout << "diff domain outside " << diff_domain << endl;

		// update chemoattractant profile 




		// internalisation
		for (int i =0;i<length_x;i++){
			for(int j = 0;j<length_y;j++){
				//go thorugh all the cells
				for (int k =0; k<particles.size();k++){
					vdouble2 x;
					x = get<position>(particles[k]);
					//intern(i,j) = intern(i,j) + exp(-((i-x[0]-update)*(i-x[0]-update)+(j-x[1])*(j-x[1]))/(2*R*R));
					intern(i,j) = intern(i,j) + exp(- (((domain_length/length_x)*i-x[0])*((domain_length/length_x)*i-x[0])+(j-x[1])*(j-x[1]))/(2*R*R));
					//intern(i,j) = intern(i,j) + exp(- ((domain_length)*(domain_length)*(i-x[0])*(i-x[0])+(j-x[1])*(j-x[1]))/(2*R*R));
				}			
			}
    		}	

		for (int i=1;i<length_x-1;i++){
				for (int j=1;j<length_y-1;j++){
					chemo_new(i,j) = dt * (D*((1/((domain_length/length_x)*(domain_length/length_x))) * (chemo(i+1,j)-2*chemo(i,j)+chemo(i-1,j))/(dx*dx) + (chemo(i,j+1)- 2* chemo(i,j)+chemo(i,j-1))/(dy*dy)) - (chemo(i,j)*lam / (2*M_PI*R*R)) * intern(i,j) + kai*chemo(i,j)*(1-chemo(i,j)) - double(domain_len_der)/double(domain_length) * chemo(i,j) ) + chemo(i,j);
					//chemo_new(i,j) = dt * (D*((1/((domain_length)*(domain_length)))* (chemo(i+1,j)-2*chemo(i,j)+chemo(i-1,j))/(dx*dx) + (chemo(i,j+1)- 2* chemo(i,j)+chemo(i,j-1))/(dy*dy)) - (chemo(i,j)*lam / (2*M_PI*R*R)) * intern(i,j) + kai*chemo(i,j)*(1-chemo(i,j)) - double(domain_len_der)/double(domain_length) *chemo(i,j) ) + chemo(i,j);
				}
				//cout << "print the internalisation term " << intern(i,j) << endl;
				//cout << "new chemo " << chemo_new(i,j) << endl;
			}
	
	// zero flux boundary conditions 


		for (int i=0;i<length_y;i++){
			chemo_new(0,i) = chemo_new(1,i);
			chemo_new(length_x-1,i) = chemo_new(length_x-2,i);
		
		}

		for (int i=0;i<length_x;i++){
			chemo_new(i,0) = chemo_new(i,1);
			chemo_new(i,length_y-1) = chemo_new(i,length_y-2);
		
		}



		chemo = chemo_new;
		

	
		/*
			Rescale x coordinates properly
		*/


			// three columns for x, y, z
	
			// form a matrix which would store x,y,z, y -same as before


			cout << "old length " << old_length << endl;
			cout << "domain length " << domain_length << endl;
			// x column
			for (int i=0;i<length_x*length_y;i++){
				chemo_3col(i,0) = chemo_3col_ind(i,0)*(domain_length/length_x);
			}
			//cout << "domain length ratio " << domain_length/length_x << endl; 
		
			// u column
			for (int i=0;i<length_x*length_y;i++){
				chemo_3col(i,3) = chemo(chemo_3col_ind(i,0),chemo_3col_ind(i,1));
			}
			
			// save data to plot chemoattractant concentration 
			ofstream output("matrix_growing_domain" + to_string(t) +".csv");
	
				output << "x, y, z, u" << "\n" << endl;


			for (int i=0;i<length_x*length_y;i++){
				for(int j=0;j<4;j++){
					output << chemo_3col(i,j) << ", ";
				}
				output << "\n" << endl;
			}
	

	/// update positions uniformly based on the domain growth

		if (t % freq_growth == 0){
			
			for (int i = 0; i< particles.size();i++){
			get<position>(particles)[i] *= vdouble2((domain_length/old_length), 1);
			}
			old_length = domain_length;
		}


	/////////////////////////////////////


		// Update positions based on the gradient




	// SEARCH MULTIPLE TIMES 

	
	//for (int j=0; j<N_steps;j++){
		for (int i=0; i < particles.size(); i++){
		

			vdouble2 x;
			x = get<position>(particles[i]);
			//cout << "particles.size " << i << endl;
			//cout << "print id " << get<id>(particles[i]) << endl;
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



			// HAVE TO CHOOSE TWO ANGLES 
			int sign_x_2, sign_y_2;
			if(sin(random_angle(rand_num_count+1))<0){
				sign_x_2=-1;
			}else{sign_x_2=1;}

			if(cos(random_angle(rand_num_count+1))<0){
				sign_y_2=-1;
			}else{sign_y_2=1;}

			
			
			
			
			//cout << "sin of the angle " << sin(random_angle(rand_num_count)) << endl;
			//cout << "cos of the angle " << cos(random_angle(rand_num_count)) << endl;

		
			// check if the x coordinates are not out of the domain, if they are, ignore that step
		
			//cout << "chemo coord x before the test " << round(x[0]+sin(random_angle(rand_num_count))+sign_x*cell_radius) << endl;
			//cout << "chemo coord y before the test " << round(x[0]+cos(random_angle(rand_num_count))+sign_y*cell_radius) << endl;

		
	//cout << "domain length " << domain_length << endl;
	//cout << "size chemo " << chemo.rows() << "x" << chemo.cols() << endl;
		// make sure that the next position is an entry of the matrix
		// have to multiply by length_x/domain_length to come back to the same matrix size
			if (round((x[0] * (length_x/domain_length)+sin(random_angle(rand_num_count))+sign_x*l_filo) )>-1 && round((x[0] * (length_x/domain_length)+sin(random_angle(rand_num_count))+sign_x*l_filo)) < length_x && round(x[1]+ cos(random_angle(rand_num_count))+sign_y*l_filo) >-1 && round(x[1]+ cos(random_angle(rand_num_count))+sign_y*l_filo)<length_y ){

			//cout << "sin of the angle (inside the loop) " << sin(random_angle(rand_num_count)) << endl;
			//cout << "cos of the angle (inside the loop) " << cos(random_angle(rand_num_count)) << endl;

			// check if the gradient in the other position is larger, if yes, move to that position, x changes by sin and y to cos, because of the the chemo is defined. 

 			
			//cout << "problem here with the coord if before 57" << endl;

			//cout << "x coord " << round(x[0]) << endl;
			//cout << "x up " << round(x[0]+sin(random_angle(rand_num_count))+sign_x*l_filo) << endl;
			//cout << "y coord " << round(x)[1] << endl;
			//cout << "y up " << round(x[1]+ cos(random_angle(rand_num_count))+sign_y*l_filo) <<endl;




			//cout << "chemo conc in current site " << chemo(round(x)[0],round(x)[1])<< endl;	
			//cout << "chemo conc in other site " << chemo(round(x[0]* (length_x/domain_length) + sin(random_angle(rand_num_count)) +sign_x*l_filo),round(x[1]+cos(random_angle(rand_num_count))+sign_y*l_filo))<< endl;


		// need that + diff_conc to make sure that the concentration is sufficiently bigger
		// map to the fixed domain to check the difference in concentration


			if (chemo(round((x)[0]* (length_x/domain_length)),round(x)[1])+diff_conc < chemo(round((x[0]* (length_x/domain_length)+sin(random_angle(rand_num_count))+sign_x*l_filo)),round(x[1] + cos(random_angle(rand_num_count))+sign_y*l_filo))){


			//cout << "compare, first one " << chemo_change_len(round(x)[0],round(x)[1]) << endl;
			//cout << "second one, should be bigger " << chemo_change_len(round(x[0]+sin(random_angle(rand_num_count))+sign_x*l_filo),round(x[1]+ cos(random_angle(rand_num_count))+sign_y*l_filo)) << endl;
			//cout << "can enter" << endl;
					//cout << "x coord " << x[0] << endl;
					//cout << "x up " << sin(random_angle(rand_num_count))+sign_x*cell_radius << endl;
					//cout << "y coord " << x[1] << endl;
					//cout << "y up " << cos(random_angle(rand_num_count))+sign_y*cell_radius << endl;	
					//cout << "chemo coonc in current site " << chemo(round(x)[0],round(x)[1])<< endl;	
					//cout << "chemo coonc in other site " << chemo(round(x[0] +sin(random_angle(rand_num_count)) +sign_x*cell_radius),round(x[1]+cos(random_angle(rand_num_count))+sign_y*cell_radius))<< endl;	
					//get<position>(particles)[i] += vdouble2(sin(random_angle(rand_num_count))+sign_x*cell_radius, cos(random_angle(rand_num_count))+sign_y*cell_radius);

		// check if there are no cells around the position where I want to move

		        /*
		         * loop over all neighbouring particles within a euclidean distance
		         * of size "diameter" after it moved
		         */
			x += vdouble2(sin(random_angle(rand_num_count)), cos(random_angle(rand_num_count)));
			//cout << "print id " << id_[x] << endl;
		

			//cout << "Position "<< x << endl;
			int count_position = 0;
			bool free_position = true; // check if the neighbouring position is free

			// if this loop is entered, it means that there is another cell where I want to move 
				for (const auto& k: euclidean_search(particles.get_query(),x,diameter)) {

				    count_position += 1; // just to check if this works
					particle_type::const_reference b = std::get<0>(k);
		    			const vdouble2& dx = std::get<1>(k);
		   			//cout << "Found a particle with dx = " << dx << " and id = " << get<id>(b) << "\n";	
	
					cout << "id of b " << get<id>(b) << endl;
					//for (int i=0; i < particles.size(); i++) {
	   					if (get<id>(b) != get<id>(particles[i])){ // check if it is not the same particle
						//cout << "reject step " << 1 << endl;
				    		free_position = false;}
						//}
				
				    //break;
				}

		
			//cout << "print position " << count_position << endl;
				if (free_position == true){
					get<position>(particles)[i] += vdouble2(sin(random_angle(rand_num_count)), cos(random_angle(rand_num_count))); // update if nothing is in the next position
				}

				}// update the position
					//else{cout << "rejected" << endl;}
			} //check if not outside the domain
		
				rand_num_count += 2; // update random number count	
				
		}// go through all the particles



		//particles.update_positions();
			// save particles after they move
				    /*
			     * on every i/o step write particle container to a vtk
			     * unstructured grid file
			     */
			    //cout << "." << flush;
		#ifdef HAVE_VTK
		vtkWriteGrid("particles",t,particles.get_grid(true));    
		#endif

		//for (int i =0;i<5;++i){cout << "koks rezas " << i <<endl;}
		//length_x_change = new_length_x_change;
		//cout << "new length " << new_length_x_change << endl;
	}//all time steps

for (int i=0; i < particles.size(); i++) {
   cout << "Positions = " << get<position>(particles[i]) << "\n";
}


}
