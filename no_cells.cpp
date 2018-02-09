/*
cells move towords high concentration of varying chemoattracted, the movement is directed,
 but if cells do not sense higher concentration, they move randomly
*/

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

double func( double lam) {


    // model parameters


    int length_x = 30;//240;//40;//4; // length of the chemoattractant vector, for fixed domain
    double domain_length = 30; //this variable is for the actual domain length
    double old_length = 30;// this is important for the update of the positions
    const int length_y = 12;//120;//20;//4;
    const double diameter = (2 * 7.5) / 10;//2 // diameter in which there have to be no cells, equivalent to size of the cell
    double cell_radius = (7.5) / 10;//0.5; // radius of a cell
    int N_steps = 200; // number of times the cells move up the gradient
    const size_t N = 4; // initial number of cells
    double l_filo = 27.5 / 10;//2; // sensing radius
    //double diff_conc = 0.5; // sensing threshold, i.e. how much concentration has to be bigger, so that the cell moves in that direction
    int freq_growth = 200; // determines how frequently domain grows (actually not relevant because it will go every timestep)
    int insertion_freq = 201;
    double speed_l = 1; // speed of a leader cell
    double speed_f = 1; // speed of a follower cell


    double domain_fraction = 0.5; // fraction of the end of the domain that we count the percentage of cells in

    /*
    double L_0 = initial_domain_length;//300;
    double a = 0.08/10;
     int L_inf = 50;//1100;//100;//16;//870;
    */

    // for logistic growth with delay

    double L_0 = 30; // will have to make this consistent with actual initial length
    double a = 0.008;
    double L_inf = 87;
    double t_s = 16;
    double constant = 30;

    double domain_len_der = 0; // for now assume linear growth

    // parameters for the dynamics of chemoattractant concentration

    double D = 1 / 1; // to 10^5 \nu m^2/h diffusion coefficient
    double t = 0; // initialise time, redundant
    double dt = 0.1 / 10; // time step
    int dx = 1; // space step in x direction
    int dy = 1; // space step in y direction
    double kai = 0.0001 / 10; // to 1 /h production rate of chemoattractant


    // parameters for internalisation

    double R = cell_radius; // \nu m cell radius
    //int lam = 100 / 10;//(100)/10; // to 1000 /h chemoattractant internalisation

    // matrix that stores the values of concentration of chemoattractant
    MatrixXf chemo(length_x, length_y), chemo_new(length_x, length_y);

    // initialise internalisation matrix
    MatrixXf intern(length_x, length_y);


    // generate initial concentration of chemoattractant
    for (int i = 0; i < length_x/2; i++) {
        for (int j = 0; j < length_y; j++) {
            chemo(i, j) = 0; // uniform concentration initially
            chemo_new(i, j) = 0; // this is for later updates
        }
    }

    for (int i = length_x/2; i < length_x; i++) {
        for (int j = 0; j < length_y; j++) {
            chemo(i, j) = i; // uniform concentration initially
            chemo_new(i, j) = i; // this is for later updates
        }
    }



    // four columns for x, y, z, u (z is necessaty for paraview)

    // form a matrix which would store x,y,z,u

    MatrixXf chemo_3col(length_x * length_y, 4), chemo_3col_ind(length_x * length_y,
                                                                2); // need for because that is how paraview accepts data, third dimension is just zeros

    // x, y coord, 1st and 2nd columns respectively
    int k = 0;
    // it has to be 3D for paraview
    while (k < length_x * length_y) {
        for (int i = 0; i < length_x; i++) {
            for (int j = 0; j < length_y; j++) {
                chemo_3col_ind(k, 0) = i;
                chemo_3col_ind(k, 1) = j;
                chemo_3col(k, 2) = 0;
                k += 1;
            }
        }
    }


    // y and x (initially) column
    for (int i = 0; i < length_x * length_y; i++) {
        chemo_3col(i, 1) = chemo_3col_ind(i, 1);
        chemo_3col(i, 0) = chemo_3col_ind(i, 0);
    }


    // u column
    for (int i = 0; i < length_x * length_y; i++) {
        chemo_3col(i, 3) = chemo(chemo_3col_ind(i, 0), chemo_3col_ind(i, 1));
    }

    // save data to plot chemoattractant concentration in MATLAB
    ofstream output("matrix_3col.csv");

    output << "x, y, z, u" << "\n" << endl;


    for (int i = 0; i < length_x * length_y; i++) {
        for (int j = 0; j < 4; j++) {
            output << chemo_3col(i, j) << ", ";
        }
        output << "\n" << endl;
    }


    /*
     * 2D domain with a few randomly placed particles
     */

    /*
     * initial cells of fixed radius
     */

    /*
         * initialise neighbour search with 2d cuboid domain,
         * periodic in x and y
         */


    for (int t = 0; t < N_steps; t++) {
        // insert new cells at the start of the domain at insertion time (have to think about this insertion time)



        /////////////////////////////////////
        // grow domain


        domain_len_der = 0;

        if (t % freq_growth == 0) {
            //cout << "are you never in " << endl;

            // no time delay
            //domain_length = ((L_inf*exp(a*t))/ (L_inf/L_0 + exp(a*t) - 1) );

            //domain_len_der = ((a*L_inf*exp(a*t))/ (L_inf/L_0 + exp(a*t) - 1) - (a*L_inf*exp(2*a*t))/(L_inf/L_0 + exp(a*t) - 1) );

            // from the paper
            //domain_length = ((L_inf*exp(a*(t-t_s)))/ (L_inf/L_0 + exp(a*(t-t_s)) - 1) ) + constant;

            //domain_len_der = ((a*L_inf*exp(a*(t-t_s)))/ (L_inf/L_0 + exp(a*(t-t_s)) - 1) - (a*L_inf*exp(2*a*(t-t_s)))/(L_inf/L_0 + exp(a*(t-t_s)) - 1) );

            // with time delay and constant to make the initial conditions consistent
            domain_length = ((L_inf * exp(a * (t - t_s))) / (L_inf / L_0 + exp(a * (t - t_s)) - 1)) + constant;

            domain_len_der = ((a * L_inf * exp(a * (t - t_s))) / (L_inf / L_0 + exp(a * (t - t_s)) - 1) -
                              (a * L_inf * exp(2 * a * (t - t_s))) / (L_inf / L_0 + exp(a * (t - t_s)) - 1));


        }
        //cout << "diff domain outside " << diff_domain << endl;

        // update chemoattractant profile





        for (int i = 1; i < length_x - 1; i++) {
            for (int j = 1; j < length_y - 1; j++) {
                chemo_new(i, j) = dt * (D * ((1 / ((domain_length / length_x) * (domain_length / length_x))) *
                                             (chemo(i + 1, j) - 2 * chemo(i, j) + chemo(i - 1, j)) / (dx * dx) +
                                             (chemo(i, j + 1) - 2 * chemo(i, j) + chemo(i, j - 1)) / (dy * dy)) +
                                        kai * chemo(i, j) * (1 - chemo(i, j)) - lam *chemo(i,j) +
                                        double(domain_len_der) / double(domain_length) * chemo(i, j)) + chemo(i, j);

            }
            //cout << "print the internalisation term " << intern(i,j) << endl;
            //cout << "new chemo " << chemo_new(i,j) << endl;
        }


        // zero flux boundary conditions


        for (int i = 0; i < length_y; i++) {
            chemo_new(0, i) = chemo_new(1, i);
            chemo_new(length_x - 1, i) = chemo_new(length_x - 2, i);

        }

        for (int i = 0; i < length_x; i++) {
            chemo_new(i, 0) = chemo_new(i, 1);
            chemo_new(i, length_y - 1) = chemo_new(i, length_y - 2);

        }


        chemo = chemo_new;



        /*
            Rescale x coordinates properly
        */


        // three columns for x, y, z

        // form a matrix which would store x,y,z, y -same as before


        //cout << "old length " << old_length << endl;
        //cout << "domain length " << domain_length << endl;
        // x column
        for (int i = 0; i < length_x * length_y; i++) {
            chemo_3col(i, 0) = chemo_3col_ind(i, 0) * (domain_length / length_x);
        }
        //cout << "domain length ratio " << domain_length/length_x << endl;

        // u column
        for (int i = 0; i < length_x * length_y; i++) {
            chemo_3col(i, 3) = chemo(chemo_3col_ind(i, 0), chemo_3col_ind(i, 1));
        }

        // save data to plot chemoattractant concentration
        ofstream output("matrix_growing_domain" + to_string(t) + ".csv");

        output << "x, y, z, u" << "\n" << endl;


        for (int i = 0; i < length_x * length_y; i++) {
            for (int j = 0; j < 4; j++) {
                output << chemo_3col(i, j) << ", ";
            }
            output << "\n" << endl;
        }




    }//all time steps



// count the number of particles that are at the last 10% of the domain


}


int main(){

    const int number_parameters = 1; // parameter range

    // define parameters that I will change
    //VectorXf slope, threshold;
    array<double, number_parameters> lam, threshold;
    //array<double,number_parameters,number_parameters>;

    MatrixXf density(number_parameters,number_parameters); // need for because that is how paraview accepts data, third dimension is just zeros



    for (int i=0; i<number_parameters; i++){

        lam[i] = 500;//*(i+1);//0.1;
        //threshold[i] = 0.5;//0.01*(i+1);// 0.01;
        threshold[i] = 1*(i+1);// 0.01;


    }


    for (int i = 0; i < number_parameters; i++){

        for (int j = 0; j < number_parameters; j++){

            density(i,j) = func(lam[i]);
            cout << "number of i " << i << endl;

        }

    }


    // save data to plot chemoattractant concentration
    ofstream output("no_cells.csv");

    MatrixXf density_3col(number_parameters * number_parameters, 4), density_3col_ind(number_parameters * number_parameters,
                                                                                      2); // need for because that is how paraview accepts data, third dimension is just zeros



    // x, y coord, 1st and 2nd columns respectively
    int k = 0;
    // it has to be 3D for paraview
    while (k < number_parameters * number_parameters) {
        for (int i = 0; i < number_parameters; i++) {
            for (int j = 0; j < number_parameters; j++) {
                density_3col_ind(k, 0) = i;
                density_3col_ind(k, 1) = j;
                density_3col(k, 2) = 0;
                k += 1;
            }
        }
    }


    // y and x (initially) column
    for (int i = 0; i < number_parameters * number_parameters; i++) {
        density_3col(i, 1) = density_3col_ind(i, 1);
        density_3col(i, 0) = density_3col_ind(i, 0);
    }


    // u column
    for (int i = 0; i < number_parameters * number_parameters; i++) {
        density_3col(i, 3) = density(density_3col_ind(i, 0), density_3col_ind(i, 1));
    }

    output << "x, y, z, u" << "\n" << endl;


    for (int i = 0; i < number_parameters * number_parameters; i++) {
        for (int j = 0; j < 4; j++) {
            output << density_3col(i, j) << ", ";
        }
        output << "\n" << endl;
    }





    // This might be useful for matlab
    ofstream output2("no_cells_matlab.csv");

    for (int i = 0; i < number_parameters; i++){

        for (int j = 0; j < number_parameters; j++){

            output2 << density(i,j) << ", ";

        }
        output2 << "\n" << endl;
    }




}
