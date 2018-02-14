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
using namespace Eigen; // objects VectorXf, MatrixXd

double func(double diff_conc, double lam, int n_seed) {


    // model parameters


    int length_x = 30;//240;//40;//4; // length of the chemoattractant vector, for fixed domain
    double domain_length = 30; //this variable is for the actual domain length
    double old_length = 30;// this is important for the update of the positions
    const int length_y = 12;//120;//20;//4;

    const double cell_radius = (7.5) / 10;//0.5; // radius of a cell
    const double diameter =  2 *cell_radius;//2 // diameter in which there have to be no cells, equivalent to size of the cell
    int N_steps = 200; // number of times the cells move up the gradient
    const size_t N = 4; // initial number of cells
    double l_filo = 27.5 / 10;//2; // sensing radius
    //double diff_conc = 0.5; // sensing threshold, i.e. how much concentration has to be bigger, so that the cell moves in that direction
    int freq_growth = 1; // determines how frequently domain grows (actually not relevant because it will go every timestep)
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
    double L_inf = 86.7;
    double t_s = 16;
    double constant = 29;


    double domain_len_der = 0; // for now assume linear growth

    // parameters for the dynamics of chemoattractant concentration

    double D = 1 / 1; // to 10^5 \nu m^2/h diffusion coefficient
    double t = 0; // initialise time, redundant
    double dt = 0.0001; // time step
    int dx = 1; // space step in x direction
    int dy = 1; // space step in y direction
    double kai = 0.0001 / 10; // to 1 /h production rate of chemoattractant


    // parameters for internalisation

    double R = cell_radius; // \nu m cell radius
    //int lam = 100 / 10;//(100)/10; // to 1000 /h chemoattractant internalisation

    // matrix that stores the values of concentration of chemoattractant
    MatrixXd chemo = MatrixXd::Zero(length_x, length_y);
    MatrixXd chemo_new = MatrixXd::Zero(length_x, length_y);

    // initialise internalisation matrix
    MatrixXd intern = MatrixXd::Zero(length_x, length_y);

    // generate initial concentration of chemoattractant
    for (int i = 0; i < length_x; i++) {
        for (int j = 0; j < length_y; j++) {
            chemo(i, j) = 1; // uniform concentration initially
            chemo_new(i, j) = 1; // this is for later updates
        }
    }

    for (int i = 0; i < length_x; i++) {
        for (int j = 0; j < length_y; j++) {
            assert(chemo(i, j) >= 0 || chemo(i, j) <= 1);
        }
    }

    // four columns for x, y, z, u (z is necessary for paraview)

    // form a matrix which would store x,y,z,u

    MatrixXd chemo_3col = MatrixXd::Zero(length_x * length_y, 4);
    MatrixXd chemo_3col_ind = MatrixXd::Zero(length_x*length_y,2); // need for because that is how paraview accepts data, third dimension is just zeros

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

//    // save data to plot chemoattractant concentration in MATLAB
//    ofstream output("matrix_3col.csv");
//
//    output << "x, y, z, u" << "\n" << endl;
//
//
//    for (int i = 0; i < length_x * length_y; i++) {
//        for (int j = 0; j < 4; j++) {
//            output << chemo_3col(i, j) << ", ";
//        }
//        output << "\n" << endl;
//    }


    /*
     * 2D domain with a few randomly placed particles
     */

    /*
     * initial cells of fixed radius
     */

    //ABORIA_VARIABLE(velocity,vdouble2,"velocity")
    ABORIA_VARIABLE(radius, double, "radius")
    typedef Particles<std::tuple<radius>, 2> particle_type;
    //typedef Particles<std::tuple<>,2,std::vector,bucket_search_serial> particle_type;
    typedef particle_type::position position;
    particle_type particles;
    std::default_random_engine gen;
    gen.seed(n_seed);
    std::uniform_real_distribution<double> uniform(2, length_y - 1);

    /*
         * initialise neighbour search with 2d cuboid domain,
         * periodic in x and y
         */

    particles.init_neighbour_search(vdouble2(0, 0), 5*vdouble2(length_x, length_y), vbool2(false, false));

    for (int i = 0; i < N; ++i) {
        bool free_position = false;
        particle_type::value_type p;
        get<radius>(p) = cell_radius;
        while (free_position == false) {
            get<position>(p) = vdouble2(cell_radius, uniform(gen)); // x=2, uniformly in y
            free_position = true;
            /*
             * loop over all neighbouring particles within "diameter=2*radius" distance
             */
            for (auto tpl: euclidean_search(particles.get_query(), get<position>(p), diameter)) {
                /*
                 * tpl variable is a tuple containing:
                 *  (0) -> neighbouring particle value_type
                 *  (1) -> relative position of neighbouring particle
                 *         from query point
                 */
                const vdouble2 &dx = std::get<1>(tpl);
                const particle_type::value_type &j = std::get<0>(tpl);
                if (dx.norm() < diameter) {
                    free_position = false;
                    break;
                }
            }
        }
        particles.push_back(p);
    }

    // save particles before they move
    vtkWriteGrid("before", 0, particles.get_grid(true));
    vtkWriteGrid("particles", t, particles.get_grid(true));


    // choose a set of random number between 0 and pi, to avoid more rejections when it goes backwords (imagine there are cells behind)
    std::default_random_engine gen1;
    std::uniform_real_distribution<double> uniformpi(0, 2 * M_PI); // can only move forward

    for (int i = 0; i < length_x; i++) {
        for (int j = 0; j < length_y; j++) {
            assert(chemo(i, j) >= 0 || chemo(i, j) <= 1);
        }
    }

    for (int t = 0; t < N_steps; t++) {
        // insert new cells at the start of the domain at insertion time (have to think about this insertion time)

        if (t % insertion_freq == 0) {
            bool free_position = false;
            particle_type::value_type p;
            get<radius>(p) = cell_radius;

            get<position>(p) = vdouble2(cell_radius, uniform(gen)); // x=2, uniformly in y
            free_position = true;
            /*
             * loop over all neighbouring particles within "dem_diameter" distance
             */
            //particle_type::value_type closest_neighbour;
            for (auto tpl: euclidean_search(particles.get_query(), get<position>(p), diameter)) {
                /*
                 * tpl variable is a tuple containing:
                 *  (0) -> neighbouring particle value_type
                 *  (1) -> relative position of neighbouring particle
                 *         from query point
                 */
                const vdouble2 &dx = std::get<1>(tpl);
                const particle_type::value_type &j = std::get<0>(tpl);
                //cout << "position from j " << get<id>(j) << endl;
                //closest_neighbour = j;
                if (dx.norm() < diameter) {
                    free_position = false;
                    break;
                }

            }
            //cout << "position from j " << get<id>(closest_neighbour) << endl;
            if (free_position == true) {
                particles.push_back(p);
            }
        }


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



        // update chemoattractant profile




        // internalisation
        for (int i = 0; i < length_x; i++) {
            for (int j = 0; j < length_y; j++) {
                //go thorugh all the cells
                for (int k = 0; k < particles.size(); k++) {
                    vdouble2 x;
                    x = get<position>(particles[k]);

                    intern(i, j) = intern(i, j) + exp(-(((domain_length / length_x) * i - x[0]) *
                                                        ((domain_length / length_x) * i - x[0]) +
                                                        (j - x[1]) * (j - x[1])) / (2 * R * R));

                }
            }
        }

        for (int i = 1; i < length_x - 1; i++) {
            for (int j = 1; j < length_y - 1; j++) {
                chemo_new(i, j) = dt * (D * ((1 / ((domain_length / length_x) * (domain_length / length_x))) *
                                             (chemo(i + 1, j) - 2 * chemo(i, j) + chemo(i - 1, j)) / (dx * dx) +
                                             (chemo(i, j + 1) - 2 * chemo(i, j) + chemo(i, j - 1)) / (dy * dy)) -
                                        (chemo(i, j) * lam / (2 * M_PI * R * R)) * intern(i, j) +
                                        kai * chemo(i, j) * (1 - chemo(i, j)) -
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

        for (int i = 0; i < length_x; i++) {
            for (int j = 0; j < length_y; j++) {
                assert(chemo(i, j) >= 0 || chemo(i, j) <= 1);
            }
        }

        // Save matrix

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



        /*
            Rescale x coordinates properly
        */


        // three columns for x, y, z

        // form a matrix which would store x,y,z, y -same as before



        /// update positions uniformly based on the domain growth

        if (t % freq_growth == 0) {
            for (int i = 0; i < particles.size(); i++) {
                get<position>(particles)[i] *= vdouble2((domain_length / old_length), 1);
            }
            old_length = domain_length;
        }


        /////////////////////////////////////

        // Update positions based on the gradient

        // SEARCH MULTIPLE TIMES


        //for (int j=0; j<N_steps;j++){
        for (int i = 0; i < particles.size(); i++) {


            vdouble2 x;
            x = get<position>(particles[i]);
            //cout << "particles.size " << i << endl;
            //cout << "print id " << get<id>(particles[i]) << endl;
            //cout << "x coord " << x[0] << endl;


            // create an array to store random directions
            std::array<double, 3> random_angle;
//            std::array<int, 3> sign_x;
//            std::array<int, 3> sign_y;
            for (int j = 0; j < 3; j++) {

                double random_angle_tem = uniformpi(gen1);
//                int sign_x_tem, sign_y_tem;

                while (round((x[0] * (length_x / domain_length) + sin(random_angle_tem) *l_filo)) < 0 ||
                       round((x[0] * (length_x / domain_length) + sin(random_angle_tem)  * l_filo)) >
                       length_x - 1 || round(x[1] + cos(random_angle_tem)  * l_filo) < 0 ||
                       round(x[1] + cos(random_angle_tem)  * l_filo) > length_y - 1) {
                    random_angle_tem = uniformpi(gen1);

//                    if (sin(random_angle_tem) < 0) {
//                        sign_x_tem = -1;
//                    } else { sign_x_tem = 1; }
//
//                    if (cos(random_angle_tem) < 0) {
//                        sign_y_tem = -1;
//                    } else { sign_y_tem = 1; }

                }

                random_angle[j] = random_angle_tem;

//                sign_x[j] = sign_x_tem;
//                sign_y[j] = sign_y_tem;


            }


            // choose which direction to move

            // store variables for concentration at new locations
            double old_chemo = chemo((round((x)[0] * (length_x / domain_length))), round(x)[1]);
            //double new_chemo_1 = chemo(round((x[0] * (length_x / domain_length) + sin(random_angle[0]) + sign_x[0] * l_filo)),
            //round(x[1] + cos(random_angle[0]) + sign_y[0] * l_filo));
            double new_chemo_1 = chemo(round((x[0] * (length_x / domain_length) + sin(random_angle[0]) * l_filo)),
                                       round(x[1] + cos(random_angle[0])* l_filo));

            //double new_chemo_2 = chemo(round((x[0] * (length_x / domain_length) + sin(random_angle[1]) + sign_x[1] * l_filo)),
            //round(x[1] + cos(random_angle[1]) + sign_y[1] * l_filo));
            double new_chemo_2 = chemo(round((x[0] * (length_x / domain_length) + sin(random_angle[1])* l_filo)),
                                       round(x[1] + cos(random_angle[1])* l_filo));


            //if both smaller, move random direction
            //absolute
            //if (new_chemo_1 - old_chemo < diff_conc && new_chemo_2 - old_chemo < diff_conc) {


            // relative
            if ((new_chemo_1 - old_chemo)/sqrt(old_chemo) < diff_conc && (new_chemo_2- old_chemo)/sqrt(old_chemo) < diff_conc){

                x += vdouble2(sin(random_angle[2]), cos(random_angle[2]));
                //cout << "print id " << id_[x] << endl;


                //cout << "Position "<< x << endl;
                int count_position = 0;
                bool free_position = true; // check if the neighbouring position is free

                // if this loop is entered, it means that there is another cell where I want to move
                for (const auto &k: euclidean_search(particles.get_query(), x, diameter)) {

                    count_position += 1; // just to check if this works
                    particle_type::const_reference b = std::get<0>(k);
                    const vdouble2 &dx = std::get<1>(k);
                    //cout << "Found a particle with dx = " << dx << " and id = " << get<id>(b) << "\n";

                    //cout << "id of b " << get<id>(b) << endl;
                    //for (int i=0; i < particles.size(); i++) {
                    if (get<id>(b) != get<id>(particles[i])) { // check if it is not the same particle
                        //cout << "reject step " << 1 << endl;
                        free_position = false;
                    }
                }


                //cout << "print position " << count_position << endl;

                // check that the position they want to move to is free and not out of bounds
                if (free_position == true && round((x[0] * (length_x / domain_length))) > 0 &&
                    round((x[0] * (length_x / domain_length))) < length_x - 1 && round(x[1]) > 0 &&
                    round(x[1]) < length_y - 1) {
                    get<position>(particles)[i] += vdouble2(sin(random_angle[2]),
                                                            cos(random_angle[2])); // update if nothing is in the next position
                }

            }


                //cout << "stops here " << endl;
                // if first direction greater, second smaller
                //absolute
            //else if (new_chemo_1 - old_chemo > diff_conc && new_chemo_2 - old_chemo < diff_conc){

                //relative
            else if ((new_chemo_1 - old_chemo)/sqrt(old_chemo) > diff_conc && (new_chemo_2 - old_chemo)/sqrt(old_chemo) < diff_conc){

                x += vdouble2(sin(random_angle[0]), cos(random_angle[0]));
                //cout << "print id " << id_[x] << endl;


                //cout << "Position "<< x << endl;
                int count_position = 0;
                bool free_position = true; // check if the neighbouring position is free

                // if this loop is entered, it means that there is another cell where I want to move
                for (const auto &k: euclidean_search(particles.get_query(), x, diameter)) {

                    count_position += 1; // just to check if this works
                    particle_type::const_reference b = std::get<0>(k);
                    const vdouble2 &dx = std::get<1>(k);
                    //cout << "Found a particle with dx = " << dx << " and id = " << get<id>(b) << "\n";

                    //cout << "id of b " << get<id>(b) << endl;
                    //for (int i=0; i < particles.size(); i++) {
                    if (get<id>(b) != get<id>(particles[i])) { // check if it is not the same particle
                        //cout << "reject step " << 1 << endl;
                        free_position = false;
                    }
                    //}

                    //break;
                }


                //cout << "print position " << count_position << endl;
                // check that the position they want to move to is free and not out of bounds
                if (free_position == true && round((x[0] * (length_x / domain_length))) > 0 &&
                    round((x[0] * (length_x / domain_length))) < length_x - 1 && round(x[1]) > 0 &&
                    round(x[1]) < length_y - 1) {
                    get<position>(particles)[i] += vdouble2(sin(random_angle[0]),
                                                            cos(random_angle[0])); // update if nothing is in the next position
                }

            }
                // if first smaller, second bigger

            //absolute
            //else if (new_chemo_1 - old_chemo < diff_conc && new_chemo_2 - old_chemo > diff_conc){

            //relative
            else if ((new_chemo_1 - old_chemo)/sqrt(old_chemo) < diff_conc && (new_chemo_2 - old_chemo)/sqrt(old_chemo) > diff_conc){



                x += vdouble2(sin(random_angle[1]), cos(random_angle[1]));
                //cout << "print id " << id_[x] << endl;


                //cout << "Position "<< x << endl;
                int count_position = 0;
                bool free_position = true; // check if the neighbouring position is free

                // if this loop is entered, it means that there is another cell where I want to move
                for (const auto &k: euclidean_search(particles.get_query(), x, diameter)) {

                    count_position += 1; // just to check if this works
                    particle_type::const_reference b = std::get<0>(k);
                    const vdouble2 &dx = std::get<1>(k);
                    //cout << "Found a particle with dx = " << dx << " and id = " << get<id>(b) << "\n";

                    //cout << "id of b " << get<id>(b) << endl;
                    //for (int i=0; i < particles.size(); i++) {
                    if (get<id>(b) != get<id>(particles[i])) { // check if it is not the same particle
                        //cout << "reject step " << 1 << endl;
                        free_position = false;
                    }
                    //}

                    //break;
                }


                //cout << "print position " << count_position << endl;
                // check that the position they want to move to is free and not out of bounds
                if (free_position == true && round((x[0] * (length_x / domain_length))) > 0 &&
                    round((x[0] * (length_x / domain_length))) < length_x - 1 && round(x[1]) > 0 &&
                    round(x[1]) < length_y - 1) {
                    get<position>(particles)[i] += vdouble2(sin(random_angle[1]),
                                                            cos(random_angle[1])); // update if nothing is in the next position
                }
                break;
            }
                // if both greater choose the bigger one

            // absolute
            //else if (new_chemo_1 - old_chemo > diff_conc && new_chemo_2 - old_chemo > diff_conc){

                //relative
            else if ((new_chemo_1 - old_chemo)/sqrt(old_chemo) > diff_conc && (new_chemo_2 - old_chemo)/sqrt(old_chemo) > diff_conc){


                // if first is greater than the second
                if (new_chemo_1 > new_chemo_2) {
                    x += vdouble2(sin(random_angle[0]), cos(random_angle[0]));
                    //cout << "print id " << id_[x] << endl;


                    //cout << "Position "<< x << endl;
                    int count_position = 0;
                    bool free_position = true; // check if the neighbouring position is free

                    // if this loop is entered, it means that there is another cell where I want to move
                    for (const auto &k: euclidean_search(particles.get_query(), x, diameter)) {

                        count_position += 1; // just to check if this works
                        particle_type::const_reference b = std::get<0>(k);
                        const vdouble2 &dx = std::get<1>(k);
                        //cout << "Found a particle with dx = " << dx << " and id = " << get<id>(b) << "\n";

                        //cout << "id of b " << get<id>(b) << endl;
                        //for (int i=0; i < particles.size(); i++) {
                        if (get<id>(b) != get<id>(particles[i])) { // check if it is not the same particle
                            //cout << "reject step " << 1 << endl;
                            free_position = false;
                        }
                        //}

                        //break;
                    }
                    //cout << "print position " << count_position << endl;
                    // check that the position they want to move to is free and not out of bounds
                    if (free_position == true && round((x[0] * (length_x / domain_length))) > 0 &&
                        round((x[0] * (length_x / domain_length))) < length_x - 1 && round(x[1]) > 0 &&
                        round(x[1]) < length_y - 1) {
                        get<position>(particles)[i] += vdouble2(sin(random_angle[0]),
                                                                cos(random_angle[0])); // update if nothing is in the next position
                    }

                }
                    // if second is greater than the first
                else if (new_chemo_1 < new_chemo_2) {
                    x += vdouble2(sin(random_angle[1]), cos(random_angle[1]));
                    //cout << "print id " << id_[x] << endl;


                    //cout << "Position "<< x << endl;
                    int count_position = 0;
                    bool free_position = true; // check if the neighbouring position is free

                    // if this loop is entered, it means that there is another cell where I want to move
                    for (const auto &k: euclidean_search(particles.get_query(), x, diameter)) {

                        count_position += 1; // just to check if this works
                        particle_type::const_reference b = std::get<0>(k);
                        const vdouble2 &dx = std::get<1>(k);
                        //cout << "Found a particle with dx = " << dx << " and id = " << get<id>(b) << "\n";

                        //cout << "id of b " << get<id>(b) << endl;
                        //for (int i=0; i < particles.size(); i++) {
                        if (get<id>(b) != get<id>(particles[i])) { // check if it is not the same particle
                            //cout << "reject step " << 1 << endl;
                            free_position = false;
                        }
                        //}

                        //break;
                    }

                    //cout << "print position " << count_position << endl;
                    // check that the position they want to move to is free and not out of bounds
                    if (free_position == true && round((x[0] * (length_x / domain_length))) > 0 &&
                        round((x[0] * (length_x / domain_length))) < length_x - 1 && round(x[1]) > 0 &&
                        round(x[1]) < length_y - 1) {
                        get<position>(particles)[i] += vdouble2(sin(random_angle[1]),
                                                                cos(random_angle[1])); // update if nothing is in the next position
                    }

                }


            }



            //cout << "sin of the angle " << sin(random_angle(rand_num_count)) << endl;
            //cout << "cos of the angle " << cos(random_angle(rand_num_count)) << endl;


            // check if the x coordinates are not out of the domain, if they are, ignore that step

            //cout << "chemo coord x before the test " << round(x[0]+sin(random_angle(rand_num_count))+sign_x*cell_radius) << endl;
            //cout << "chemo coord y before the test " << round(x[0]+cos(random_angle(rand_num_count))+sign_y*cell_radius) << endl;

            //rand_num_count += 3; // update random number count

        }// go through all the particles

        particles.update_positions();

        // save particles after they move
        /*
     * on every i/o step write particle container to a vtk
     * unstructured grid file
     */
        //cout << "." << flush;
#ifdef HAVE_VTK
        vtkWriteGrid("particles", t + 1, particles.get_grid(true));
#endif

        //for (int i =0;i<5;++i){cout << "koks rezas " << i <<endl;}
        //length_x_change = new_length_x_change;
        //cout << "new length " << new_length_x_change << endl;
    }//all time steps

    for (int i = 0; i < particles.size(); i++) {
        cout << "Positions = " << get<position>(particles[i]) << "\n";
    }

// count the number of particles that are at the last 10% of the domain

    /*
    double last_10_domain = domain_length - domain_length*domain_fraction;
    cout << "last 10 percent " << last_10_domain << endl;
    int number_of_cells = 0; // number of cells in the last 10% of the domain

    for (int p =0; p< particles.size(); p++){
        vdouble2 x = get<position>(particles[p]);
        if(x[0] > last_10_domain )
            number_of_cells += 1;
    }


    double proportion_cells_last = double(number_of_cells)/double(particles.size());

    cout << "proportion of cells " << proportion_cells_last << endl;

    return proportion_cells_last;
    */

    // return the furthest distance travelled by the cells
    double furthest_distance = 0 ;
    vdouble2 dist; // variable for positions

    for (int i = 0; i < particles.size(); i++) {

        dist = get<position>(particles[i]);

        if (furthest_distance < dist[0]){
            furthest_distance = dist[0];
        }

    }



    return furthest_distance;



}


int main(){

    const int number_parameters = 100; // parameter range

    const int sim_num = 40;

    MatrixXf all_distances = MatrixXf::Zero(number_parameters,sim_num); //matrix over which I am going to average

//n would correspond to different seeds
    for (int n = 0; n < sim_num; n++) {

        // define parameters that I will change
        //VectorXf slope, threshold;
        array<double, 1> lam;
        array<double, number_parameters> threshold;
        //array<double,number_parameters,number_parameters>;

        //MatrixXd density = MatrixXd::Zero(number_parameters, number_parameters); // need for because that is how paraview accepts data, third dimension is just zeros
        MatrixXf density = MatrixXf::Zero(number_parameters, 1);


        for (int i = 0; i < number_parameters; i++) {

            //lam[i] = 10 * (i + 1);//0.1;
            //threshold[i] = 0.5;//0.01*(i+1);// 0.01;
            //threshold[i] = 1*(i+1);// 0.01;
            threshold[i] = 0.005 * (i + 1);// 0.01;

        }

        lam[0] = 100;

        //lam[1] = 500;
        //cout << "lam " << lam[0] << endl;
        //cout << "threshold " << threshold[0] << endl;

        //cout << "value " << func(0.05,410.0) << endl;

        for (int i = 0; i < number_parameters; i++) {
            for (int j = 0; j < 1; j++) {
                cout << "iteration i " << i << endl;
                density(i, j) = func(threshold[i], lam[j]);


            }
            all_distances(i, n) = density(i, 0);
        }


        /*
        // save data to plot chemoattractant concentration
        ofstream output("density_matrix_cell_induced.csv");

        MatrixXd density_3col(number_parameters * number_parameters, 4), density_3col_ind(number_parameters * number_parameters,
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
        */




        // This might be useful for matlab
        ofstream output2("density_matrix_matlab_cell_induced.csv");

        for (int i = 0; i < number_parameters; i++) {

            for (int j = 0; j < 1; j++) {

                output2 << density(i, j) << ", ";

            }
            output2 << "\n" << endl;
        }


        ofstream output3("simulations_cell_induced.csv");

        for (int i = 0; i < number_parameters; i++) {

            for (int j = 0; j < sim_num; j++) {

                output3 << all_distances(i, j) << ", ";

            }
            output3 << "\n" << endl;
        }

    }

}