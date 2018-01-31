/*
cells move towords high concentration of varying chemoattracted, the movement is directed,
 cells only move if they sense higher concentration. Leaders and followers separately
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

int main() {


    // model parameters

    int length_x = 30;//240;//40;//4; // length of the chemoattractant vector, for fixed domain
    double domain_length = 30; //this variable is for the actual domain length
    double old_length = 30;// this is important for the update of the positions
    const int length_y = 12;//120;//20;//4;
    const double diameter = (2*7.5)/10;//2 // diameter in which there have to be no cells, equivalent to size of the cell
    double cell_radius = (7.5)/10;//0.5; // radius of a cell
    int N_steps = 200; // number of times the cells move up the gradient
    const size_t N = 4; // initial number of cells
    double l_filo = 27.5/10;//2; // sensing radius
    double diff_conc = 0.05; // sensing threshold, i.e. how much concentration has to be bigger, so that the cell moves in that direction
    int freq_growth = 1; // determines how frequently domain grows (actually not relevant because it will go every timestep)
    int insertion_freq = 1;
    double speed_l = 0.5; // speed of a leader cell
    double speed_f = 1; // speed of a follower cell

    // domain growth parameters


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

    double domain_len_der = 0; // initialise derivative of the domain growth function

    // parameters for the dynamics of chemoattractant concentration


    double D = 1/10; // to 10^5 \nu m^2/h diffusion coefficient
    double t = 0; // initialise time, redundant
    double dt = 0.1/10; // time step
    int dx = 1; // space step in x direction
    int dy = 1; // space step in y direction
    double kai = 1/10;//0.0001/10; // to 1 /h production rate of chemoattractant


    // parameters for internalisation

    double R = 7.5/10; // \nu m cell radius
    int lam = 100/10;//(100)/10; // to 1000 /h chemoattractant internalisation


    // matrix that stores the values of concentration of chemoattractant
    MatrixXf chemo(length_x, length_y), chemo_new(length_x,length_y);

    // initialise internalisation matrix
    MatrixXf intern(length_x,length_y);

    // generate initial chemoattractant concentration
    for (int i = 0;i<length_x;i++){
        for (int j = 0; j< length_y; j++){
            chemo(i,j) = 1; // uniform concentration initially
            chemo_new(i,j) = 1; // this is for later updates
        }
    }


    // four columns for x, y, z, u (z is necessary for paraview)

    // form a matrix which would store x,y,z,u

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
    ABORIA_VARIABLE(type,int,"type") // 1 if a cell is a leader, 0 if it is attached to another cell, 2 if dettached
    ABORIA_VARIABLE(direction,vdouble2,"direction")
    ABORIA_VARIABLE(attached_to_id,int,"attached_to_id")
    ABORIA_VARIABLE(attached_to_type,int,"attached_to_type")
    ABORIA_VARIABLE(chain,int,"chain") // 1 if attached to a chain with a leader in front, 0 if not attached
    typedef Particles<std::tuple<radius, direction>,2> particle_type; // 2 stands for dimension
    typedef Particles<std::tuple<radius, attached_to_id, attached_to_type, direction, chain>,2> followers_type;
    /*
     * if attached to a leader attached_to_type = 2;
     * if attached to a follower attached_to_type = 1;
     * if dettached attached_to_type = 0;
     * */
    //typedef Particles<std::tuple<radius, attached_to, >,2> attached_followers_type;
    //typedef Particles<std::tuple<radius,type>,2> detached_followers_type;

    //typedef Particles<std::tuple<>,2,std::vector,bucket_search_serial> particle_type;
    typedef particle_type::position position;
    typedef followers_type::position position;
    //typedef attached_followers_type::position position;
    //typedef detached_followers_type::position position;

    particle_type particles;
    followers_type followers;
    //attached_followers_type attached_followers;
    //detached_followers_type detached_followers;

    std::default_random_engine gen;
    std::uniform_real_distribution<double> uniform(2,length_y-1);

    /*
         * initialise neighbour search with 2d cuboid domain,
         * periodic in x and y
         */

    particles.init_neighbour_search(vdouble2(0,0), vdouble2(length_x,length_y), vbool2(false,false));
    followers.init_neighbour_search(vdouble2(0,0), vdouble2(length_x,length_y), vbool2(false,false));


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
    vtkWriteGrid("particles",t,particles.get_grid(true));


    // choose a set of random number between 0 and 2*pi, to avoid more rejections when it goes backwords (it would always be rejected)
    std::default_random_engine gen1;
    std::uniform_real_distribution<double> uniformpi(0, 2*M_PI);


    for (int t = 0; t < N_steps; t++){
        // insert new cells at the start of the domain at insertion time (have to think about this insertion time)

        if (t % insertion_freq == 0 ){
            bool free_position = false;
            followers_type::value_type f;
            get<radius>(f) = cell_radius;


            get<position>(f) = vdouble2(cell_radius,uniform(gen)); // x=2, uniformly in y
            free_position = true;
            /*
             * loop over all neighbouring leaders within "dem_diameter" distance
             */
            for (auto tpl: euclidean_search(particles.get_query(),get<position>(f),diameter)) {
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

            /*
             * loop over all neighbouring leaders within "diameter" distance
             */

            for (auto tpl: euclidean_search(followers.get_query(),get<position>(f),diameter)) {
                /*
                 * tpl variable is a tuple containing:
                 *  (0) -> neighbouring particle value_type
                 *  (1) -> relative position of neighbouring particle
                 *         from query point
                 */
                const vdouble2& dx = std::get<1>(tpl);
                const followers_type::value_type& j = std::get<0>(tpl);
                if (dx.norm() <  diameter) {
                    free_position = false;
                    break;
                }
            }


            // find the closest leader or other follower and set the attach to it

            // variables to keep track of distances
            particle_type::value_type closest_neighbour;
            followers_type::value_type closest_neighbour_follower;
            double distance = l_filo; // initially set to l_filo
            double distance_follower = l_filo;

            // closest leader
            for (auto tpl: euclidean_search(particles.get_query(),get<position>(f),l_filo)) {
                /*
                 * tpl variable is a tuple containing:
                 *  (0) -> neighbouring particle value_type
                 *  (1) -> relative position of neighbouring particle
                 *         from query point
                 */
                const vdouble2& dx = std::get<1>(tpl);
                const particle_type::value_type& j = std::get<0>(tpl);

                if (dx.norm() <  l_filo) {
                    // if the new distance is shorter than the previous shortest, set these new values
                    if (distance > dx.norm()) {
                        closest_neighbour = j;
                        distance = dx.norm();
                    }
                }
            }

            // closest follower

            for (auto tpl: euclidean_search(followers.get_query(),get<position>(f),l_filo)) {
                /*
                 * tpl variable is a tuple containing:
                 *  (0) -> neighbouring particle value_type
                 *  (1) -> relative position of neighbouring particle
                 *         from query point
                 */
                const vdouble2& dx = std::get<1>(tpl);
                const followers_type::value_type& j = std::get<0>(tpl);


                // make sure I do not look at the same particle, to avoid distance 0
                if (get<id>(j) != get<id>(f)) { // check if it is not the same particle
                    if (dx.norm() <  l_filo) {
                        // if the new distance is shorter than the previous shortest, set these new values
                        if (distance > dx.norm()) {
                            closest_neighbour_follower = j;
                            distance_follower = dx.norm();
                        }
                    }
                }
            }

            // compare if the closest leader or follower is the nearest neighbour to choose correct type
            cout << " distance " << distance << endl;
            cout << " follower distance " << distance_follower << endl;

            if (distance < distance_follower && distance != 0){
                get<attached_to_id>(f) = get<id>(closest_neighbour); // the id of a leader
                get<attached_to_type>(f) = 2; // attached to a leader
                get<chain>(f) = 1;
                cout << "enters here and sets correctly " << endl;
            }
            else if (distance > distance_follower && distance_follower != 0){
                get<attached_to_id>(f) = get<id>(closest_neighbour_follower); // the id of a leader
                get<attached_to_type>(f) = 1; // attached to a follower
                if (get<chain>(closest_neighbour_follower) == 1){
                    get<chain>(f) = 1;
                }
                else
                {
                    get<chain>(f) = 0; // not part of a chain, this is actually equivalent as detached
                }
            }
            else{
                get<attached_to_id>(f) = 0; // the id of a leader
                get<attached_to_type>(f) = 0; // dettached
                get<chain>(f) = 0;
            }
            cout << "attached to " << get<attached_to_type>(f) << endl;

            if (free_position == true){
                followers.push_back(f);}

        }


        /////////////////////////////////////
        // grow domain


        if (t % freq_growth == 0){

            //domain_length = domain_length + 1.0;

            // no time delay
            //domain_length = ((L_inf*exp(a*t))/ (L_inf/L_0 + exp(a*t) - 1) );

            //domain_len_der = ((a*L_inf*exp(a*t))/ (L_inf/L_0 + exp(a*t) - 1) - (a*L_inf*exp(2*a*t))/(L_inf/L_0 + exp(a*t) - 1) );

            // from the paper
            //domain_length = ((L_inf*exp(a*(t-t_s)))/ (L_inf/L_0 + exp(a*(t-t_s)) - 1) ) + constant;

            //domain_len_der = ((a*L_inf*exp(a*(t-t_s)))/ (L_inf/L_0 + exp(a*(t-t_s)) - 1) - (a*L_inf*exp(2*a*(t-t_s)))/(L_inf/L_0 + exp(a*(t-t_s)) - 1) );

            // with time delay and constant to make the initial conditions consistent
            domain_length = ((L_inf*exp(a*(t-t_s)))/ (L_inf/L_0 + exp(a*(t-t_s)) - 1) ) + constant;

            domain_len_der = ((a*L_inf*exp(a*(t-t_s)))/ (L_inf/L_0 + exp(a*(t-t_s)) - 1) - (a*L_inf*exp(2*a*(t-t_s)))/(L_inf/L_0 + exp(a*(t-t_s)) - 1) );

        }
        //cout << "diff domain outside " << diff_domain << endl;

        // update chemoattractant profile




        // internalisation
        for (int i =0;i<length_x;i++){
            for(int j = 0;j<length_y;j++){
                //go through all the cells
                // leaders
                for (int k =0; k<particles.size();k++){
                    vdouble2 x;
                    x = get<position>(particles[k]);

                    intern(i,j) = intern(i,j) + exp(- (((domain_length/length_x)*i-x[0])*((domain_length/length_x)*i-x[0])+(j-x[1])*(j-x[1]))/(2*R*R)); // mapping to fixed domain
                }
                //followers
                for (int k =0; k<followers.size();k++){
                    vdouble2 x;
                    x = get<position>(followers[k]);

                    intern(i,j) = intern(i,j) + exp(- (((domain_length/length_x)*i-x[0])*((domain_length/length_x)*i-x[0])+(j-x[1])*(j-x[1]))/(2*R*R)); // mapping to fixed domain
                }
            }
        }

        for (int i=1;i<length_x-1;i++){
            for (int j=1;j<length_y-1;j++){


                // logistic production rate
                //chemo_new(i,j) = dt * (D*((1/((domain_length/length_x)*(domain_length/length_x))) * (chemo(i+1,j)-2*chemo(i,j)+chemo(i-1,j))/(dx*dx) + (chemo(i,j+1)- 2* chemo(i,j)+chemo(i,j-1))/(dy*dy)) - (chemo(i,j)*lam / (2*M_PI*R*R)) * intern(i,j) + kai*chemo(i,j)*(1-chemo(i,j)) - double(domain_len_der)/double(domain_length) * chemo(i,j) ) + chemo(i,j);


                //different production rate, linear
                //chemo_new(i,j) = dt * (D*((1/((domain_length/length_x)*(domain_length/length_x))) * (chemo(i+1,j)-2*chemo(i,j)+chemo(i-1,j))/(dx*dx) + (chemo(i,j+1)- 2* chemo(i,j)+chemo(i,j-1))/(dy*dy)) - (chemo(i,j)*lam / (2*M_PI*R*R)) * intern(i,j) + kai*chemo(i,j)*(1-chemo(i,j)) - double(domain_len_der)/double(domain_length) * chemo(i,j) ) + chemo(i,j);

                // different production rate, threshold value

                /*if (chemo(i,j)<=1){
                    chemo_new(i,j) = dt * (D*((1/((domain_length/length_x)*(domain_length/length_x))) * (chemo(i+1,j)-2*chemo(i,j)+chemo(i-1,j))/(dx*dx) + (chemo(i,j+1)- 2* chemo(i,j)+chemo(i,j-1))/(dy*dy)) - (chemo(i,j)*lam / (2*M_PI*R*R)) * intern(i,j) + kai - double(domain_len_der)/double(domain_length) * chemo(i,j) ) + chemo(i,j);
                }
                else{
                    chemo_new(i,j) = dt * (D*((1/((domain_length/length_x)*(domain_length/length_x))) * (chemo(i+1,j)-2*chemo(i,j)+chemo(i-1,j))/(dx*dx) + (chemo(i,j+1)- 2* chemo(i,j)+chemo(i,j-1))/(dy*dy)) - (chemo(i,j)*lam / (2*M_PI*R*R)) * intern(i,j) - double(domain_len_der)/double(domain_length) * chemo(i,j) ) + chemo(i,j);
                }*/


                // no source
                chemo_new(i,j) = dt * (D*((1/((domain_length/length_x)*(domain_length/length_x))) * (chemo(i+1,j)-2*chemo(i,j)+chemo(i-1,j))/(dx*dx) + (chemo(i,j+1)- 2* chemo(i,j)+chemo(i,j-1))/(dy*dy)) - (chemo(i,j)*lam / (2*M_PI*R*R)) * intern(i,j)  - double(domain_len_der)/double(domain_length) * chemo(i,j) ) + chemo(i,j);

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


        //cout << "old length " << old_length << endl;
        //cout << "domain length " << domain_length << endl;
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
            for (int i = 0; i< followers.size();i++){
                get<position>(followers)[i] *= vdouble2((domain_length/old_length), 1);
            }
            old_length = domain_length;
        }


        /////////////////////////////////////


        // Update positions based on the gradient




        // SEARCH MULTIPLE TIMES

        // move all the leaders

        for (int i = 0; i < particles.size(); i++) {


            vdouble2 x;
            x = get<position>(particles[i]);
            //cout << "particles.size " << i << endl;
            //cout << "print id " << get<id>(particles[i]) << endl;
            //cout << "x coord " << x[0] << endl;


            // create an array to store random directions
            std::array<double, 3> random_angle;
            std::array<int, 3> sign_x;
            std::array<int, 3> sign_y;
            for (int i = 0; i < 3; i++) {

                double random_angle_tem = uniformpi(gen1);
                int sign_x_tem, sign_y_tem;

                while (round((x[0] * (length_x / domain_length) + sin(random_angle_tem) + sign_x_tem * l_filo)) < 0 ||
                       round((x[0] * (length_x / domain_length) + sin(random_angle_tem) + sign_x_tem * l_filo)) >
                       length_x - 1 || round(x[1] + cos(random_angle_tem) + sign_y_tem * l_filo) < 0 ||
                       round(x[1] + cos(random_angle_tem) + sign_y_tem * l_filo) > length_y - 1) {
                    random_angle_tem = uniformpi(gen1);

                    if (sin(random_angle_tem) < 0) {
                        sign_x_tem = -1;
                    } else { sign_x_tem = 1; }

                    if (cos(random_angle_tem) < 0) {
                        sign_y_tem = -1;
                    } else { sign_y_tem = 1; }

                }

                random_angle[i] = random_angle_tem;

                sign_x[i] = sign_x_tem;
                sign_y[i] = sign_y_tem;


            }


            // choose which direction to move

            // store variables for concentration at new locations


            double old_chemo = chemo((round((x)[0] * (length_x / domain_length))),round(x)[1]);

            double new_chemo_1 = chemo(round((x[0] * (length_x / domain_length) + sin(random_angle[0]) + sign_x[0] * l_filo)),
                                       round(x[1] + cos(random_angle[0]) + sign_y[0] * l_filo));

            double new_chemo_2 = chemo(round((x[0] * (length_x / domain_length) + sin(random_angle[1]) + sign_x[1] * l_filo)),
                                       round(x[1] + cos(random_angle[1]) + sign_y[1] * l_filo));


            //if both smaller, move random direction
            //absolute
            if (new_chemo_1 - old_chemo < diff_conc && new_chemo_2 - old_chemo < diff_conc) {


                // relative
                //if ((new_chemo_1 - old_chemo)/sqrt(old_chemo) < diff_conc && (new_chemo_2- old_chemo)/sqrt(old_chemo) < diff_conc){

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

                    //for (int i=0; i < particles.size(); i++) {
                    if (get<id>(b) != get<id>(particles[i])) { // check if it is not the same particle
                        //cout << "reject step " << 1 << endl;
                        free_position = false;
                    }
                }

                for (const auto &k: euclidean_search(followers.get_query(), x, diameter)) {

                    followers_type::const_reference b = std::get<0>(k);
                    const vdouble2 &dx = std::get<1>(k);
                    //cout << "Found a particle with dx = " << dx << " and id = " << get<id>(b) << "\n";

                    //for (int i=0; i < particles.size(); i++) {
                    //if (get<id>(b) != get<id>(followers[i])) { // check if it is not the same particle
                        //cout << "reject step " << 1 << endl;
                        free_position = false;
                    //}
                }


                //cout << "print position " << count_position << endl;

                // check that the position they want to move to is free and not out of bounds
                if (free_position == true && round((x[0] * (length_x / domain_length) + sin(random_angle[2]) )) > 0 &&
                        round((x[0] * (length_x / domain_length) + sin(random_angle[2]))) < length_x - 1 && round(x[1] + cos(random_angle[2])) > 0 &&
                        round(x[1] + cos(random_angle[2])) < length_y - 1) {
                    get<position>(particles)[i] += speed_l*vdouble2(sin(random_angle[2]),
                                                            cos(random_angle[2])); // update if nothing is in the next position
                    get<direction>(particles)[i] = speed_l*(sin(random_angle[2]),
                                                            cos(random_angle[2]));
                }

            }
                //cout << "stops here " << endl;
                // if first direction greater, second smaller
                //absolute
            else if (new_chemo_1 - old_chemo > diff_conc && new_chemo_2 - old_chemo < diff_conc){

                //relative
                //else if ((new_chemo_1 - old_chemo)/sqrt(old_chemo) > diff_conc && (new_chemo_2 - old_chemo)/sqrt(old_chemo) < diff_conc){

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

                    //for (int i=0; i < particles.size(); i++) {
                    if (get<id>(b) != get<id>(particles[i])) { // check if it is not the same particle
                        //cout << "reject step " << 1 << endl;
                        free_position = false;
                    }

                }

                for (const auto &k: euclidean_search(followers.get_query(), x, diameter)) {

                    followers_type::const_reference b = std::get<0>(k);
                    const vdouble2 &dx = std::get<1>(k);
                    //cout << "Found a particle with dx = " << dx << " and id = " << get<id>(b) << "\n";

                    //for (int i=0; i < particles.size(); i++) {
                    //if (get<id>(b) != get<id>(followers[i])) { // check if it is not the same particle
                    //cout << "reject step " << 1 << endl;
                    free_position = false;
                    //}
                }



                //cout << "print position " << count_position << endl;
                // check that the position they want to move to is free and not out of bounds
                if (free_position == true && round((x[0] * (length_x / domain_length) + sin(random_angle[0]) )) > 0 &&
                    round((x[0] * (length_x / domain_length) + sin(random_angle[0]))) < length_x - 1 && round(x[1] + cos(random_angle[0])) > 0 &&
                    round(x[1] + cos(random_angle[0])) < length_y - 1){
                    get<position>(particles)[i] += speed_l*vdouble2(sin(random_angle[0]),
                                                            cos(random_angle[0])); // update if nothing is in the next position
                    get<direction>(particles)[i] = speed_l*(sin(random_angle[0]),
                            cos(random_angle[0]));
                }

            }
                // if first smaller, second bigger

                //absolute
            else if (new_chemo_1 - old_chemo < diff_conc && new_chemo_2 - old_chemo > diff_conc){

                //relative
                //else if ((new_chemo_1 - old_chemo)/sqrt(old_chemo) < diff_conc && (new_chemo_2 - old_chemo)/sqrt(old_chemo) > diff_conc){



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

                    //for (int i=0; i < particles.size(); i++) {
                    if (get<id>(b) != get<id>(particles[i])) { // check if it is not the same particle
                        //cout << "reject step " << 1 << endl;
                        free_position = false;
                    }

                }


                for (const auto &k: euclidean_search(followers.get_query(), x, diameter)) {

                    followers_type::const_reference b = std::get<0>(k);
                    const vdouble2 &dx = std::get<1>(k);
                    //cout << "Found a particle with dx = " << dx << " and id = " << get<id>(b) << "\n";

                    //for (int i=0; i < particles.size(); i++) {
                    //if (get<id>(b) != get<id>(followers[i])) { // check if it is not the same particle
                    //cout << "reject step " << 1 << endl;
                    free_position = false;
                    //}
                }



                //cout << "print position " << count_position << endl;
                // check that the position they want to move to is free and not out of bounds
                if (free_position == true && round((x[0] * (length_x / domain_length) + sin(random_angle[1]) )) > 0 &&
                     round((x[0] * (length_x / domain_length) + sin(random_angle[1]))) < length_x - 1 && round(x[1] + cos(random_angle[1])) > 0 &&
                     round(x[1] + cos(random_angle[1])) < length_y - 1){
                    get<position>(particles)[i] += speed_l*vdouble2(sin(random_angle[1]),
                                                            cos(random_angle[1])); // update if nothing is in the next position
                    get<direction>(particles)[i] = speed_l*(sin(random_angle[1]),
                            cos(random_angle[1]));
                }
                break;
            }
                // if both greater choose the bigger one

                // absolute
            else if (new_chemo_1 - old_chemo > diff_conc && new_chemo_2 - old_chemo > diff_conc){

                //relative
                //else if ((new_chemo_1 - old_chemo)/sqrt(old_chemo) > diff_conc && (new_chemo_2 - old_chemo)/sqrt(old_chemo) > diff_conc){


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

                        //for (int i=0; i < particles.size(); i++) {
                        if (get<id>(b) != get<id>(particles[i])) { // check if it is not the same particle
                            //cout << "reject step " << 1 << endl;
                            free_position = false;
                        }

                    }


                    for (const auto &k: euclidean_search(followers.get_query(), x, diameter)) {

                        followers_type::const_reference b = std::get<0>(k);
                        const vdouble2 &dx = std::get<1>(k);
                        //cout << "Found a particle with dx = " << dx << " and id = " << get<id>(b) << "\n";

                        //for (int i=0; i < particles.size(); i++) {
                        //if (get<id>(b) != get<id>(followers[i])) { // check if it is not the same particle
                        //cout << "reject step " << 1 << endl;
                        free_position = false;
                        //}
                    }





                    //cout << "print position " << count_position << endl;
                    // check that the position they want to move to is free and not out of bounds
                    if (free_position == true && round((x[0] * (length_x / domain_length) + sin(random_angle[0]) )) > 0 &&
                        round((x[0] * (length_x / domain_length) + sin(random_angle[0]))) < length_x - 1 && round(x[1] + cos(random_angle[0])) > 0 &&
                        round(x[1] + cos(random_angle[0])) < length_y - 1){
                        get<position>(particles)[i] += speed_l*vdouble2(sin(random_angle[0]),
                                                                cos(random_angle[0])); // update if nothing is in the next position
                        get<direction>(particles)[i] = speed_l*(sin(random_angle[0]),
                                cos(random_angle[0]));
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

                        //for (int i=0; i < particles.size(); i++) {
                        if (get<id>(b) != get<id>(particles[i])) { // check if it is not the same particle
                            //cout << "reject step " << 1 << endl;
                            free_position = false;
                        }

                    }

                    for (const auto &k: euclidean_search(followers.get_query(), x, diameter)) {

                        followers_type::const_reference b = std::get<0>(k);
                        const vdouble2 &dx = std::get<1>(k);
                        //cout << "Found a particle with dx = " << dx << " and id = " << get<id>(b) << "\n";

                        //for (int i=0; i < particles.size(); i++) {
                        //if (get<id>(b) != get<id>(followers[i])) { // check if it is not the same particle
                        //cout << "reject step " << 1 << endl;
                        free_position = false;
                        //}
                    }



                    //cout << "print position " << count_position << endl;
                    // check that the position they want to move to is free and not out of bounds
                    if (free_position == true && round((x[0] * (length_x / domain_length) + sin(random_angle[1]) )) > 0 &&
                        round((x[0] * (length_x / domain_length) + sin(random_angle[1]))) < length_x - 1 && round(x[1] + cos(random_angle[1])) > 0 &&
                        round(x[1] + cos(random_angle[1])) < length_y - 1) {
                        get<position>(particles)[i] += speed_l*vdouble2(sin(random_angle[1]),
                                                                cos(random_angle[1])); // update if nothing is in the next position
                        get<direction>(particles)[i] = speed_l*(sin(random_angle[1]),
                                cos(random_angle[1]));
                    }
                }
            }
        }// go through all the particles-leaders


        // update the positions of all followers

        for (int i = 0; i < followers.size(); i++) {

            cout << "next time attached to " << get<attached_to_id>(followers[i]) << endl;
            vdouble2 x;
            x = get<position>(followers[i]);

            // check what the closest neighbour is

            // variables to keep track of distances
            particle_type::value_type closest_neighbour;
            followers_type::value_type closest_neighbour_follower;
            double distance = l_filo; // initially set to l_filo
            double distance_follower = l_filo;

            // closest leader
            for (auto tpl: euclidean_search(particles.get_query(),get<position>(followers[i]),l_filo)) {
                /*
                 * tpl variable is a tuple containing:
                 *  (0) -> neighbouring particle value_type
                 *  (1) -> relative position of neighbouring particle
                 *         from query point
                 */
                const vdouble2& dx = std::get<1>(tpl);
                const particle_type::value_type& j = std::get<0>(tpl);

                if (dx.norm() <  l_filo) {
                    // if the new distance is shorter than the previous shortest, set these new values
                    if (distance > dx.norm()) {
                        closest_neighbour = j;
                        distance = dx.norm();
                    }
                }
            }

            // closest follower

            for (auto tpl: euclidean_search(followers.get_query(),get<position>(followers[i]),l_filo)) {
                /*
                 * tpl variable is a tuple containing:
                 *  (0) -> neighbouring particle value_type
                 *  (1) -> relative position of neighbouring particle
                 *         from query point
                 */
                const vdouble2& dx = std::get<1>(tpl);
                const followers_type::value_type& j = std::get<0>(tpl);

                // make sure I do not consider the same particle
                if (get<id>(j) != get<id>(followers[i])) {
                    if (dx.norm() < l_filo) {
                        // if the new distance is shorter than the previous shortest, set these new values
                        if (distance > dx.norm()) {
                            closest_neighbour_follower = j;
                            distance_follower = dx.norm();
                        }
                    }
                }
            }

            // compare if the closest leader or follower is the nearest neighbour to choose correct type
           int id_main = get<id>(followers[i]);
           int id_nbh = get<id>(closest_neighbour);

            cout << "id of a particle  " << id_main << endl;
            cout << "id of the closest neighbour " << id_nbh << endl;
            cout << " distance " << distance << endl;
            cout << " follower distance " << distance_follower << endl;


            if (distance < distance_follower && distance != 0  ){
                get<attached_to_id>(followers[i]) = get<id>(closest_neighbour); // the id of a leader
                get<attached_to_type>(followers[i]) = 2; // attached to a leader
                get<chain>(followers[i]) = 1;
                cout << "cia irgi niekada?????????????????" << endl;
            }
            else if (distance > distance_follower && distance_follower != 0){
                get<attached_to_id>(followers[i]) = get<id>(closest_neighbour_follower); // the id of a leader
                get<attached_to_type>(followers[i]) = 1; // attached to a follower

                if (get<attached_to_id>(closest_neighbour_follower) == 2){
                    get<chain>(followers[i]) = 1;
                }
                else
                {
                    get<chain>(followers[i]) = 0; // not part of a chain, this is actually equivalent as detached
                }
            }
            else{ // if it was greater than t
                get<attached_to_id>(followers[i]) = 0; //
                get<attached_to_type>(followers[i]) = 0; // dettached
                get<chain>(followers[i]) = 0;
            }

            cout << "attached to type " << get<attached_to_type>(followers[i]) << endl;
            cout << "distance later " << distance << endl;

            // take the direction of the leader it is attached to

            if (get<attached_to_type>(followers)[i] == 2) {

                int id_lead = get<attached_to_id>(followers)[i]; // id of the leader it is attached to

                vdouble2 direction_cur; // current direction of a cell

                //access the direction of the particle it is attached to
                for (int i = 0; i < particles.size(); i++) {
                    if (get<id>(particles)[i] == id_lead) {
                        direction_cur = get<direction>(particles)[i];
                    }
                }

                // temporarily update the position

                x += speed_f*direction_cur;

                // now check if there is no other particle/follower in the other position
                // check if it is not out of bounds
                int count_position = 0;
                bool free_position = true; // check if the neighbouring position is free

                // if this loop is entered, it means that there is another cell where I want to move
                for (const auto &k: euclidean_search(particles.get_query(), x, diameter)) {

                    count_position += 1; // just to check if this works
                    particle_type::const_reference b = std::get<0>(k);
                    const vdouble2 &dx = std::get<1>(k);
                    //cout << "Found a particle with dx = " << dx << " and id = " << get<id>(b) << "\n";

                    //for (int i=0; i < particles.size(); i++) {
                    // if (get<id>(b) != get<id>(particles[i])) { // check if it is not the same particle
                        //cout << "reject step " << 1 << endl;
                        free_position = false;
                    //}
                }

                for (const auto &k: euclidean_search(followers.get_query(), x, diameter)) {

                    count_position += 1; // just to check if this works
                    followers_type::const_reference b = std::get<0>(k);
                    const vdouble2 &dx = std::get<1>(k);
                    //cout << "Found a particle with dx = " << dx << " and id = " << get<id>(b) << "\n";


                    //for (int i=0; i < particles.size(); i++) {
                    if (get<id>(b) != get<id>(followers[i])) { // check if it is not the same particle
                        //cout << "reject step " << 1 << endl;
                        free_position = false;
                    }
                }

                cout << " ar cia kada buna " << endl;
                //cout << "print position " << count_position << endl;
                // check that the position they want to move to is free and not out of bounds
                if (free_position == true && round((x[0] * (length_x / domain_length))) > 0 &&
                    round((x[0] * (length_x / domain_length))) < length_x - 1 && round(x[1]) > 0 &&
                    round(x[1]) < length_y - 1) {
                    get<position>(followers)[i] += speed_f*direction_cur; // update if nothing is in the next position
                    get<direction>(followers)[i] = speed_f*direction_cur;
                }

            }

            // if attahced to a follower which is a part of a chain
            else if (get<attached_to_type>(followers)[i] == 1 && get<chain>(followers)[i] == 1) {

                int id_fol = get<attached_to_id>(followers)[i]; // id of the follower it is attached to

                vdouble2 direction_cur; // current direction of a cell

                //access the direction of the follower it is attached to
                for (int i = 0; i < followers.size(); i++) {
                    if (get<id>(followers)[i] == id_fol) {
                        direction_cur = get<direction>(followers)[i];
                    }
                }

                // temporarily update the position

                x += speed_f*direction_cur;

                // now check if there is no other particle/follower in the other position
                // check if it is not out of bounds
                int count_position = 0;
                bool free_position = true; // check if the neighbouring position is free

                // if this loop is entered, it means that there is another cell where I want to move
                for (const auto &k: euclidean_search(particles.get_query(), x, diameter)) {

                    count_position += 1; // just to check if this works
                    particle_type::const_reference b = std::get<0>(k);
                    const vdouble2 &dx = std::get<1>(k);
                    //cout << "Found a particle with dx = " << dx << " and id = " << get<id>(b) << "\n";


                    //for (int i=0; i < particles.size(); i++) {
                    if (get<id>(b) != get<id>(particles[i])) { // check if it is not the same particle
                        //cout << "reject step " << 1 << endl;
                        free_position = false;
                    }
                    //}

                    //break;
                }

                for (const auto &k: euclidean_search(followers.get_query(), x, diameter)) {

                    count_position += 1; // just to check if this works
                    followers_type::const_reference b = std::get<0>(k);
                    const vdouble2 &dx = std::get<1>(k);
                    //cout << "Found a particle with dx = " << dx << " and id = " << get<id>(b) << "\n";


                    //for (int i=0; i < particles.size(); i++) {
                    if (get<id>(b) != get<id>(followers[i])) { // check if it is not the same particle
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
                    get<position>(followers)[i] += speed_f * direction_cur; // update if nothing is in the next position
                    get<direction>(followers)[i] = speed_f * direction_cur;
                }

            }

                // if the cell is not attached to anything, choose a random direction
            else {

                double random_angle = uniformpi(gen1);
                int sign_x, sign_y;

                while (round((x[0] * (length_x / domain_length) + sin(random_angle) + sign_x * l_filo)) < 0 ||
                       round((x[0] * (length_x / domain_length) + sin(random_angle) + sign_x * l_filo)) >
                       length_x - 1 || round(x[1] + cos(random_angle) + sign_y * l_filo) < 0 ||
                       round(x[1] + cos(random_angle) + sign_y * l_filo) > length_y - 1) {
                    random_angle = uniformpi(gen1);

                    if (sin(random_angle) < 0) {
                        sign_x = -1;
                    } else { sign_x = 1; }

                    if (cos(random_angle) < 0) {
                        sign_y = -1;
                    } else { sign_y = 1; }

                }

                x += vdouble2(sin(random_angle), cos(random_angle));

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

                for (const auto &k: euclidean_search(followers.get_query(), x, diameter)) {

                    count_position += 1; // just to check if this works
                    followers_type::const_reference b = std::get<0>(k);
                    const vdouble2 &dx = std::get<1>(k);
                    //cout << "Found a particle with dx = " << dx << " and id = " << get<id>(b) << "\n";


                    //for (int i=0; i < particles.size(); i++) {
                    if (get<id>(b) != get<id>(followers[i])) { // check if it is not the same particle
                        //cout << "reject step " << 1 << endl;
                        free_position = false;
                    }
                }


                //cout << "print position " << count_position << endl;

                // check that the position they want to move to is free and not out of bounds
                if (free_position == true && round((x[0] * (length_x / domain_length))) > 0 &&
                    round((x[0] * (length_x / domain_length))) < length_x - 1 && round(x[1]) > 0 &&
                    round(x[1]) < length_y - 1) {
                    cout << "how frequently come in here " << endl;
                    get<position>(followers)[i] += speed_f * vdouble2(sin(random_angle), cos(random_angle)); // update if nothing is in the next position
                    get<direction>(followers)[i] = speed_f * vdouble2(sin(random_angle), cos(random_angle));
                }
            }


        }


    #ifdef HAVE_VTK
            vtkWriteGrid("particles",t,particles.get_grid(true));
    #endif
    #ifdef HAVE_VTK
            vtkWriteGrid("followers",t,followers.get_grid(true));
    #endif
    }


}
