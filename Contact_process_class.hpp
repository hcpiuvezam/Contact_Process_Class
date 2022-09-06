/*
Programmer: Helena Christina Piuvezam de Albuquerque Bastos
Descrição: Program to reproduce figures from the article 
"Griffiths phases and the stretching of criticality in brain networks"
	from Paolo Moretti and Miguel Angel Muñoz 
*/

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>
#include <bits/stdc++.h>
#include <sstream>

using namespace std;

class Contact_Process{
	/* Creates an object that holds the simulation itself 
		The first parameter determines the size of the network
			if HMN1 it is the s in 2^s;
			otherwise it is the amount of sites.
		The second parameter determines the control parameter
			if model A lambda is any number
			if model B lambda is between 0 and 1
		The third parameter defines the model
		And, the forth parameter defines the network
			if 0 it's a random network with a fixec number of neighbors
			if 1 it's a fully connected network
			if 2 it opens the human connectome file
			if 3 it's a HMN-1 network
			and if 4 it opens the ER network file */

	/* neurons are identified with nodes of the network and are endowed 
		with a binary state-variable σ=0.1 */

	private:
		/* Private functions and variables */
		// Probabilities
		double prob_creation,
			   prob_inact,
			   r;
		// indices
		int x,
			y,
			a,
			j = 0;
		// Network name
		char name_net[80];

		vector< vector<int> > HMN1_ind(double p, double alpha, int s);
		vector< vector<int> > random_network(int N_sites, int N_neighbours);
		vector< vector <int> > open_file_to_matrix (char *name);

	public:
		/* Public functions and variables*/
		// Parameters
		double lambda,
		   	   mu = 1;
		// time
		double t;
		// counters
		int size;
		// Size of Simulation
		int L,
			S = 14,
			N_neighbours,
			s,
			active_num;
		// Network, model, and acquisition
		int model,
			fullnet,
			data_graph;
		// Vectors
		vector<int> cells,
					active_indexes,
					Network_line;
	
		std::vector< std::vector<int> > Network;

		// Constructor
		Contact_Process (int L_new, double lambda_new, int model_new, int fullnet_new);
		// Functions
		void advance_time();
		void simulation();
};

class Density_Decay : public Contact_Process{
	private:
		// time
		double dt,
			   t_max;
		bool time_log;
		// counters
		int	d,
			N_s = 0,
			N_trials,
			steps,
			counter;
		double N_inv = 0.;
		// auxiliar
		double dev_m = 0.,
			   rho_m = 0.,
			   aux_rho,
			   rho_mold = 0.;

		void time_loop();
		void trial_loop();
		vector<double> vector_log_time(double t_max);

	public:
		std::vector<double> rho_med,
							rho_std,
							t_med,
							rho_aux;

		Density_Decay (double t_max, int N_trials, bool time_log, int L_new, double lambda_new, int model_new, int fullnet_new);
};

class Avalanche_dynamics : public Contact_Process{};

class Dynamic_Range : public Contact_Process{};
