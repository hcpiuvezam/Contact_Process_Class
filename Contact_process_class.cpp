#include "Contact_process_class.hpp"

using namespace std;

/* NETWORK FUNCTIONS */

vector< vector< int > > Contact_Process::HMN1_ind(double p, double alpha, int s){
	
    /* Function that builds a HMN-1 (Hierarchical Modular Network) */

	int M0 = 2,
		b = 2,
		N = M0 * pow(b, s),
		l = 1,
		conections = 1,
		block_size = 1;
	
	double r,
		   prob;
			
	std::vector<int> aux,
					 aux_pair;
	std::vector< std::vector<int> > Network,
									Pairs,
									Pairs_new;
	
	for(int i = 0; i < N; i += b){
		
		for(int j = b-1; j >= 0; --j){
		aux.push_back(i + j);
		aux_pair.push_back(i + j);
		Network.push_back(aux);
		aux.clear();
		}
		
		Pairs.push_back(aux_pair);
		aux_pair.clear();
	}
	
	while(l <= s){
		prob = ((double) alpha * pow(p, l));
		
		for(size_t i = 0; i < Pairs.size(); i += b){
			for(int j = 0; j < b; ++j)
				for(size_t m = 0; m < Pairs[i + j].size(); ++m)
					aux_pair.push_back(Pairs[i + j][m]);
					
			Pairs_new.push_back(aux_pair);
			aux_pair.clear();
			
			int a = 0;
			
			while(a == 0){
				for(size_t j = 0; j < Pairs[i].size(); ++j){
					for(size_t k = 0; k < Pairs[i+1].size(); ++k){
						if(not(std::find(Network[Pairs[i][j]].begin(), Network[Pairs[i][j]].end(), Pairs[i+1][k]) != Network[Pairs[i][j]].end())){
							r = ((double) rand()/RAND_MAX);
							if ( r < prob ){
								a = 1;
								Network[Pairs[i][j]].push_back(Pairs[i+1][k]);
								Network[Pairs[i+1][k]].push_back(Pairs[i][j]);
							}
						}
					}
				}
			}
			
		}
		Pairs.clear();
		Pairs = Pairs_new;
		Pairs_new.clear();
		
		l++;
	}
	
	/*for(size_t i = 0; i < Network.size(); ++i){
        cout << i << '\t';
		for(size_t j = 0; j < Network[i].size(); ++j)
			cout << Network[i][j] << ' ';
		cout << '\n';
	}*/
	
	return Network;
}

vector< vector<int> > Contact_Process::random_network(int N_sites, int N_neighbours){

	/* Function that builds a Random Neural Network with k neighbours */

	int r = 0, m = 0;
	vector<int> indices (N_neighbours, 0);
	vector< std::vector<int> > Network;
	
	if(N_neighbours >= N_sites){
		cout << "Try Again!\n";
		exit(0);
	}
	
	for (int i = 0; i < N_sites; ++i){
		
		for (int j = 0; j < N_neighbours; ++j){
			
			m = 0;
			
			do {
				r = rand() % N_sites;
				for(int k = 0; k <= j; ++k){
					if(r == indices[k]){
						m = 1;
						}
				}
			} while (r == i and m == 1);
			
			indices[j] = r;
		}
		Network.push_back(indices);
	}
	
	return Network;
}

vector< vector < int > > Contact_Process::open_file_to_matrix (char *name){

	/* Function that opens a file with the neighbours indices and save it to a vector of vectors  */
	
	string line;
	int data;
	
	vector<int> Network_line;
	vector< std::vector<int> > Network;
	
	ifstream file(name);
	
	if (file.is_open())
	{
	  while ( getline (file,line, ' ') )
	  {
	  	if(line.size() > 0){
	      istringstream(line) >> data;
	      Network_line.push_back(data);
	      
	    	for (size_t i = 0; i < line.length(); i++)
	   			if (line[i] == '\n'){
	   				Network.push_back(Network_line);
	   				Network_line.clear();
	   			}
	    }
	  }
	  file.close();
	}

	else cout << "Unable to open file";
	
	return Network;
}

/* Constructors */

Contact_Process::Contact_Process(int L_new, double lambda_new, int model_new, int fullnet_new)
: L(L_new), lambda(lambda_new), model(model_new), fullnet(fullnet_new) {
    switch(fullnet_new){
        case 0:
            Network = random_network(L, N_neighbours);
            break;
        case 1:
            /* Fully Connected Netwrok*/
			for(int i = 0; i < L; ++i){
				for(int l = 0; l < L; ++l)
					Network_line.push_back(l);
                Network.push_back(Network_line);
                Network_line.clear();
			}
            break;
        case 2:
            /* Human Connectome */	
			sprintf(name_net,"Networks/adj_ind_Hagmann_1000.txt");
			Network = open_file_to_matrix(name_net);
            break;
        case 3:
            /* 2 ** S Nodes in a Hyerarchical Modular Network */
            S = L;
			Network = HMN1_ind(.25, 1., S - 1);//open_file_to_matrix (name2); // 
			L = ((int) 2 * pow(2, S - 1));
            break;
        case 4:
			/* 10 ** 6 Nodes in an Erdos-Rényi Network */
			sprintf(name_net,"Networks/adj_ind_ER_10_6_n_20.txt");
			Network = open_file_to_matrix(name_net);
            L = pow(10, 6);
            break;
    }
    for(int i = 0; i < L; ++i)
        cells.push_back(0);

    prob_creation = (double) lambda / (mu + lambda);
	prob_inact = (double) mu / (mu + lambda); /* becomes inactive σ=0 with probability μ/(λ+μ) */

    t = 0.;
}

Density_Decay::Density_Decay (double t_max, int N_trials, bool time_log, int L_new, double lambda_new, int model_new, int fullnet_new)
: Contact_Process(L_new, lambda_new, model_new, fullnet_new){
	for(int i = 0; i < L; ++i)
		/* Initializes the network at all active */
        cells[i] = 1;

	active_num = L;
	
	if (time_log){
		t_med = vector_log_time(t_max);
		steps = (int) t_med.size();
	}
	else{
		dt = 0.1;
		steps = (int) t_max / dt;
		for(int n = 0; n < steps + 1; n++)
			t_med[n] = (double) n * dt;
	}

	for(std::size_t n = 0; n < steps + 1; n++){
		rho_med.push_back(0.);
		rho_aux.push_back(0.);
		rho_std.push_back(0.);
	}
}


/* Simulation functions */

void Contact_Process::advance_time(){
    t += (double) 1 / active_num; // Rejection-free Monte Carlo time
	simulation();
}

void Contact_Process::simulation(){
    /* Pick from active sites */
    
    a = (int) rand() % active_num;
    x = active_indexes[a]; /* At each step, an active node is selected */
    
    N_neighbours = (int) Network[x].size(); // Amount of neighbours the site x has.

    r = (double) rand() / RAND_MAX;

    if(r < prob_inact){
        /*Deactivate site x with probability mu / (mu + lambda)*/
        
        cells[x] = 0;

        active_num -= 1;
        active_indexes.erase(active_indexes.begin() + a); //becomes inactive σ=0 with probability μ/(λ+μ)
        
    }
    else{
        switch(model){
            case 0:
            /* Activation model A */
            
                y = (int) rand() % N_neighbours;
                y = Network[x][y]; /* it activates one randomly chosen nearest neighbour provided it was inactive */

                if(cells[y] == 0){
                    size += 1;
                    active_indexes.push_back(y);
                    active_num += 1;
                }

                cells[y] = 1;
                break;

            case 1:
                /* Activation model B */
                
                for(int k = 0; k < N_neighbours; ++k){ /* it checks all of its nearest neighbours */
                    
                    y = (int) Network[x][k];
                        
                    if(cells[y] == 0){ /* provided it was inactive */
                    
                        r = ((double) rand()/RAND_MAX);
                        
                        if(r < lambda){ /* activating each of them with probability 0<λ<1 */

                            cells[y] = 1;
                            
                            /*if(data_graph == 0)
                                outfile2 << dt << ' ' << y << '\n';*/
                            
                            active_indexes.push_back(y);
                            active_num += 1;
                            size += 1;
                        }
                    }
                }

                cells[x] = 0;
                
                active_num -= 1;
                active_indexes.erase(active_indexes.begin()+a); /* then it deactivates */
                
                break;
        }
    }
    s = active_num;
}

/* Density Functions */

void Density_Decay::time_loop(){
	do{
		advance_time();

		while( t > t_med[d] and d <= steps){
			rho_aux[d] = (double) s / L;
			++d;
		}

	}while(t < t_max and s > 0);
}

void Density_Decay::trial_loop(){
	int j = 0;
	while(j < N_trials){
		t = 0.;

		std::fill(cells.begin(), cells.end(), 1);

		for(int i = 0; i <  L; ++i)
			active_indexes.push_back(i);

		s = L;
		active_num = L;
		rho_aux[0] = 1.;

		d = 1;

		time_loop();

		++N_s;
		N_inv = (double) 1 / N_s;

		for(std::size_t n = 0; n < t_med.size(); ++n){
			aux_rho = N_inv * rho_aux[n] + (1 - N_inv) * rho_med[n];
			rho_std[n] = rho_std[n] + (rho_aux[n] - aux_rho) * (rho_aux[n] - rho_med[n]);
			rho_med[n] = aux_rho;
		} /* Mean value of rho at each selected time step */

		j++;

		std::fill(rho_aux.begin(), rho_aux.end(), 0.);
	}
	
	if(N_s > 1)
		for(std::size_t n = 0; n < t_med.size(); ++n)
			rho_std[n] = rho_std[n] / (N_s - 1);
}

vector<double> Density_Decay::vector_log_time(double t_max){
	double aux;
	int p_t = 0,
		d = 1,
		j = 0,
		steps;

	aux = t_max;
		
	while(aux > 1.){
		aux = aux / 10;
		p_t = p_t + 1;
	}
	steps = 1 + 9 * p_t;

	vector<double> time_vector(steps + 1, 0.);
	// time acquisition
	for(std::size_t n = 0; n < steps + 1; n++){

		time_vector[n] = j * d;
		
		++j;
		
		if (j == 10){
			d = d * 10;
			j = 1;
		}
	}

	return time_vector;
}