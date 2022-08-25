# Contact Process Class
Creates an object that holds the simulation itself
1. The first parameter determines the size of the network
	- if HMN1 it is the s in 2^s^
	- otherwise it is the amount of sites.
3. The second parameter determines the control parameter
	- if model A lambda is any number
	- if model B lambda is between 0 and 1
4. The third parameter defines the model
5. And, the forth parameter defines the network
	- if 0 it's a random network with a fixec number of neighbors
	- if 1 it's a fully connected network
	- if 2 it opens the human connectome file
	- if 3 it's a HMN-1 network
	- and if 4 it opens the ER network file
