
//
// ---------------------------------------------------------------------------------
// 
// This program simulates a 2D Ising magnet using single- and cluster moves.
// 
// Things that should be fixed or coniddered:
//		
//		x Ambiguous ordering of neighbouring nodes
//		x Using vectors or single value types for neighbours
//		x Energy reconsideration rather than recalculation for cluster moves
//		x The boolean free variable could be eliminate with immiate flips
//
// ---------------------------------------------------------------------------------
//

#include <iostream>				// for basic inout output
#include <fstream>				// for read/write to files
#include <vector>				// for vector template functionallity
#include <random>				// for random distributions and generators
#include <cmath>				// for mathematical functions
#include <sstream>				// for creating namestrings
#include <iomanip>				// for formating output strings

using namespace std;

//double kb = 1.38064852e-23;	// Boltzmann's constant in SI-units
double kb = 8.6173303e-5;		// Boltzmann's constant in eV
int pingcount = 0;				// for debugging
double gbug = 0;

void ping() {pingcount++; cout << "ping" << pingcount << endl;}	// Simple debugging

class magnet					// class for the magnet
{	
	// private section
	
	// object members
	class node					// class for one spin-element
	{
		// private data members
		int index;				// number
		int pox;				// x position
		int poy;				// y position
		int spin;				// spin number
		double energy;			// energy of node
		double entop = 0;		// energy to top neighbour
		double enlft = 0;		// energy to left neighbour
		double enrgt = 0;		// energy to right neighbour
		double enbot = 0;		// energy to bottom neighbour
		node* top;				// pointer to top neigh
		node* lft;				// pointer to left neigh
		node* rgt;				// pointer to right neigh
		node* bot;				// pointer to bot neigh
		vector <node*> nears;	// a list of references to all neighbours
		bool free;				// was node already visited by algorithm

		void nodeenergy() {energy = entop+enlft+enrgt+enbot;} // calcs node energy

		// public section
		public:
		// accessors
		void setpos(int,int,int);			// for setting position
		void setspin(int s){spin = s;}		// for setting spin
		int getspin() {return spin;}		// for getting spin

		friend class magnet;				// granting magnet friendship privileges
	};
	
	//data members
	int xsize;								// number of nodes in x direction
	int ysize;								// number of nodes in y direction
	int size;								// total size
	double energy;							// energy of magnet
	double M = 0;							// mangentization
	double J;								// interaction energy
	double T;								// temperature
	
	vector <node> nodes;					// vector of all nodes
	vector <double> corrvec;				// for correlations

	// public section
	public:
	
	// function members
	void makegrid();					// set up grid of nodes
	void findnear();					// find all neighbours
	void magnetize(mt19937&,float);		// set magnetization
	void magnetizemax();				// set maximum energy magnetization
	int metropolis(mt19937&);			// metropolis algorithm
	double pairenergy(node*,node*);		// returns energy between two nodes
	void calcenergy();					// calculates total energy of magnet
	void wolffcluster(mt19937&);		// wolff's algorithm
	void addnode(node*,vector<node*>*,
		 mt19937&,bernoulli_distribution&,int);	// add node to cluster
	void updatenode(node*,double,double,double,double);	// update energy on node
	void correlation(int);				// calculate correlation 
	void calcmagnetization();			// calculate magnetization
	void zeroenergies();				// set all energies to zero
	void nodedist();					// prond node distribution
	int getspin(int);					// get spin of node
	int getmag();						// returns magnetization
	void setcorrvec(int);				// build correlation vector
	void printpos(int);					// accessor for getting position
	void printnear(int);				// accessor for getting neigh
	void printenergy();					// prints energy
	void printenergynode(int);			// prints node energy
	void writeenergy(ofstream&);		// writes energy to file
	void writenodes(ofstream&);			// writes node distribution to file
	void writecluster(vector<node*>*,ofstream&);	// writes cluster to file
	void writeimage(int,int);			// write image of node distribution
	void writecorr(int);				// write correlations to file
	void writemag(ofstream&);			// write magnetization to file

	// constructor
	magnet(int ixsize, int iysize, double iJ, double iT)
		: xsize(ixsize), ysize(iysize), J(iJ), T(iT)
	{
		size = xsize*ysize;						// calculate number of nodes
		nodes.reserve(size);					// set size of node vector	
	}
};


//-------------------------------------------------------------
// out of line declarations for node class
//-------------------------------------------------------------


void magnet::node::setpos(int x, int y, int i)
{
	pox = x; poy = y; index = i;
}


//-------------------------------------------------------------
// out of line declarations of magnet class
// ------------------------------------------------------------

int magnet::getspin(int n)
{
	return nodes[n].spin;
}

int magnet::getmag()
{
	return M;
}

void magnet::setcorrvec(int n)
{
	corrvec = vector <double> (n);
}

void magnet::printnear(int n)
{
	node anode = nodes[n];
	cout << "top: " << anode.top->index << " left: " << anode.lft->index << " right: " << anode.rgt->index << " bot: " << anode.bot->index << endl;
}

void magnet::printpos(int n)					
{	
	node anode = nodes[n];
	cout << "index: " << anode.index << endl <<  "position: (" << anode.pox << ","	<< anode.poy << ")" << endl;
}

void magnet::printenergynode(int n)
{
	node anode = nodes[n];
	cout << "index " << anode.index << endl << "energy " << anode.energy << endl << "spin: " << anode.spin << endl;
}

void magnet::printenergy()
{
	cout << "total energy: " << energy << endl;
}

void magnet::writeenergy(ofstream& fs)
{
	fs << energy << endl;
}

void magnet::writenodes(ofstream& fs)
{
	for (auto& el : nodes)
		fs <<  el.pox << "," << el.poy << "," << el.spin << endl;
}

void magnet::writecluster(vector <node*>* cluster, ofstream& fs)
{
	for (auto& el : (*cluster))
		fs <<  (*el).pox << "," << (*el).poy << "," << (*el).spin << endl;
}

void magnet::writecorr(int steps)
{
	ofstream corrstream;
	corrstream.open("correlation.dat");
	
	for (auto& el : corrvec)
		corrstream << el/steps << endl;
	
	corrstream.close();
}

void magnet::writemag(ofstream& fs)
{
	fs << M << endl;
}

void magnet::writeimage(int n,int m)
{
	// create file name of fixed size
	ostringstream namestream;
	namestream << "image" << setfill('0') << setw(m) << n << ".pbm";
	string filename = namestream.str();

	ofstream imgstream;
	imgstream.open(filename);

	// set pbm header
	imgstream << "P1" << endl << xsize << " " << ysize << endl;	

	// build pbm image
	int index = 0;
	for(int irow = 0; irow < ysize; irow++)	
		for (int icol = 0; icol < xsize; icol++)
		{
			index = irow*xsize + icol;
			
			int spin = nodes[index].spin;
			if (spin == 1 )
				imgstream << 1 << " ";
			else 
				imgstream << 0 << " ";

			if (icol == xsize - 1)
				imgstream << endl;
		}

}

void magnet::makegrid()
{
	for (int x = 0; x < xsize; x++)
		for (int y = 0; y < ysize; y++)
		{
			node anode;
			int i = x*xsize + y;
			anode.setpos(x,y,i);
			nodes.push_back(anode);
		}
}

void magnet::findnear()
{	
	// to increase readability
	int rows = ysize;
	int cols = xsize;

	for (int x = 0; x < rows; x++)
		for (int y = 0; y < cols; y++)
		{   
			int i = x*cols + y;

		   	//periodic boundaries
          	int periodicup      = 0;  
 			int periodicdown    = 0;
			int periodicright   = 0;
			int periodicleft    = 0;

			if (x == 0)
				periodicup = rows*(cols);
			if (x == rows - 1)
				periodicdown = -rows*cols;
			if (y == 0)
				periodicleft = cols;
			if (y == cols - 1)
				periodicright = -cols;

			// for debugging
			//cout << "node  " << i << endl;
			//cout << "left  " << i-1 + periodicleft << endl;
			//cout << "right " << i+1 + periodicright << endl;
			//cout << "up    " << i-cols + periodicup << endl;
			//cout << "down  " << i+cols + periodicdown << endl << endl;
		
			nodes[i].lft = (&nodes[i-1 + periodicleft]);
			nodes[i].rgt = (&nodes[i+1 + periodicright]);
			nodes[i].top = (&nodes[i-cols + periodicup]);
			nodes[i].bot = (&nodes[i+cols + periodicdown]);
		}

	for (auto& el : nodes)
	{
		el.nears = {el.top,el.lft,el.rgt,el.bot};
	}
}

void magnet::magnetize(mt19937& gen, float bias)
{
	bernoulli_distribution randspin(bias);
	for (auto& el : nodes)
	{
		int spin = randspin(gen)*2-1;
		el.spin = spin;
	}
}

void magnet::magnetizemax()
{
	for (auto& el : nodes)
	{
		el.index % 2 ? el.spin = 1 : el.spin = -1;
	}
}

void magnet::calcmagnetization()
{
	M = 0;

	for (auto& el : nodes)
		M += el.spin;
	
	M /= ((double) nodes.size());
}

void magnet::nodedist()
{
	int ups = 0;
	int downs = 0;
	for (auto& el : nodes)
	{
		if (el.spin == 1)
			ups++;
		else
			downs++;
	}

	cout << "number of spin up " << ups << endl << "number of spin down " << downs << endl;
}

void magnet::zeroenergies()
{
	for (auto& el : nodes)
	{
		el.energy = 0;
		el.entop = 0;
		el.enlft = 0;
		el.enrgt = 0;
		el.enbot = 0;
	}

}

double magnet::pairenergy(node* node1, node* node2)
{
	double nodeen = 0;
	
	nodeen = - J * node1->spin * node2->spin;
	
	//cout << "node energy " << nodeen << endl << "spin 1 " << node1->spin << "spin 2 " << node2-> spin << endl;

	return nodeen;
}

void magnet::calcenergy()
{
	// calculate energy with forward neighbours
	for (auto& el : nodes)
	{
		double enrgt = pairenergy(&el,el.rgt);
		double enbot = pairenergy(&el,el.bot);

		el.enrgt = enrgt;
		el.enbot = enbot;
	}

	for (auto& el : nodes)
		el.nodeenergy();

	double accen = 0;
	
	for (auto& el : nodes)
		accen += el.energy;

	energy = accen;
}

	
void magnet::updatenode(node* anode, double entop, double enbot, double enlft, double enrgt)
{
	anode->entop = entop; 
	anode->enbot = enbot; 
	anode->enlft = enlft; 
	anode->enrgt = enrgt; 
	anode->top->entop = enbot;
	anode->bot->enbot = entop;
	anode->rgt->enlft = enrgt;
	anode->lft->enrgt = enlft;
	
	anode->nodeenergy();
	anode->top->nodeenergy();
	anode->bot->nodeenergy();
	anode->lft->nodeenergy();
	anode->rgt->nodeenergy();
}

int magnet::metropolis(mt19937& gen)
{
	uniform_int_distribution<> randnode (0,nodes.size()-1);
	uniform_int_distribution<> randspin(0,1);
	
	int n = randnode(gen);
	int spin = randspin(gen)*2-1;
	
	node* rand = &nodes[n];
	node* top = rand->top;
	node* bot = rand->bot;
	node* lft = rand->lft;
	node* rgt = rand->rgt;

	int oldspin = rand->spin;
	double olden = energy;
	double oldnodeen = rand->energy;
	double oldentop = rand->entop;
	double oldenbot = rand->enbot;
	double oldenlft = rand->enlft;
	double oldenrgt = rand->enrgt;

	rand->spin = spin;

	double entop = pairenergy(rand,top);
	double enbot = pairenergy(rand,bot);
	double enlft = pairenergy(rand,lft);
	double enrgt = pairenergy(rand,rgt);

	updatenode(rand,entop,enbot,enlft,enrgt);

	double newnodeen = (entop + enbot + enlft + enrgt);
	double nodediff = (newnodeen - oldnodeen );

	energy = olden + nodediff;

	double prob = 1.0 - exp(-nodediff/(kb*T));
	bernoulli_distribution acc(prob);

	if (nodediff < 0)
		return 1;
	else if (acc(gen))
		return 1;
	else
	{
		energy = olden;
		rand->spin = oldspin;
		updatenode(rand,oldentop,oldenbot,oldenlft,oldenrgt);
		return 0;
	}
}


void magnet::addnode(node* anode, vector<node*>* cluster, mt19937& gen, bernoulli_distribution& dist, int spin)
{
	if (anode->free && anode->spin == spin && dist(gen))
	{
		anode->free = false;
		(*cluster).push_back(anode);
		
		for (auto& el : anode->nears)
			addnode(el,cluster,gen,dist,spin);
	}
	else
		anode->free = false;
}

void magnet::wolffcluster(mt19937& gen)
{
	// reset cluster state
	for (auto& el : nodes)
		el.free = true;
	
	uniform_int_distribution<> randnode (0,nodes.size()-1);	

	unsigned int n = randnode(gen);		
	node* rand = &nodes[n];					
	int spin = rand->spin;					

	// construct temporary cluster vector
	vector <node*> cluster;					
	cluster.reserve(static_cast<unsigned int>(sqrt(nodes.size())));
	
	double accprob = 1 - exp(-2*J/(kb*T));

	bernoulli_distribution  dist (accprob);
	
	addnode(rand,&cluster,gen,dist,spin);
	
	// uncomment to write cluster to file
	//ofstream fs;
	//fs.open("cluster.csv");
	//writecluster(&cluster,fs);
	//fs.close();

	for (auto& el : cluster)
		el->spin = -spin;
}

void magnet::correlation(int n)
{
	for (unsigned int k = 0; k < corrvec.size(); k++)
	{
		corrvec[k] += nodes[n].spin * nodes[n+k].spin;
	}
}

int main()
{	
	int xsize = 512;				// horezontal atoms
	int ysize = 512;				// vertical atoms
	double T = 300;					// temperature
	double J = 0.0175;				// interaction energy of the ising model	
	float bias = 0.5;				// statistical bias in the spin distribution
	int steps = 10000;				// timesteps
	int count = 0;					// for acceptance probability
	int corrindex = 0;				// usd by correlation function
	int digits = to_string(steps+1).length();	// for file name
	bool usemetro = false;			// for using metropolis
	bool usewolff = true;			// for using wolff's

	magnet ising(xsize,ysize,J,T);		// initialize magnet

	ising.setcorrvec(floor(xsize/2.0));	// create correlation vector

	random_device rd;				// using unix random device for seeding
	mt19937 gen(rd());				// Mersenne Twister as RNG

	ising.makegrid();				// set up grid
	ising.findnear();				// identify neibghbours
	ising.magnetize(gen,bias);		// magnetize with bias
	//ising.magnetizemax();			// maximize energy

	ising.zeroenergies();			// initialize energies as zero
	ising.calcenergy();				// calculate initial enrgies
	ising.calcmagnetization();		// calculate initial magnetization

	ising.writeimage(0,digits);		// write initial magnetization image

	ofstream energystream;			
	energystream.open("energy.dat");
	ising.writeenergy(energystream);	// write initial energies
	
	ofstream magstream;
	magstream.open("magnetization.dat");
	ising.writemag(magstream);			// write initial magnetization


	if (usemetro)						// run metropolis algorithm
	{
		for (int k = 0; k < steps; k++)
		{
			count += ising.metropolis(gen);
			ising.metropolis(gen);
		  
			if(!(k % 100))
			{
				ising.writeenergy(energystream);
				ising.writeimage(k+1,digits);
				ising.correlation(corrindex);
				ising.calcmagnetization();
				ising.writemag(magstream);
			}
		}

		double accprob = (float)count / (float)steps;
		cout << "acceptance prob: " << accprob << endl;
	}

	if (usewolff)					// run wolff's algorithm
	{
		for (int k = 0; k < steps; k++)
		{
			ising.wolffcluster(gen);
			ising.calcenergy();
			ising.writeenergy(energystream);
			ising.writeimage(k+1,digits);
			ising.correlation(corrindex);
			ising.calcmagnetization();
			ising.writemag(magstream);
		}
	}

	//close streams
	energystream.close();
	magstream.close();

	ising.writecorr(steps);

	// for diagnostic purposes
	cout << "final nodedist " << endl;
	ising.nodedist();
	
	//write final node distribution
	ofstream nodestream;
	nodestream.open("nodes.csv");
	ising.writenodes(nodestream);
	nodestream.close();
}



























