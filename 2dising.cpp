
//
// ---------------------------------------------------------------------------------
// 
// This program simulates a 2D Ising magnet using single- and cluster moves.
// 
// Things that should be fixed:
//		
//		x Ambiguous ordering of neighbouring nodes
//		x Using vectors or single value types for neighbours
// 
// ---------------------------------------------------------------------------------
//

#include <iostream>
#include <fstream>	
#include <vector>
#include <random>
#include <cmath>
#include <sstream>
#include <iomanip>

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
		double entop = 0;
		double enlft = 0;
		double enrgt = 0;
		double enbot = 0;
		node* top;
		node* lft;
		node* rgt;
		node* bot;
		vector <node*> nears;	// a list of references to all neighbours, declared in findnear();
		bool free;

		void nodeenergy() {energy = entop+enlft+enrgt+enbot;}

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
	double M = 0;
	double J;								// interaction energy
	double T;								// temperature
	
	vector <node> nodes;					// vector of all nodes
	vector <double> corrvec;				// for correlations

	// public section
	public:
	
	// function members
	void makegrid();					// set up grid of nodes
	void findnear();					// find all neighbours
	void magnetize(mt19937&,float);
	void magnetizemin();
	int metropolis(mt19937&);
	double pairenergy(node*,node*);
	void calcenergy();
	void wolffcluster(mt19937&);
	void addnode(node*,vector<node*>*,mt19937&,bernoulli_distribution&,int);
	void updatenode(node*,double,double,double,double);
	void correlation(int);
	void calcmagnetization();
	void zeroenergies();
	void nodedist();
	int getspin(int);
	int getmag();
	void setcorrvec(int);
	void printpos(int);					// accessor for getting position
	void printnear(int);
	void printenergy();
	void printenergynode(int);
	void writeenergy(ofstream&);
	void writenodes(ofstream&);
	void writecluster(vector<node*>*,ofstream&);
	void writeimage(int,int);
	void writecorr(int);
	void writemag(ofstream&);

	// constructor
	magnet(int ixsize, int iysize, double iJ, double iT)
		: xsize(ixsize), ysize(iysize), J(iJ), T(iT)
	{
		size = xsize*ysize;
		nodes.reserve(size);					// vector of all nodes	
	}
};


//-------------------------------------------------------------
// out of line declarations for node class
//-------------------------------------------------------------


void magnet::node::setpos(int x, int y, int i)	// accessor for setting position
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

void magnet::printpos(int n)						// accessor for getting position
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
	ostringstream namestream;
	namestream << "image" << setfill('0') << setw(m) << n << ".pbm";
	string filename = namestream.str();

	ofstream imgstream;
	imgstream.open(filename);

	imgstream << "P1" << endl << xsize << " " << ysize << endl;

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
		//el.nears = {el.rgt,el.bot};
	}
}

void magnet::magnetize(mt19937& gen, float bias)
{
	//random_device rdd;
	//mt19937 genn(rdd());
	
	bernoulli_distribution randspin(bias);
	for (auto& el : nodes)
	{
		int spin = randspin(gen)*2-1;
		el.spin = spin;
	}
}

void magnet::magnetizemin()
{
	//random_device rdd;
	//mt19937 genn(rdd());
	
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
	for (auto& el : nodes)
	{
		double enrgt = pairenergy(&el,el.rgt);
		double enbot = pairenergy(&el,el.bot);

		//el.energy += enrgt;
		//el.energy += enbot;
		//el.rgt->energy += enrgt;
		//el.bot->energy += enbot;
		
		el.enrgt = enrgt;
		el.enbot = enbot;
		//el.rgt->enlft = enrgt;
		//el.bot->entop = enbot;
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

	//updatenode(rand,0,enbot,0,enrgt);
	updatenode(rand,entop,enbot,enlft,enrgt);

	double newnodeen = (entop + enbot + enlft + enrgt);
	//double newnodeen = (enrgt + enbot);
	double nodediff = (newnodeen - oldnodeen );

	energy = olden + nodediff;
	//calcenergy();

	//double endiff = energy - olden;
		
	double prob = 1.0 - exp(-nodediff/(kb*T));
	bernoulli_distribution acc(prob);

	//cout << "endiff " << endiff << endl;
	//cout << "prob: " << prob << endl;
	
	if (nodediff < 0)
		return 1;
		//{cout << "accepted " << prob << endl; return 1;}
	else if (acc(gen))
		return 1;
		//{cout << "metrepted " << prob << endl; return 1;}
	else
	{
		energy = olden;
		rand->spin = oldspin;
		updatenode(rand,oldentop,oldenbot,oldenlft,oldenrgt);
		return 0;
	}
		//{cout << "rejected " << prob << endl; energy = olden; rand->spin = oldspin; return 0;}
}


void magnet::addnode(node* anode, vector<node*>* cluster, mt19937& gen, bernoulli_distribution& dist, int spin)
{
	//bool p = dist(gen);
	
	//cout << "free? " << anode->free << endl << "spin: " << anode->spin << " " << spin << endl << "prob " << p << endl;
	
	//if (!(anode->free) && anode->spin == spin && p)
	
	if (anode->free && anode->spin == spin && dist(gen))
	{
		//gbug++;
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
	for (auto& el : nodes)					// reset cluster
		el.free = true;
	
	uniform_int_distribution<> randnode (0,nodes.size()-1);		// dist for chosing random node

	unsigned int n = randnode(gen);			// a random node
	node* rand = &nodes[n];					// increase readability
	int spin = rand->spin;					// store spin

	//cout << "randomly picked node " << n << " with spin " << spin << endl;

	vector <node*> cluster;					// vector to hold cluster
	cluster.reserve(static_cast<unsigned int>(sqrt(nodes.size())));	// estimate size

	double accprob = 1 - exp(-2*J/(kb*T));
	//double accprob = 0.4;

	//double tmp = pairenergy(&nodes[n],&nodes[n+10]);
	//cout << nodes[n].spin << " " << nodes[n+10].spin << " " << tmp <<  endl;
	
	//double accprob = 0.9;
	//cout << "accprob " << accprob << endl;

	bernoulli_distribution  dist (accprob);
	
	//gbug = 0;
	addnode(rand,&cluster,gen,dist,spin);
	
	ofstream fs;
	fs.open("cluster.csv");
	writecluster(&cluster,fs);
	fs.close();

	for (auto& el : cluster)
		el->spin = -spin;
	

	//nodedist();
	//printenergy();

	//for (auto& el : cluster)
	//	printpos(el->index);
	
	//cout << "in cluster: " << gbug << endl;
	//gbug = 0;
}

void magnet::correlation(int n)
{
	for (unsigned int k = 0; k < corrvec.size(); k++)
	{
		corrvec[k] += nodes[n].spin * nodes[n+k].spin;
	}
}

int main()
//int main(int argc, char** argv)
{	
	//argc++; argc--;					// ugly hack to placiate clang...
	//int xsize = atoi(argv[1]);		// horezontal atoms
	//int ysize = atoi(argv[1]);		// horezontal atoms
	//int corrindex = atoi(argv[2]);

	//corrindex < (double)xsize/2 ? cout << "calculating correlation at " << corrindex << endl
	//	: cout << "warning: correlation check longer than width " << endl;
	
	int xsize = 512;		// horezontal atoms
	int ysize = 512;		// vertical atoms
	double T = 300;			// temperature
	//double J = kb*T;		// interaction energy of the ising model
	double J = 0.0175;		// interaction energy of the ising model		// 0.0175
	float bias = 0.5;		// statistical bias in the spin distribution
	int steps = 1;
	//int count = 0;
	int corrindex = 0;
	int digits = to_string(steps+1).length();

	magnet ising(xsize,ysize,J,T);	// initialize magnet

	ising.setcorrvec(floor(xsize/2.0));

	random_device rd;		// using unix random device as seed
	mt19937 gen(rd());		// Mersenne Twister as RNG

	ising.makegrid();		// 
	ising.findnear();
	ising.magnetize(gen,bias);
	//ising.magnetizemin();

	ising.zeroenergies();
	ising.calcenergy();
	ising.calcmagnetization();

	//cout << "energy before " << endl;
	//ising.printenergy();
	//cout << "initial nodedist " << endl;
	//ising.nodedist();
	
	ising.writeimage(0,digits);

	//ofstream nodestream;
	//nodestream.open("nodes.csv");
	//ising.writenodes(nodestream);
	//nodestream.close();

	ofstream energystream;
	energystream.open("energy.dat");
	ising.writeenergy(energystream);
	
	ofstream magstream;
	magstream.open("magnetization.dat");
	ising.writemag(magstream);

	//for (int k = 0; k < steps; k++)
	//{
		//count += ising.metropolis(gen);
	//	ising.metropolis(gen);
		
	//	if(!(k % 100))
	//	{
	//		ising.writeenergy(energystream);
	//		ising.writeimage(k+1,digits);
	//		ising.correlation(corrindex);
	//		ising.calcmagnetization();
	//		ising.writemag(magstream);
	//	}
	//}

	//double accprob = (float)count / (float)steps;
	//cout << "acceptance prob: " << accprob << endl;

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

	energystream.close();
	magstream.close();

	ising.writecorr(steps);

	//cout << "energy after " << endl;
	//ising.printenergy();
	cout << "final nodedist " << endl;
	ising.nodedist();
	
	ofstream nodestream;
	nodestream.open("nodes.csv");
	ising.writenodes(nodestream);
	nodestream.close();
}



























