#ifndef NETWORK01_H
#define NETWORK01_H

using namespace std;

struct Point
{
	double xpos, ypos, zpos;
};

struct Pore
{
	int index;
	Point p;
	double radius;
	double length;
	double volume;
	int coorno;
	int currentCn;
	double clayvol;
	double shapefac;
	int status;
	int inletstatus;
	int outletstatus;
	vector<int> thrcoon;
	vector<int> porcoon;
	vector<int> porcoon2;
	double aspectRat;
	double radiusvote;
	double effthradius;
	pair<double, int> pore2Allocate;
	
 };

struct Throat
{
	int index;
	Pore n1, n2;
	double radius;
	double totlength;
	double length;
	double volume;
	double clayvol;
	double shapefac;
	double tradiusvote;
};


struct poreData
{
	int coordNo;
	double length;
	double volume;
	double radius;
	double shapefac;
};

class throatData
{
public:
	double thrRadius;
	double totlength;
	double length;
	double volume;
	double shapefac;
};

class compare
{
public:
	bool operator()(const throatData& pt1, const throatData& pt2) const {return pt1.thrRadius < pt2.thrRadius;}
};
class Network
{	
    public:
		Network() {}
		~Network(){}

	/*	Network()
		{
			Data = new poreData[numpore];
			Tdata = new throatData[throatNum]; 
		}
		~Network()
		{
			delete [] Data;
			delete [] Tdata; 
		} */

		double xsize;
		double ysize;
		double zsize;
		double inletOutletFactor;
		int numpore, cnumpore, throatNum;
		int numthroat;
		void netwsize(double, double, double, int, int);
		Pore n1, n2;
		Throat t1, t2;

		vector<poreData> Data;
		vector<throatData> Tdata; 
		double weibull(double, double, double, double) const;
		double totThlen(const Pore&, const Pore&);
		double thlen(const Throat& t1);
		double thvol(const Throat& t1);
		vector<Pore> pores;
		vector<Throat> throats;
		vector<int> thrpor;
		vector<int> thrporcomp;
		vector<int> thrporcoorno;
		vector <double> shapefactor;
		vector <int> inletporesindex;
		vector <int> outletporesindex;
		vector < pair <double, int> > allocation;
		vector<double> sortedTradius;
	//	vector<double> sortedTr;
		vector<double> lengthc;
		vector<double> poresfac;
		vector<double> throatsfac;
		double max (double a, double b) {return a > b ? a : b;}

 };

#endif