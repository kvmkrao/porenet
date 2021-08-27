
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>
#include <algorithm>

#include "MersenneTwister.h"

#include "network01.h"

using namespace std;

double Network::totThlen(const Pore& pore1, const Pore& pore2)
{
	double mx = pore2.p.xpos - pore1.p.xpos;
	double my = pore2.p.ypos - pore1.p.ypos;
	double mz = pore2.p.zpos - pore1.p.zpos;

	double totalLength = sqrt(mx * mx + my * my + mz * mz);
	return totalLength;
} 

void Network::netwsize(double xp, double yp, double zp, int np, int thn)
{
	 xsize = xp;
	 ysize = yp;
	 zsize = zp;
	 numpore = np;
	 throatNum = thn;
	 Data.resize(numpore);
	 Tdata.resize(throatNum);
	 inletOutletFactor = .022 * xsize;
	 double clayper(0.25),  minporerad(3.62E-6), maxporerad(7.35E-5); 
	 double diff = maxporerad - minporerad;

	cout << endl << "Input data read in from default data file." << endl
		 << "Dimension of equivalent network desired:	" << xsize << "m X " 
		 << ysize << "m X " << zsize << "m" << endl
		 << "Number of pores =				" << numpore << endl
		 << "Desired number of throats =			" << throatNum << endl;

	MTRand mtrand1;
	double PI = 3.14159;

	ifstream ins, throatin;
	ins.open ("PoreData.txt");
    if (ins.fail())
        {
                       cout << "Opening pore coordination number distrubtion data failed" << endl;
                       system ("pause");
                       exit (1);
        }

	for (int ww = 0; ww < numpore; ww++)
	{
		ins >> Data[ww].coordNo >> Data[ww].volume >> Data[ww].radius >> Data[ww].shapefac >> Data[ww].length; 
	}

	throatin.open("throatData.txt"); 
	if (throatin.fail())
			{
                       cout << "opening throat radius input file failed" << endl;
                       system ("pause");
                       exit (1);
			}
	
	int yh = throatNum; 	
	for (int tww = 0; tww < yh; tww++)
	{
		throatin >> Tdata[tww].thrRadius >> Tdata[tww].shapefac >> Tdata[tww].totlength >>Tdata[tww].length 
				 >> Tdata[tww].volume; 
	}

	sort(Tdata.begin(), Tdata.end(), compare());

	cout << endl << "Equivalent network generation in progress: " << endl;

	int i =0;
	cnumpore = numpore - 1; //indexing of pore coordination number

	for (int i = 0; i < numpore; i++)
		{
			n1.index = i + 1;
	 //randomly generate x, y, z position of each pore
			n1.p.xpos = mtrand1() * xp;
			n1.p.ypos = mtrand1() * yp;
			n1.p.zpos = mtrand1() * zp;
			n1.status = 0;
			n1.inletstatus = 0;
			n1.outletstatus = 0;
			n1.radius = Data[i].radius;
			n1.coorno = Data[i].coordNo;
			n1.currentCn = n1.coorno;
			n1.radiusvote = (n1.radius - minporerad) / diff;
			n1.shapefac = Data[i].shapefac; 
			n1.length =  Data[i].length; 
			n1.volume = Data[i].volume; 
			n1.clayvol =  ((n1.volume * clayper)/(1.0 - clayper)) ;
			pores.push_back(n1);
	}

	int sumthroat = 0;
	int suminlet = 0;
	int sumoutlet = 0;
	
	for (int j = 0; j < numpore; j++)
	{
		sumthroat = sumthroat + pores[j].coorno;

// determination of pore status as inlet
		if (pores[j].p.xpos < inletOutletFactor)
		{
			pores[j].status = 1;
			pores[j].inletstatus = 1;
			inletporesindex.push_back(pores[j].index);
			++suminlet;
		}

// determination of pore status as outlet
		if (pores[j].p.xpos > xsize - inletOutletFactor) 
		{
			pores[j].status = 2;
			pores[j].outletstatus = 1;
			outletporesindex.push_back(pores[j].index);
			++sumoutlet;
		}
	}


	int inletoutletsum = suminlet + sumoutlet; 
	
//Addition of inlet/outlet reservoir to the pores vector. Length and radius are randomly selected from pore data distribution
	n1.index = numpore + 1; 
	n1.length = 1.42705E-6;
	n1.radius = 1.92E-5;
	n1.radiusvote = (n1.radius -minporerad) / diff;
	pores.push_back(n1);
	n1.index = numpore +2; 
	n1.length = 1.42705E-6;
	n1.radius = 1.92E-5;
	n1.radiusvote = (n1.radius -minporerad) / diff;
	pores.push_back(n1);



//creation of inlet pores throat to inlet reservoir. Shape factor, length, total length and length are randomly selected 
// from throat data distribution
		for (int e = 0; e < suminlet; e++)
		{
			t1.index = e + 1;
			t1.n1.index = numpore + 1 ; 
			t1.n2.index = inletporesindex[e];
			pores[t1.n2.index - 1].currentCn = pores[inletporesindex[e] - 1].coorno - 1;
			pores[t1.n2.index - 1].porcoon.push_back(-1); 
			pores[t1.n2.index - 1].thrcoon.push_back(t1.index);
			t1.tradiusvote = pores[t1.n2.index - 1].radiusvote;
			sortedTradius.push_back(t1.tradiusvote);
			t1.shapefac = 1.23e-2;
			t1.volume = 0.0; 
			t1.clayvol = 0.0; 
			t1.length = 1.37E-5;
			t1.totlength = 1.16E-4; 
			t1.n1.length = 3.61E-5;
			t1.n2.length = pores[inletporesindex[e] - 1].length;
			 
			throats.push_back(t1);
			
		}

// creation of outlet pores throat to outlet reservoir.  Shape factor, length, total length and length are randomly selected 
// from throat data distribution
		for (int d = 0; d < sumoutlet; d++)
		{
			t1.index = suminlet + d + 1 ;
			t1.n1.index = outletporesindex[d];
			pores[outletporesindex[d] - 1].currentCn = pores[outletporesindex[d] - 1].coorno - 1;
			t1.n2.index = numpore + 2; 
			pores[t1.n1.index - 1].porcoon.push_back(0); 
			pores[t1.n1.index - 1].thrcoon.push_back(t1.index);
			t1.tradiusvote = pores[t1.n1.index - 1].radiusvote ;
			sortedTradius.push_back(t1.tradiusvote);
			t1.shapefac = 1.23e-2; 
			t1.volume = 0.0; 
			t1.clayvol = 0.0; 
			t1.length = 1.37E-5;
			t1.totlength = 1.16E-4; 
			t1.n1.length = pores[outletporesindex[d] - 1].length;
			t1.n2.length = 3.61E-5;
			
			throats.push_back(t1);
		}
		
	numthroat = (sumthroat - inletoutletsum) / 2 + inletoutletsum;
	double avgcoorno = double (sumthroat) / double (numpore);

	cout << endl 
		<<"Total number of inlet pores =			" << suminlet << endl
		<<"Total number of outlet pores =			" << sumoutlet << endl << endl;
	
//Determination of throat length based on the nearest pores to a pore of interest. 
	int ay = inletoutletsum + 1;
	int ay2 = ay;
	for (int re = 0; re < numpore; re++)
	{
		for (int ff = re + 1; ff < numpore; ff++)
		{
			double length = totThlen(pores[re], pores[ff]);
			if (length <= 0.15 * ysize)
			{
			pores[re].pore2Allocate.first = length;
			pores[re].pore2Allocate.second = pores[ff].index;
			allocation.push_back(pores[re].pore2Allocate);
			}
		}

		vector < pair <double, int> > ::iterator start = allocation.begin();          
		vector < pair <double, int> > ::iterator end = allocation.end();
	    sort (start, end), allocation[re].first; 
		
		int loopcounter(0);
		
		for ( size_t jk = 0; jk < allocation.size(); jk++)
		{			
			double overlap = allocation[jk].first - pores[re].radius - pores[allocation[jk].second - 1].radius;
			size_t uuu = allocation.size()-loopcounter;
			if ( pores[re].currentCn >= 1 && uuu >= 1)
				{
					if (pores[allocation[jk].second - 1].currentCn > 0 && overlap >= 0.0) 
						{
						t1.index = ay;			
						ay = ay + 1;
						t1.n1.index = pores[re].index;
						t1.n2.index = allocation[jk].second;
						pores[t1.n2.index - 1].currentCn = pores[t1.n2.index - 1].currentCn - 1;
						t1.tradiusvote = (pores[t1.n1.index - 1].radiusvote + pores[t1.n2.index - 1].radiusvote) / 2.0;
						sortedTradius.push_back(t1.tradiusvote);
						pores[re].porcoon.push_back(t1.n2.index); 
						pores[re].thrcoon.push_back(t1.index); 
						pores[t1.n2.index - 1].porcoon.push_back(t1.n1.index);
						pores[t1.n2.index - 1].thrcoon.push_back(t1.index);
						t1.n1.length = pores[t1.n1.index -1].length;
						t1.n2.length = pores[t1.n2.index -1].length;								
						throats.push_back(t1);
						loopcounter = loopcounter + 1;
						pores[re].currentCn = pores[re].currentCn - 1;

						}
				}
		}
		allocation.clear();
	}
	
	cout << endl <<"Total throat number connected =			" << throats.size() << endl
		<< "Percentage of unconnected throats =		" << static_cast<double>((numthroat - throats.size()) * 100)/ 
		static_cast<double>(numthroat) << endl;

// Determination of coordination number 
	for (int m = 0; m < numpore; m++)
	{
		int sumpore = pores[m].coorno - pores[m].currentCn;
		thrporcoorno.push_back(sumpore);
	}

		vector < double > ::iterator comm = sortedTradius.begin();          
		vector < double > ::iterator stopp = sortedTradius.end();
		sort (comm, stopp); 
		
		for (size_t j6 = 0; j6 < throats.size(); j6++)
		{
			size_t j7 = -1; 
			do
			{
				j7 = j7 + 1;
			}while ((throats[j6].tradiusvote != sortedTradius[j7]) && j7 < throats.size());

			if(Tdata[j7].thrRadius <= pores[throats[j6].n1.index - 1].radius && 
				Tdata[j7].thrRadius <= pores[throats[j6].n2.index - 1].radius )	
			{
			throats[j6].radius = Tdata[j7].thrRadius;
		
						sortedTradius[j7] = 0.0;
						throats[j6].shapefac = Tdata[j7].shapefac; 
						throats[j6].totlength = Tdata[j7].totlength;
						throats[j6].length = Tdata[j7].length;
						throats[j6].volume = Tdata[j7].volume;
						throats[j6].clayvol =  ((throats[j6].volume * clayper)/(1.0 - clayper)) ;			
			}
			else
			{
				do
				{
					--j7 ;

				}while((Tdata[j7].thrRadius > pores[throats[j6].n1.index - 1].radius || 
					Tdata[j7].thrRadius > pores[throats[6].n2.index - 1].radius ) && j7 > 0);
				throats[j6].radius = Tdata[j7].thrRadius;
						sortedTradius[j7] = 0.0;
					throats[j6].shapefac = Tdata[j7].shapefac; 
					throats[j6].totlength = Tdata[j7].totlength;
					throats[j6].length = Tdata[j7].length;
					throats[j6].volume = Tdata[j7].volume;
					throats[j6].clayvol =  ((throats[j6].volume * clayper)/(1.0 - clayper)) ;				
			}
		}



//	calculation of total volumes of pores and clay volume
	double sumporevol(0);
	double sumporeclayvol(0);
	for (int ad = 0; ad < numpore; ad++) 
	{
		sumporevol = sumporevol + pores[ad].volume;
		sumporeclayvol = sumporeclayvol + pores[ad].clayvol;
	} 

	//determination of effective throat radius
	for (int t = 0; t < numpore; t++)
	{
		double sumthroatradius(0), maxradius(0);
		for (size_t z = 0; z < pores[t].thrcoon.size(); z++)
		{
			sumthroatradius = sumthroatradius + throats[pores[t].thrcoon[z] - 1].radius;
			maxradius = max (maxradius, throats[pores[t].thrcoon[z] - 1].radius);
		}
		pores[t].effthradius = sumthroatradius / double(pores[t].coorno);
		pores[t].aspectRat = pores[t].radius / pores[t].effthradius;
	}

	// determination of total throat volume and clay volume
		double sumthroatvol(0);
		double sumthroatclayvol(0);
		for (size_t as = 0; as < throats.size(); as++)
		{
			sumthroatvol = sumthroatvol + throats[as].volume; 
			sumthroatclayvol = sumthroatclayvol + throats[as].clayvol;
		}

	//Changing of inlet and outlet pore indexes to -1 and 0 to conform with 2-phase code network input data format;
		for (size_t j13= 0; j13< throats.size(); j13++)
		{
			if (throats[j13].n1.index == numpore + 1) throats[j13].n1.index = -1;
			if (throats[j13].n2.index == numpore + 2) throats[j13].n2.index = 0;
		}

		cout << endl << "Writing network data to output files." << endl;

		ofstream outn1, outn2, outt1, outt2;
	
	if (outn1.fail())
        {
                       cout << "opening pore 1 output file failed" << endl;
                       system ("pause");
                       exit (1);
		}

	outn1.open ("EqBerea_node1.dat");
	outn1.setf(ios::showpoint);
	outn1.setf(ios::scientific);
	outn1.precision(4);
	outn1 << setw(20) << numpore << setw(15) << xsize << setw(15) << ysize << setw(15) << zsize << endl;
	
		for (int jjq = 0; jjq < numpore; jjq++)
	{
			outn1 <<setw(20)<< pores[jjq].index << setw(20) <<pores[jjq].p.xpos <<  setw(20)<< pores[jjq].p.ypos
				  <<  setw(20) << pores[jjq].p.zpos <<  setw(20) << thrporcoorno[jjq]
				  <<  setw(20);
			for (size_t h = 0; h < pores[jjq].porcoon.size(); h++) outn1 << pores[jjq].porcoon[h] <<  setw(15) ;
			outn1 << pores[jjq].inletstatus <<  setw(15) << pores[jjq].outletstatus <<  setw(15);
			 for (size_t w = 0; w < pores[jjq].thrcoon.size(); w++) outn1 << pores[jjq].thrcoon[w] <<  setw(15);
			 outn1 << endl;
	}

	outn2.open ("EqBerea_node2.dat");
		outn2.setf(ios::showpoint);
		outn2.setf(ios::scientific);
		outn2.precision(4);
	
	for (int jc = 0; jc < numpore; jc++)
	{
		outn2 << pores[jc].index << setw(20) << pores[jc].volume << setw(20) << pores[jc].radius 
			<< setw(20)<< pores[jc].shapefac<< setw(20) << pores[jc].clayvol <<endl;
	}

	outt1.open ("EqBerea_link1.dat");
	outt1.setf(ios::showpoint);
	outt1.setf(ios::scientific);
	outt1.precision(4);
	outt1 << throats.size() << endl;
	
	for (size_t qs = 0; qs < throats.size(); qs++)
		{
			outt1 <<  throats[qs].index << setw(15) <<throats[qs].n1.index << setw(15)<< throats[qs].n2.index
				<< setw(15) << throats[qs].radius << setw(15) << throats[qs].shapefac	<< setw(15) 
				<< throats[qs].totlength << endl; 
		}

		outt2.open ("EqBerea_link2.dat");
		outt2.setf(ios::showpoint);
		outt2.setf(ios::scientific);
		outt2.precision(4);
	
		for (size_t qa = 0; qa < throats.size(); qa++)
		{
			outt2 <<  throats[qa].index << setw(15)<<throats[qa].n1.index << setw(15) << throats[qa].n2.index
				<< setw(15)<< throats[qa].n1.length << setw(15) << throats[qa].n2.length <<setw(15)
				<< throats[qa].length << setw(15) << throats[qa].volume <<setw(15) << throats[qa].clayvol
				<< endl; 
		}

		cout << endl << "Equivalent network generation succesfully completed." << endl;
}




