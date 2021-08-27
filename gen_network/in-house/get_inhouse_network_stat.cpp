/*
 * Author: V Kotteda
 * Date  : July 29, 2021 
 */

#include<cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <sstream>
#include <vector>
#include <iomanip>      // std::setiosflags
#include <stdlib.h>     /* srand, rand */

using namespace std;
const double PI  =3.141592653589793238463;
double nonle = 1e-9;
double m2nm  = 1.0/nonle;
double pnond = 1e6; 
const int maxrand = 50000;
double dt     = 5e-9 ;
int nsteps    = 1000; //1e-7/dt;  // number of time steps 
double nondmass = nonle*nonle*nonle*pnond ; // non-dimensional mass

string fnode ="/home1/vkotteda/muskat/pore_networks/toy_networks/toy_networks/murali_networks/10x10x10_uniform/node.dat"; 
string flink ="/home1/vkotteda/muskat/pore_networks/toy_networks/toy_networks/murali_networks/10x10x10_uniform/link.dat"; 
/***********************************************
Table 4.2: Berea network statistics. (Dr. Piri thesis) 

Item   Throats Pores Total
Number 26146  12349  38495
//Porosity excl. clay (%) 4.562 13.746 18.309
//Porosity incl. clay (%) 6.238 17.785 24.024
Average shape factor 0.035 0.033 0.034
Triangular cross-sections (%) 90.729 95.506 92.261
Square cross-sections (%) 7.542 4.324 6.510
Circular cross-sections (%) 1.729 0.170 1.229
//ASCHA (deg.)(1) 15.235 13.744 14.751
//AWCHA (deg.)(2) 48.828 49.215 48.954
Minimum radius („m) 0.903 3.623 0.903
Maximum radius („m) 56.850 73.539 73.539
Average radius („m) 10.970 19.167 13.60
Connected to the inlet 254 0 254
Connected to the outlet 267 0 267
//Isolated clusters - - 3
Isolated 3 6 9
Minimum coordination number - 1 -
Maximum coordination number - 19 -
Average coordination number - 4.192 -
//Clay volume (%) 1.676 4.039 5.715
//Kabs (cal. box: 0.05-0.95) (mD) - - 3055
X dimension (mm) - - 3
Y dimension (mm) - - 3
Z dimension (mm) - - 3
***********************************************/

void readnode(double* dia, std::vector<int>& ni, double* shfrn, double* voln) {
	size_t nodes, id;
	size_t position, shapf; 
	double nr, vol, shapeFactor, clay; 
	double dd1, dd2, dd3,dd4; 
	fstream nodeloc;
	nodeloc.open(fnode.c_str()); 
        nodeloc >> nodes >>  dd1 >>  dd2 >> dd3 >> dd4;
        cout.precision(3);
        string str1, str2, str3, str4, str5, str6;
	for(int i=0; i<nodes ; i++) {  
		// units are in nano meters
		nodeloc >> id >> str1 >> str2 >> str3 >> position >> ni[i] >> str4 >> str5 >> str6 >> clay;
                //xn[i]  = std::stod(str1);  
                //yn[i]  = std::stod(str2);  
                //zn[i]  = std::stod(str3);  
                dia[i]   = std::stod(str4); 
		voln[i]  = std::stod(str5);
		shfrn[i] = std::stod(str6);
        }  
	nodeloc.close();
	return;
}

void readlink(int* l1, int* l2, double* ld, double* ll, int* lid, double *volth, double *shfnl) {
	int links,lp,id;
        double vol, shapeFactor, clay, n1l, n2l, distance; 
        
	fstream linknn;
	linknn.open(flink.c_str()); 
	linknn >> links;
        for(int i=0; i<links ; i++) {     //( id = 0(inside) 1(inlet) 2(outlet) 
		// units are in micro meters
        	linknn >>    id >> l1[i] >> l2[i]  >>      lid[i]     >> ld[i]  >> volth[i] >> shfnl[i] >> clay >> n1l >>  n2l >> ll[i] >> distance;
	}

        //exit(0); 
	linknn.close();
	return;
}

int main(int argc, char *argv[])
{
      double sigma= 5.0e-2; //kg/s^2
      int n,nodes,links,count, i1, i2,inode;
      double xt, yt, zt;

      fstream rnode;
      rnode.open(fnode.c_str());
      rnode >> nodes >> xt >> yt >> zt;
      
      double* nd       = (double*)malloc(nodes * sizeof(double*));
      double* shfrn    = (double*)malloc(nodes * sizeof(double*));
      double* voln     = (double*)malloc(nodes * sizeof(double*));
      
      std::vector<int> ni;       // # of connections
      ni.reserve(nodes);
      // read node locations and diameters
      readnode(nd, ni, shfrn, voln);
       
      double maxd = 10.0;  // maximum diameter
      double maxl = 10.0;  // maximum diameter
      int maxcon = 0; 
      int mincon = 100; 
      double avgcon = 0;  
      //maxd = 100.0;  // center to center distance =  100 nm

      //read link information
      fstream rlink;
      rlink.open(flink.c_str());
      rlink >> links;   // number of links
      
      int* l1      = (int*)malloc(links * sizeof(int*));       // node at one end of link 
      int* l2      = (int*)malloc(links * sizeof(int*));       // node at the other end of the link 
      double* ld   = (double*)malloc(links * sizeof(double*)); // link diameter 
      double* ll   = (double*)malloc(links * sizeof(double*)); // link length  
      double* voll = (double*)malloc(links * sizeof(double*)); // link length  
      double* shfrl = (double*)malloc(links * sizeof(double*)); // link length  
      int* lid      = (int*)malloc(links * sizeof(int*));       // link id 

      double linkl;

      double sumcond,tmp,tmpi,tmpj,tmpl,tmps;
      double rho,cnstc, pavg;
      double len = 100.0; //1e-7;         // link length
      int k1, k2,ic,id; 
      
      int maxIters=90000;    
      double resid;
      int iters;
      
      double massin=0.0,massout=0.0,tmass;
      double vcomp,effr, voli;
      double rid, rjd, qid, qjd;
      double klef, ckkb = 1e7;
      //read link information
      readlink(l1, l2, ld, ll, lid, voll, shfrl);

      //calculate Average coordination number
      int ison = 0;  
      for (int i=0; i<nodes; i++) { 
	      if(ni[i] < mincon) mincon = ni[i];   
	      if(ni[i] > maxcon) maxcon = ni[i];
	      if(ni[i] <= 1) ison = ison + 1;  
	      avgcon = avgcon + ni[i]; 
      }

      double avgshfn = 0.0; 
      //calculate Average shafe factor of pore bodies 
      for (int i=0; i<nodes; i++) { 
	      avgshfn = avgshfn + shfrn[i]; 
      }
 
      // pore body minimum radius 
      double minradn = 1000.0; 
      double maxradn = -1000.0; 
      double avgradn = 0.0; 

      for (int i=0; i<nodes; i++) {
	     if(nd[i] < minradn)  minradn = nd[i];  
	     if(nd[i] > maxradn)  maxradn = nd[i]; 
	     avgradn  = avgradn + nd[i]; 
      }  

      cout <<"----------------------------------------" << endl; 
      cout << " Number of pore bodies " << nodes << endl; 
      cout <<"----------------------------------------" << endl; 
      cout << " Number of throats     " << links << endl; 
      // pore body diameter information
      cout <<"----------------------------------------" << endl;  
      cout <<" pore body:  minimum radius (\mu m) " << minradn  << endl; 
      cout <<" pore body:  maximum radius (\mu m) " << maxradn  << endl; 
      cout <<" pore body:  average radius (\mu m) " << avgradn/double(nodes) << endl; 
      cout <<"----------------------------------------"  << endl; 
      double minradl = 1000.0; 
      double maxradl = -1000.0; 
      double avgradl = 0.0; 

      for (int i=0; i<links; i++) {
	     if(ld[i] < minradl)  minradl = ld[i];  
	     if(ld[i] > maxradl)  maxradl = ld[i]; 
	     avgradl  = avgradl + ld[i]; 
      }  

      // pore throat diameter information 
      cout <<" pore throat:  minimum radius (\mu m) " << minradl  << endl; 
      cout <<" pore throat:  maximum radius (\mu m) " << maxradl  << endl; 
      cout <<" pore throat:  average radius (\mu m) " << avgradl/double(links) << endl; 

      int inletn  = 0;
      int outletn = 0; 
      double avgshfl = 0; 
     //# of nodes connected to the inlet  or outlet 254 0 254
      for (int i=0; i<links; i++) { 
	      k1 = l1[i]; 
	      k2 = l2[i]; 
	      if( lid[i] == 1 )  inletn = inletn  + 1; 
	      if( lid[i] == 2)   outletn = outletn + 1; 
	      avgshfl = avgshfl + shfrl[i]; 
      }

      double totvn=0.0,totvl=0.0; 
      for (int i=0; i<links; i++) {
              totvl    = totvl +  PI*pow(ld[i],2)*ll[i]; 
      }
      for (int i=0; i<nodes; i++) {
              totvn    = totvn +  (4.0/3.0)*PI*pow(ld[i],3); 
      }



      cout <<"----------------------------------------" << endl; 
      cout << " Average shape factor of pore bodies  "  << avgshfn/double(nodes)  << endl; 
      cout << " Average shape factor of pore throats " << avgshfl/double(links)  << endl; 

      cout <<"----------------------------------------" << endl; 

      cout << " Number of links connected to the intlet " << inletn   << endl; 
      cout << " Number of links connected to the outlet " << outletn  << endl;  
      cout <<"----------------------------------------" << endl;

      cout << " Minimum coordination number " << mincon << endl; 
      cout << " Maximum coordination number " << maxcon << endl; 
      cout << " Average coordination number " << avgcon/double(nodes) << endl; 
      cout << " Isolated nodes              " << ison << endl; 
      //cout << " Porosity                    " << (totvl+totvn)/(xt*yt*zt) << endl; 
      cout << " Porosity                    " << (totvl)/(xt*yt*zt) << endl; 
      cout <<"----------------------------------------" << endl;
      cout <<" X dimension (mm)" << xt/1e3 << endl; 
      cout <<" Y dimension (mm)" << yt/1e3 << endl;
      cout <<" Z dimension (mm)" << zt/1e3 << endl;    

      free(nd);
      free(shfrn);

      free(l1);
      free(l2);
      free(ld);
      free(ll);
      free(voll);
      free(shfrl);
      free(lid);
      
      return 0; 
} 
