
/*
 * Author: V Kotteda
 * Date  : October 28, 2021
 *
 * References
 * 1. Huang, X., Bandilla, K. W., & Celia, M. A. (2016).
 * Multi-physics pore-network modeling of two-phase shale matrix flows. Transport in Porous Media, 111(1), 123-141.
 * https://doi.org/10.1007/s11242-015-0584-8
 *
 * 2. JOEKAR-NIASAR, V., HASSANIZADEH, S., & DAHLE, H. (2010).
 * Non-equilibrium effects in capillarity and interfacial area in two-phase flow: Dynamic pore-network modelling. Journal of Fluid Mechanics, 655, 38-71.
 * https://doi.org/10.1017/S0022112010000704
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

// petsc header files 
#include <petscksp.h>
#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscviewer.h>

// include linear solver header file 
#include "linear_solver_petsc.h"
#include "write_vtk.h"
#include "sphase_unsteady_comp.hpp"

using namespace std;
const double PI  =3.141592653589793238463;
double nonle = 1e-9;
double m2nm  = 1.0/nonle;
double pnond = 1e6; 
const int maxrand = 50000;  

class Node
{
        public:
                Node()
                {}
                ~Node()
                {}
                void readdata()
                {
                        ifstream nodeloc("node.dat");
                        nodeloc >> nn; 
                        ifstream noderad("node_dia.dat");
                        for(int i=0; i<nn; i++){ 
                            nodeloc >> xx >> yy >> zz; 
                            noderad >> radius;
                            rad.push_back(radius);
                            xc.push_back(xx);
                            yc.push_back(yy);
                            zc.push_back(zz);
                        }
                        nodeloc.close(); 
                        noderad.close(); 
                }
    
                void DisplayVectorContents()
                {
                        for( unsigned int i = 0; i < rad.size(); i++ )
                        {
                                std::cout << "Element[" << i << "] = " << rad[i] << std::endl;
                        }
                        std::cout << std::endl;
                }

        //private:
                std::vector<double> rad;
                std::vector<double> xc;
                std::vector<double> yc;
                std::vector<double> zc;

         private:
                double xx, yy, zz, radius;
                int nn; 
};


class Link
{
        public:
                Link() {}
                ~Link() {}
                void readdata()
                {
                        ifstream linkinfo("link.dat");
                        linkinfo >> nn;
                        ifstream linkdia("link_dia.dat");
                        linkdia >> nl;
                        for(int i=0; i<nn; i++){
                            linkinfo >> l1 >> l2; 
                            linkdia >> radius; 
                            ends.push_back(std::make_pair(l1, l2)); 
                            rad.push_back(radius);
                        }
                        linkinfo.close();
                        linkdia.close();
                }

                void DisplayVectorContents()
                {
                        for( unsigned int i = 0; i < ends.size(); i++ )
                        {
                                std::cout << "Element[" << i << "] = " << ends.at(i).first <<"\t" << ends.at(i).second << std::endl;
                        }
                        std::cout << std::endl;
                }

        //private:
                std::vector<pair<int,int>> ends;
                std::vector<double> rad;
         private:
                int l1, l2 ;
                int nn, nl; 
                double radius; 
};


/*
void readnode(double* xn, double* yn, double* zn, double* dia) {
	int nn;
	ifstream nodeloc("node.dat");
	nodeloc >> nn;
	
	ifstream nodedia("node_dia.dat");
	for(int i=0; i<nn ; i++) {
	     nodeloc >> xn[i] >> yn[i]>> zn[i];	
	     nodedia >> dia[i]; // in  nm 
	}
	
	nodeloc.close();
	nodedia.close(); 
	return;
}

void readlink(int* l1, int* l2, double* ld) {
	int nl;
	ifstream linknn("link.dat");
	linknn >> nl;
	
	ifstream linkdia("link_dia.dat");
	linkdia >> nl;
	for(int i=0; i<nl ; i++) {
	     linknn >> l1[i] >> l2[i]; 
             linkdia >> ld[i];
	     //ld[i] = 100.0;  //1e-7;  // link diameter 
	}
	
	linknn.close();
	return;
}

*/
double conductance( double radius, double vis, double len) 
{
	//double effrad = sqrt(4.0/PI)*radius;  //effective radius 
	double effrad = radius; //sqrt(4.0/PI)*radius;  //effective radius 
	return (PI*pow(effrad,4)/(8.0*vis*len)); 
}

double dv_sorption( double radius, double cnstc, double rho) 
{
	double cmg    = 2.0e-7; //kg/m^2
	return (6.0*pow(2*radius,2)*cmg*cnstc/rho);
}

double calrho(double pres, double nonpres, double temp)
{
	// pressure in atms 
        //methane molecular weight 16.04 g/mol
	return (16.04*1e-3*(pres*nonpres)/(8.314*temp) ); //kg/m^3	
}
		
double calvol(double rad){
	return ((4.0/3.0)*PI*pow(rad,3));
}


int genwet(int nodes, int* wet, int nx, int ny, int nz) {
	int count = 0, sum=0; 
	int rx, ry, inode; 

	for (int i=0; i< nodes; i++) 
		wet[i] = 0; 

        for(int i=0; i< maxrand; i++) {
	      rx  = rand() % 8 + 2; // random number between 2 and 9 
	      ry  = rand() % 8 + 2; // random number between 2 and 9
              for(int k=2; k< nz-1; k++) {
                      inode = k*nx*ny + rx*ny + ry; 
		      wet[inode] = 10; 
	      }
	      for(int j=0; j<nodes ; j++)
		      sum = sum + wet[j]; 
      
	}

return 0; 
}

int main(int argc, char *argv[])
{
      static char help[] = "Solves a linear system in parallel with KSP.\n";
      
      PetscErrorCode ierr;
      ierr = PetscInitialize(&argc,&argv,(char*)0,help);
      
      double sigma= 5.0e-2; //kg/s^2
      int n,nodes,links,count, i1, i2,inode,nsteps;
      double xt, yt, zt, dt,zleft, zright;

      double pinlet ; // Pressure at the inlet  boundary 
      double poutlet; // Pressure at the outlet boundary 
      double mu     ; // viscosity
      double tempe  ;
      int    uniform; // ==1 true 
      double ll     ; //link length 
      double linkr  ; // link radius
      int nx, ny, nz;  

      fstream infile;
      infile.open("input.txt",ios::in); //open a file to perform read operation using file object
      infile >> pinlet ;   // in MPa 
      infile >> poutlet ;  // in MPa 
      infile >> mu ;       // viscosity (Pa. s)
      infile >> dt ;       // time step 
      infile >> nsteps;    // number of time steps
      infile >> tempe;     // temperature in K 
      infile >> uniform;  // uniform 
      infile >> linkr;    // link radius
      infile >> ll;       // link length
      infile >> nx;
      infile >> ny;
      infile >> nz;

      infile >> zleft; 
      infile >> zright; 
      infile.close();      //close the file object.
      double pmin = poutlet, pmax = pinlet;  //MPa 
      ifstream rnode("node.dat");
      rnode >> nodes; 
      
      //double* xc    = (double*)malloc(nodes * sizeof(double*));
      //double* yc    = (double*)malloc(nodes * sizeof(double*));
      //double* zc    = (double*)malloc(nodes * sizeof(double*));
      //double* nd    = (double*)malloc(nodes * sizeof(double*));
      double* nreff = (double*)malloc(nodes * sizeof(double*));
      double* nflow = (double*)malloc(nodes * sizeof(double*));
      
      // read node locations and diameters 
      //readnode(xc, yc, zc, nd); 
      Node Body;  
      Body.readdata();  
      //Body.DisplayVectorContents();  
      
      //exit(0); 
 
      double maxd = 10.0;  // maximum diameter   
      for(int i=0; i<nodes ; i++) {
	      // radius of pore bodies 
              Body.rad[i] = Body.rad[i]/2.0; 
	      //nd[i] = nd[i]/2.0; 
      }
      maxd = 100.0;  // center to center distance =  100 nm 

      //read link information 
      ifstream rlink("link.dat");
      rlink >> links;   // number of links 
     
      Link Throat;  
      Throat.readdata(); 
      //int* l1    = (int*)malloc(links * sizeof(int*));       // node at one end of link 
      //int* l2    = (int*)malloc(links * sizeof(int*));       // node at the other end of the link 
      //double* ld = (double*)malloc(links * sizeof(double*)); // link diameter 

      int* wet     = (int*)malloc(nodes * sizeof(int*));     // hydrocarbon-wet 
      double linkl;  

      std::vector<double> den;   // density 
      den.reserve(nodes); 
      
      std::vector<int> ni;       // # of connections 
      ni.reserve(nodes);
     
      std::vector<double> rates; // flow rate
      rates.reserve(nodes);

      std::vector<double> val;   // non-zero values of matrix 
      val.reserve(links);
      
      std::vector<double> rhs;  // rhs vector
      rhs.reserve(nodes);
      
      std::vector<double> pressures; // unknown vector 
      pressures.reserve(nodes); 
      
      std::vector<double> pold;    // old unknown vector 
      pold.reserve(nodes);
      
      std::vector<int> row;       // non-zero row index of matrix 
      row.reserve(links); 
      
      std::vector<int> col;       // non-zero column index of matrix 
      col.reserve(links); 
      
      double sumcond,tmp,tmpi,tmpj,tmpl,tmps,position,vol,shapeFactor,clay,distance;
      double rho,cnstc, pavg;
      double len = 100.0; //1e-7;         // link length 
      int k1, k2,ic,id;
      
      int maxIters=90000;    
      double resid;
      int iters, nsec=nx*ny;
      
      double length, area;
      double domainvol = Body.xc[nodes]*Body.yc[nodes]*Body.zc[nodes];  
      double massin=0.0,massout=0.0,tmass; 
      double vcomp,effr, voli;
      double rid, rjd, qid, qjd;  
      double porevol = 0, klef, ckkb = 1e7;
      double permeability, flowr, flow;  

      int conn[nodes][10];
      
      //read link information 
      //readlink(l1, l2, ld);
      
      // initialize rhs, and pressure 
      for(int i=0; i<nodes; i++) {
         rhs.push_back(0);
         pressures.push_back(poutlet); // initial pressure 
      }
      
      // apply pressure bc 
      for(int i=0; i<nodes; i++) {
         if(Body.zc[i] < zleft)  {
		 pressures[i] = pinlet; 
	 } 
	 else if(Body.zc[i] > zright) { 
		 pressures[i] = poutlet;
	 }
      }
      
      pold = pressures ;
      //loop through the nodes and make connectivity array 
      for(int i=0; i<nodes; i++) 
      	      ni[i] = 0; 

      // loop through the links
      for(int j=0; j<links; j++) {
	      k1 = Throat.ends.at(j).first;  //l1[j];
              k2 = Throat.ends.at(j).second; //l2[j];
              // check the link end nodes with the node
	      i1 = ni[k1];
	      conn[k1][i1] = j ;
	      ni[k1] = ni[k1]+ 1;

	      i2 = ni[k2]; 
	      conn[k2][i2] = j ; 
	      ni[k2] = ni[k2]+ 1; 
      }

      area   = (ny*ll)*(nz*ll)*nonle*nonle ;  // not in use 
      length = (nx*ll)*nonle ;  // not in use 
      
      ofstream myfile ("pressure_out.txt");
      ofstream outmass ("time_massflow.txt");
     
      // loop to define hydrocarbon/water wet
      genwet(nodes, wet, nx, ny, nz); 

      count = 0; 
      for(int i=0; i<nodes; i++) {
	      if(wet[i]==1) count = count + 1; 
      }
      cout << " % of organic pores " << float(count)/nodes  << endl;  

      //exit(0);  
      // constant in adsorption term  
      cnstc = 1.0/(pmax-pmin)*(0.1*pmax/(1.0+0.1*pmax) - 0.1*pmin/(1.0+0.1*pmin));

      outmass   << 0 <<"\t"<< dt*0 <<"\t" << massin*dt*1e-9*0 <<"\t"<< -massout*dt*1e-9*0 <<"\t" << std::endl;

      // Non-equilibrium effects in capillarity and interfacial area in two-phase flow: dynamic pore-network modelling 
      // https://doi.org/10.1017/S0022112010000704
      // calculate inscribed radius of pore throat 

      if(uniform ==1) {  
	      for(int j=0; j<links; j++) { 
		      Throat.rad[j] = linkr; 
                      // ld[j] = linkr;  
	      }
      }
      else 
      {
             for(int j=0; j<links; j++) {
                 k1     = Throat.ends.at(j).first;  //l1[j];
	         k2     = Throat.ends.at(j).second; //l2[j];
	         rid    = Body.rad[k1]/maxd ;  //nd[k1]/maxd;
	         rjd    = Body.rad[k2]/maxd ;  //nd[k2]/maxd;
                 qid    = rid*sin(PI/4.0)/pow((1.0-rid*cos(PI/4.0)),0.2);
                 qjd    = rjd*sin(PI/4.0)/pow((1.0-rjd*cos(PI/4.0)),0.2); 
                 Throat.rad[j] = maxd*qid*qjd*pow((pow(qid,1.0/0.2)+pow(qjd,1.0/0.2)),-0.2);
                 //ld[j] = maxd*qid*qjd*pow((pow(qid,1.0/0.2)+pow(qjd,1.0/0.2)),-0.2);
              }


              for(int j=0; j<links; j++) {
                      k1 = Throat.ends.at(j).first; //l1[j]; 
                      k2 = Throat.ends.at(j).second; //l2[j]; 
                      //if(ld[j] > nd[k1]) 
                      if(Throat.rad[j] > Body.rad[k1]) 
                	      Throat.rad[j] = Body.rad[k1]/3.0; 
                	      //ld[j] = nd[k1]/3.0; 
                      else if (Throat.rad[j] > Body.rad[k2]) // nd[k2]) 
                	      Throat.rad[j] = Body.rad[k2]/3.0; 
                	      //ld[j] = nd[k2]/3.0; 
              }
      }                                        
      // calculate initial density 
      for(int i=0; i<nodes; i++) {
	      den[i]    = calrho(pressures[i],pnond, tempe); 
      }

      // calculate porosity
      double porosity = 0.0;
      for(int i=0; i<nodes; i++)
              porosity = porosity + (4.0/3.0)*PI*pow(Body.rad[i],3); // nd[i],3); // add nodes volume 
      for(int i=0; i<links; i++)
              porosity = porosity + PI*pow(Throat.rad[i],2)*100.0; //ld[i],2)*100.0;     // add links volume 

      porosity =  porosity/(1000*1000*1000);
      cout << "porosity" << porosity << endl;
      ofstream outavgpfile ("z_avg_pressure_out.txt");
      double zpresavg[nodes/nsec];
      int indx;
      /*******************************************************************************/
      // loop through time steps 
      for(int k=1; k<nsteps; k++) {
	      printf("time step loop # %i\n", k);
             //calculate effective radius (update it at hydro-carbon wet pores)  
	     for(int i=0; i<nodes; i++) {
		     if(wet[i] == 1) {
                     // gas sorption - modification to effective pore radius 
                     effr = 0.0; //0.4*0.1*pressures[i]/(1.0+0.1*pressures[i]);
                     nreff[i]  = Body.rad[i] - effr; //nd[i] - effr; 
	 	     }
	 	     else {
	 		     nreff[i] = Body.rad[i] ; //nd[i];
	             }
	     }

             //loop through the nodes
             for(int i=0; i<nodes; i++) {
		     sumcond  = 0;
                     rhs[i]   = 0;
	             nflow[i] = 0;
      	             voli     = calvol(Body.rad[i]); //calvol(Body.rad[i]); // nd[i]); 
                     // loop through the links
                     if(ni[i] ==0) {  //shape factor = 0  
	   		     val.push_back(1.0);
	   		     row.push_back(i);
			     col.push_back(i);
			     rhs[i] = poutlet; 
		     }
		     else if(Body.zc[i] <zleft) { 
			     rhs[i] = pinlet; 
			     val.push_back(1.0);
			     row.push_back(i);
			     col.push_back(i);
		     }
		     else if(Body.zc[i] > zright) {
			     rhs[i] = poutlet;
			     val.push_back(1.0);
			     row.push_back(i);
			     col.push_back(i);
		     }
		     else { 
	   		     tmps = 0.0; 
		      	     //den[i]    = calrho(pressures[i],pnond); 
			     //if(wet[i] == 1) tmps = dv_sorption(nreff[i], cnstc, den[i])/(m2nm);
			     for(int j=0; j<ni[i]; j++) {
				     id  = conn[i][j]; // link id 
				     k1  = Throat.ends.at(id).first;   //l1[id];     // node at one end of the link 
				     k2  = Throat.ends.at(id).second;  //l2[id];     // node on the other end of the link 
				     {
					     linkl = ll- nreff[k1] - nreff[k2]; // link length 	
					     //effr  = 0.4*0.1*pressures[i]/(1+0.1*pressures[i]); // pressure in MPa
					     //only pore throats offer resistence to flow but they do not have volume 
					     //pore bodies possess volumes but they do not provide resistance to flow   
					     effr  = 0.0; //  0.4*0.1*pavg/(1.0+0.1*pavg);                 // link effective radius 
					     tmpl  = conductance(Throat.rad[id]-effr, mu,linkl)*((pressures[k1]+pressures[k2])/(2.0*pressures[i])) ; // compressible 
					     //tmpl  = conductance(ld[id]-effr, mu,linkl)*((pressures[k1]+pressures[k2])/(2.0*pressures[i])) ; // compressible 
					     tmp   = tmpl*dt; //1.0/(1.0/tmpi + 1.0/tmpl + 1.0/tmpj);
					     //klef = (1.0+2.0*(ckkb/pnond)/(pressures[k1]+pressures[k2])); 

					     sumcond = sumcond + tmp; //*klef;
          	           
					     // Applying pressure inlet  condition 
					     if(Body.zc[k1] < zleft || Body.zc[k2] < zleft) {
						     //  if any end connected to inlet 
						     rhs[i] = rhs[i] - pinlet*tmp;
					     }
					     // Applying pressure outlet condition 
					     else if(Body.zc[k1] > zright || Body.zc[k2] > zright) {
						     // if any end connected to outlet 
						     rhs[i] = rhs[i] - poutlet*tmp;
					     }
					     else {
						     if(i==k1) {   // k1 -node number  >0  internal node 
							     val.push_back(tmp);
							     row.push_back(i); 
							     col.push_back(k2);
							     nflow[i]  = nflow[i] + tmp*(pressures[i] - pressures[k2]);   
						     }
						     else if (i==k2) { // k2 - node number  > 0 
							     val.push_back(tmp);
							     row.push_back(i);
							     col.push_back(k1);
							     nflow[i]  = nflow[i] + tmp*(pressures[i] - pressures[k1]);   
						     }
                           
					     }
		      
				     } // if condition 
                  
			     }  // connection loop

			     vcomp  = voli/pnond; // compressibility term  
			     if(sumcond < 1.e-200  && sumcond > -1.e-200) sumcond = 1.0;
			     sumcond = sumcond + tmps;     // add sorption term  		  
			     val.push_back(-sumcond-vcomp/pressures[i]); // add compressibility term 
			     rhs[i] = rhs[i]-pressures[i]*tmps- vcomp ; //modify right hand side vector to consider sorption and compressibility 
			     row.push_back(i);
			     col.push_back(i);

			     }  // if condition loop end 
       	     }   //for node loop 
      
	     pold = pressures; 
	     linear_solver_petsc(val,
                            rhs,
                            pressures,
                            row,
                            col,
                            nodes,
                            maxIters,
                            resid,
                            iters); 

        
	     //density calculations
	     for(int i=0; i<nodes; i++) {
		     voli    = calvol(Body.rad[i]);  //nd[i]); // (4.0/3.0)*PI*pow(nd[i],3);
		     tmp     = den[i]; //kg/m^3  // density
		     tmps    = 0.0;
		     if(wet[i] == 1) tmps = 0.0; // dv_sorption(nreff[i], cnstc, den[i])/(m2nm);
		     // saturation = 1.0 
		     den[i]  = (tmp*(nflow[i] + tmps*(pold[i]-pressures[i]))*pnond + tmp*voli*1.0)/ (1.0*voli) ;
     	     }   //for node loop


        
	     // calculate massflow at the inlet 
	     // massin = 0.0; 
	     for(int i=0; i<nsec; i++) {
		     for(int j=0; j<ni[i]; j++) {
			     id     = conn[i][j];
			     k1     = Throat.ends.at(id).first ; //l1[id];
			     k2     = Throat.ends.at(id).second ;; //l2[id];
			     linkl  = ll - nreff[k1] - nreff[k2]; 
			     pavg   = (pressures[k1] + pressures[k2])/2.0;  
			     //effr = 0.4*0.1*pressures[i]/(1+0.1*pressures[i]); // pressure in MPa
			     effr   = 0.0; //0.4*0.1*pavg/(1.0+0.1*pavg); // pressure in MPa
			     tmpl   = conductance(Throat.rad[id]-effr, mu,linkl)*((pressures[k1]+pressures[k2])/(2.0*pressures[i])) ; // compressible
			     //tmpl   = conductance(ld[id]-effr, mu,linkl)*((pressures[k1]+pressures[k2])/(2.0*pressures[i])) ; // compressible
			     tmp    =  tmpl*dt; //1.0/(1.0/tmpi + 1.0/tmpl + 1.0/tmpj);
			     if(i==k1) massin = massin + tmp * den[i]*(pressures[k1] - pressures[k2]);
                             if(i==k1) flowr  = flowr  + tmpl *(pressures[k1] - pressures[k2]);
			     if(i==k2) massin = massin + tmp * den[i]*(pressures[k2] - pressures[k1]);
                             if(i==k2) flowr  = flowr  + tmpl *(pressures[k2] - pressures[k1]);
			     
		     }
	     }
	     // calculate massflow at the outlet 
	     // massout = 0.0; 
	     for(int i=nodes-nsec; i<nodes; i++) {
		     if(ni[i] != 0) {
			     for(int j=0; j<ni[i]; j++) {
				     id    = conn[i][j];
				     k1    = Throat.ends.at(id).first;    //l1[id];
				     k2    = Throat.ends.at(id).second;  //l2[id];
				     linkl = ll- nreff[k1] - nreff[k2]; 
				     pavg  = (pressures[k1] + pressures[k2])/2.0;  
				     //effr  = 0.4*0.1*pressures[i]/(1+0.1*pressures[i]); // pressure in MPa
				     effr  = 0.0; //0.4*0.1*pavg/(1.0+0.1*pavg); // pressure in MPa
				     tmpl  = conductance(Throat.rad[id]-effr, mu,linkl)*((pressures[k1]+pressures[k2])/(2.0*pressures[i])) ; // compressible
				     //tmpl  = conductance(ld[id]-effr, mu,linkl)*((pressures[k1]+pressures[k2])/(2.0*pressures[i])) ; // compressible
				     tmp   = dt*tmpl; //1.0/(1.0/tmpi + 1.0/tmpl + 1.0/tmpj);
				     if(i==k1) massout = massout + tmp * den[i]*(pressures[k1] - pressures[k2]); 
				     if(i==k2) massout = massout + tmp * den[i]*(pressures[k2] - pressures[k1]); 
			     }
		     }
	     }
	     cout << "time step "<< k <<" "<<k*dt<<" "<< massin << " "<<massout << endl; 
	     //cout << "avg flow rate  "<< flow <<" "<< flow*mu*length/area/(pinlet-poutlet) << std::endl;
	     outmass   << k <<"\t"<< k*dt <<"\t" << massin*1e-21 <<"\t"<< massout*1e-21 <<"\t" << std::endl;

             permeability =  flowr*nonle*nonle*mu*length/(area*(pinlet-poutlet)*pnond);
//             cout << "permeability " <<  flowr*1e-12  <<"\t"<< mu <<"\t"<< area <<"\t"<< length <<"\t"<< (pinlet-poutlet)*pnond <<"\t" << permeability << endl;
             cout << "permeability " <<  permeability << endl;   ///(0.987*1e-18) << endl;
      } // time step loop 
      
     
      for(int i=0; i<nodes; i++) {
	      myfile << pressures[i] <<" "<<Body.xc[i] <<" "<<Body.yc[i]<<" "<<Body.zc[i]<<" "<<  std::endl;
      }
      myfile.close();
      
      //printNodesVTK(nodes, Body.xc, yc, zc, nd, pressures); 
      //printLinksVTK(nodes, links, xc, yc, zc, l1, l2, ld, pressures); 

      cout << "print vtp files: completed" << endl; 

      outavgpfile.close();

      system("gnuplot plot_mass.gnu");
      //system("gnuplot -p -e plot_mass.gnu");
//      free(xc);
//      free(yc);
//      free(zc);
//      free(nd);
      free(nreff);
      free(nflow);
      free(wet); 
      
//      free(l1);
//      free(l2);
//      free(ld);
      
      ierr = PetscFinalize();
      if (ierr) {return ierr;}
} 
