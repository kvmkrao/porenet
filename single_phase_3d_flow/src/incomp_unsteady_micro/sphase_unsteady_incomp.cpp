/*
 * Author: V Kotteda
 * Date  : Aug 19, 2021 
 *
 */

#include<cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <sstream>
#include <vector>
#include <iomanip>    
#include <stdlib.h>   
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

using namespace std;
const double PI  =3.141592653589793238463;
double nonle = 1e-6;
double m2nm  = 1.0/nonle;
double pnond = 1e6; 
double dt    = 5e-6 ;
int nsteps   = 1000; //1e-7/dt;  // number of time steps 
double nondmass = nonle*nonle*nonle*pnond ; // non-dimensional mass

void readnode(string fnode, double* dia, int* ni) {
	size_t nodes, id;
	size_t position, shapf; 
	double nr, vol, shapeFactor, clay; 
	double dd1, dd2, dd3,dd4; 
	fstream nodeloc;
	nodeloc.open(fnode.c_str()); 
        nodeloc >> nodes >>  dd1 >>  dd2 >> dd3 >> dd4;
        cout.precision(3);
        string str1, str2, str3, str4,str5;
	for(int i=0; i<nodes ; i++) {  
		// units are in nano meters
	        nodeloc >> id >> str1 >> str2 >> str3 >> position >> str5 >> str4 >> vol >> shapeFactor >> clay;
                //xn[i] = stod(str1);  
                //yn[i] = stod(str2);  
                //zn[i] = stod(str3);  
                dia[i] = stod(str4);  
		ni[i]  = stoi(str5); 
        }
	nodeloc.close();
	return;
}

void readlink(string flink, int* l1, int* l2, double* ld, double* ll, int* lid) {
	int links,lp,id;
        double vol, shapeFactor, clay, n1l, n2l, distance; 
        
	fstream linknn;
	linknn.open(flink.c_str()); 
	linknn >> links;
        for(int i=0; i<links ; i++) {     //( id = 0(inside) 1(inlet) 2(outlet) 
		// units are in nano meters
        	linknn >>    id >> l1[i] >> l2[i]  >>  lid[i]  >> ld[i]  >> vol >> shapeFactor >> clay >> n1l  >>  n2l  >> ll[i]  >> distance;
	}
	linknn.close();
	return;
}

double conductance( double radius, double vis, double len) 
{
	double effrad = radius; // sqrt(4.0/PI)*radius;  //effective radius 
	return (PI*pow(effrad,4)/(8.0*vis*len)); 
}

double dv_sorption( double radius, double cnstc, double rho) 
{
	double cmg    = 2.0e-7; //kg/m^2
	return (6.0*pow(2*radius,2)*cmg*cnstc/rho);
}

double calvol(double rad){
	return ((4.0/3.0)*PI*pow(rad,3));
}

int main(int argc, char *argv[])
{
      static char help[] = "Solves a linear system in parallel with KSP.\n";
      
      PetscErrorCode ierr;
      ierr = PetscInitialize(&argc,&argv,(char*)0,help);
      
      double sigma= 5.0e-2; //kg/s^2
      int n,nodes,links,count, i1, i2,inode;
      double xt, yt, zt;
      double dt   ; // time step
      int nsteps  ; // number of time steps

      double pinlet ; // Pressure at the inlet  boundary 
      double poutlet; // Pressure at the outlet boundary 
      double mu     ; // viscosity
      string fnode ; 
      string flink ; 

      fstream infile;
      infile.open("input.txt",ios::in); //open a file to perform read operation using file object
      infile >> pinlet ;   // in MPa 
      infile >> poutlet ;  // in MPa 
      infile >> mu ;       // viscosity (Pa. s)
      infile >> dt ;       // time step 
      infile >> nsteps;    // number of time steps 
      infile >> fnode;     // node.dat path 
      infile >> flink;     // link.dat path 
      infile.close();      //close the file object.

      double pmin = poutlet, pmax = pinlet; //MPa 

      fstream rnode;
      rnode.open(fnode.c_str()); 
      rnode >> nodes >> xt >> yt >> zt; 
      
      double* nd    = (double*)malloc(nodes * sizeof(double*));
      int* ni       = (int*)malloc(nodes * sizeof(int*));
      
      // read node locations and diameters 
      readnode(fnode, nd, ni); 
       
      double maxd = 10.0;  // maximum diameter   
      double maxl = 10.0;  // maximum diameter  
      int maxcon = 0;  
      //maxd = 100.0;  // center to center distance =  100 nm 

      //read link information 
      fstream rlink;
      rlink.open(flink.c_str()); 
      rlink >> links;   // number of links 
      
      int* l1    = (int*)malloc(links * sizeof(int*));       // node at one end of link 
      int* l2    = (int*)malloc(links * sizeof(int*));       // node at the other end of the link 
      double* ld = (double*)malloc(links * sizeof(double*)); // link diameter 
      double* ll = (double*)malloc(links * sizeof(double*)); // link length  
      int* lid   = (int*)malloc(links * sizeof(int*));       // link id 

      double linkl;

      std::vector<double> val;   // non-zero values of matrix 
      val.reserve(links);
      
      std::vector<double> rhs;  // rhs vector
      rhs.reserve(nodes);
      
      std::vector<double> pressures; // unknown vector 
      pressures.reserve(nodes); 
      
      std::vector<double> pold;   // old unknown vector 
      pold.reserve(nodes);

      std::vector<double> perr;   // old unknown vector 
      perr.reserve(nodes);
      std::vector<int> row;       // non-zero row index of matrix 
      row.reserve(links); 
      
      std::vector<int> col;       // non-zero column index of matrix 
      col.reserve(links); 

      double sumcond,tmp,tmpi,tmpj,tmpl,tmps;
      double rho,cnstc, pavg;
      double len = 100.0; //1e-7;         // link length 
      int k1, k2,ic,id;
      
      int maxIters=90000;    
      double resid;
      int iters;

      double length, area;
      double permeability, flowr, flow;
      double massin=0.0,massout=0.0,tmass; 
      double vcomp,effr, voli;
      double rid, rjd, qid, qjd; 
      double klef, ckkb = 1e7;
      double l2norm; 
      //read link information 
      readlink(flink,l1, l2, ld, ll, lid);

      ofstream tout ("output.txt");
 
      // initialize rhs, and pressure 
      for(int i=0; i<nodes; i++) {
         rhs.push_back(0.0);
         pressures.push_back(poutlet);
      }
      //loop through the nodes and make connectivity array 
      for(int i=0; i<nodes; i++) 
	      if(ni[i] > maxcon) maxcon = ni[i];   
      tout << "maximum co-ordination number " << maxcon <<endl; 
   
// Declare memory block of size nodes
      int** conn = new int*[nodes];
      for (int i=0; i<nodes; ++i) 
	    conn[i] = new int[ni[i]];   // Declare a memory block of size ni[i]

      area = PI*(yt/2.0)*(yt/2.0); 
      length = xt; 
      ofstream myfile ("pressure_out.txt");
      ofstream outmass ("time_massflow.txt");
      // constant in adsorption term  
      cnstc = 1.0/(pmax-pmin)*(0.1*pmax/(1.0+0.1*pmax) - 0.1*pmin/(1.0+0.1*pmin));

      outmass   << 0 <<"\t"<< dt*0 <<"\t" << massin*dt*1e-9*0 <<"\t"<< -massout*dt*1e-9*0 <<"\t" << std::endl;
/*      
 * JOEKAR-NIASAR, V., HASSANIZADEH, S., & DAHLE, H. (2010). 
 * Non-equilibrium effects in capillarity and interfacial area in two-phase flow: Dynamic pore-network modelling. Journal of Fluid Mechanics, 655, 38-71. 
 * https://doi.org/10.1017/S0022112010000704 
 */
      // calculate inscribed radius of pore throat 

      /*
      for(int j=0; j<links; j++) {
              if(k1 != 0) { 
              k1     = l1[j]-1;
	      k2     = l2[j]-1;
	      rid    = nd[k1]/maxd;
	      rjd    = nd[k2]/maxd;
              qid    = rid*sin(PI/4.0)/pow((1.0-rid*cos(PI/4.0)),0.2);
              qjd    = rjd*sin(PI/4.0)/pow((1.0-rjd*cos(PI/4.0)),0.2); 
              ld[j]  = maxd*qid*qjd*pow((pow(qid,1.0/0.2)+pow(qjd,1.0/0.2)),-0.2);
	      }
      }
      */

      pold = pressures ;
      //loop through the nodes and make connectivity array 
      for(int i=0; i<nodes; i++)
              ni[i] = 0;

      //loop through the links
      for(int j=0; j<links; j++) {
              k1 = l1[j];
              k2 = l2[j];
              // check the link end nodes with the node
              if(k1 > 0) { //node1 ==0 : inlet or outlet 
                      i1 = ni[k1-1];
                      conn[k1-1][i1] = j ;
                      ni[k1-1] = ni[k1-1]+ 1;
              }

              i2 = ni[k2-1];
              conn[k2-1][i2] = j ;
              ni[k2-1] = ni[k2-1]+ 1;
      }

      // set the link length to 100nm 
/*      for(int j=0; j<links; j++) {
              ll[j] = 100.0;
      }
*/
      // calculate porosity
      double porosity = 0.0;
      for(int i=0; i<nodes; i++)
              porosity = porosity + (4.0/3.0)*PI*pow(nd[i],3); // add nodes volume 
      for(int i=0; i<links; i++)
              porosity = porosity + PI*pow(ld[i],2)*ll[i];     // add links volume 

      porosity =  porosity/(PI*xt*(yt/2.0)*(zt/2.0));
      tout << "porosity" << porosity << endl;
      /*******************************************************************************/
      // loop through time steps 
      for(int k=1; k<nsteps; k++) {
	      tout << "time step loop # " << k << endl; 
             //calculate effective radius (update it at hydro-carbon wet pores)  
             // loop through the nodes
             for(int i=0; i<nodes; i++) {
		     sumcond  = 0.0;
                     rhs[i]   = 0.0;
      	             voli     = calvol(nd[i]);
                     // loop through the links
                     if(ni[i] ==0) {  //shape factor = 0  // isolated pore body  
	   		     val.push_back(1.0);
	   		     row.push_back(i);
			     col.push_back(i);
			     rhs[i] = poutlet; 
		     }
		     else {       // # of connections >  0 
	   		     tmps = 0.0; 
			     for(int j=0; j<ni[i]; j++) {
				     id  = conn[i][j]; // link id 
				     k1  = l1[id];     // node at one end of the link 
				     k2  = l2[id];     // node on the other end of the link 
				     if(k1 == 0 && i == k2-1){      // check whether it is connected to inlet/outlet 
					     if(lid[id] == 1) { //  connected to inlet 
						     linkl = ll[id]; // -  nd[k2-1]; // link length 
						     effr  = 0.0; 
						     tmpl  = conductance(ld[id], mu,linkl);
						     tmp   = tmpl*dt;  
						     rhs[i] = rhs[i] - pinlet*tmp;
						     sumcond = sumcond + tmp; 
					     }
					     else if (lid[id] == 2) { // connected to outlet 
						     linkl = ll[id]; 
						     effr  = 0.0; 
						     tmpl  = conductance(ld[id], mu,linkl);
						     tmp   = tmpl * dt;
						     rhs[i] = rhs[i] - poutlet*tmp;
						     sumcond = sumcond + tmp;
					     }     
				     }
				     else {  // internal node 
					     linkl = ll[id]; 
					     effr  = 0.0; //  0.4*0.1*pavg/(1.0+0.1*pavg);                 // link effective radius 
					     tmpl  = conductance(ld[id], mu,linkl); 
					     tmp   = tmpl*dt; 
					     sumcond = sumcond + tmp; //*klef;
					     if(i==k1-1) {   // k1 -node number  >0  internal node 
						     val.push_back(tmp);
						     row.push_back(i); 
						     col.push_back(k2-1);
					     } 
					     else if (i==k2-1) { // k2 -node number  > 0  
						     val.push_back(tmp); 
						     row.push_back(i); 
						     col.push_back(k1-1); 
                                             }
		      
				     } // if condition 
                  
			     }  // connection loop

			     if(sumcond < 1.e-200  && sumcond > -1.e-200) sumcond = 1.0;
			     val.push_back(-sumcond);
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
			    ni, 
                            nodes,
			    maxcon, 
                            maxIters,
                            resid,
                            iters);

             for(int i=0; i<nodes; i++) 
			perr[i] = abs(pold[i] - pressures[i]); 

             l2norm = -1e-6; 
             for(int i=0; i<nodes; i++)
                     if(perr[i] > l2norm) l2norm = perr[i]; 

             flowr = 0.0;  
             for(int i=0; i<links; i++) {
                     id     = i;
                     k1     = l1[id];
                     k2     = l2[id];
                     if(k1 == 0) { 
                	     if(lid[id] == 1) {   // inlet 
                		     linkl  = ll[id]; // - nd[k2-1]; 
                	             effr   = 0.0; //0.4*0.1*pavg/(1.0+0.1*pavg); // pressure in MPa
                		     tmpl   = conductance(ld[id], mu,linkl); //*((pinlet+pressures[k2-1])/(2.0*pressures[k2-1])) ; // compressible 
                		     tmp    = tmpl*dt; 
                		     massin = massin + tmp *(pinlet - pressures[k2-1]);
                		     flowr =  flowr  + tmpl *(pinlet - pressures[k2-1]);
                	     }
                	     else if(lid[id] == 2) {    // outlet	            	     // calculate massflow at the outlet 
                		     linkl = ll[id]; //  - nd[k2-1]; 
                		     effr  = 0.0; //0.4*0.1*pavg/(1.0+0.1*pavg); // pressure in MPa
                		     tmpl  = conductance(ld[id]-effr, mu,linkl); //*((poutlet+pressures[k2-1])/(2.0*pressures[k2-1])) ; // compressible
                		     tmp   = tmpl*dt; 
                		     massout = massout + tmp *(pressures[k2-1] - poutlet);
                		     flowr   = flowr   + tmpl*(pressures[k2-1] - poutlet);
                	     }
                     }
             }    // mass flow 

             permeability =  flowr*1e-12*mu*length/(area*(pinlet-poutlet)*pnond);
             tout << "permeability " <<  flowr*1e-12  <<"\t"<< mu <<"\t"<< area <<"\t"<< length <<"\t"<< (pinlet-poutlet)*pnond <<"\t" << permeability << endl;
             tout << "permeability " <<  permeability << endl;   ///(0.987*1e-18) << endl;

             tout << "time step "<< k <<" "<<k*dt<<" "<< massin << " "<<massout << "\t"<<  l2norm << endl; 
             outmass   << k <<"\t"<< k*dt <<"\t" << massin*nondmass <<"\t"<< massout*nondmass <<"\t" <<  permeability <<  std::endl;
             //cout << "avg flow rate  "<< flow <<" "<< flow*mu*length/area/(pinlet-poutlet) << std::endl;
             //outmass   << k <<"\t"<< k*dt <<"\t" << massin*1e-21 <<"\t"<< -massout*1e-21 <<"\t" << std::endl;

             if(l2norm < poutlet/1e4 || k== nsteps-1 ) { 
                      tout <<" reached steady state" << std::endl; 
        
	              // calculate massflow at the inlet 
	              // massin  = 0.0;
	              // massout = 0.0;
      
     
                      for(int i=0; i<nodes; i++) {
                              myfile << pressures[i] <<" "<<i << std::endl;
                      }
                      myfile.close();
                      
                      printNodesVTK(fnode, nodes, nd, pressures); 
                      printLinksVTK(fnode, nodes, links, l1, l2, ld, pressures); 

                      tout << "print vtp files: completed" << endl; 
             }  // if condition loop 

     } // time step loop

      free(nd);
      //Delete the array created
      for(int i=0;i<nodes;++i)  //delete the inner arrays
          delete [] conn[i];   
      delete [] conn;          //delete the outer array
      free(ni);

      free(l1);
      free(l2);
      free(ld);
      free(ll);
      free(lid);
      
      ierr = PetscFinalize();
      if (ierr) {return ierr;}
} 
