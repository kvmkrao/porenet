/*
 * Author: V Kotteda
 * Date  : Dec 8, 2021
 * Find the permeability of a 3D network 
 * flow type: incompressible or compressible
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
#include "write_data.h"
#include "node.hpp"
#include "link.hpp"
#include "constants.h"

#include<omp.h> 

using namespace std;
using namespace constants; 

static double conductance( double radius, double vis, double len) 
{
	double effrad = radius; //sqrt(4.0/PI)*radius;  //effective radius 
	return (PI*pow(effrad,4)/(8.0*vis*len)); 
}

static double dv_sorption( double radius, double cnstc, double rho) 
{
	double cmg    = 2.0e-7; //kg/m^2
	return (6.0*pow(2*radius,2)*cmg*cnstc/rho);
}

static double calrho(double pres, double nonpres, double temp)
{
        //methane molecular weight 16.04 g/mol
	return (16.04*1e-3*(pres*nonpres)/(8.314*temp) ); //kg/m^3	
}
		
static double calvol(double rad){
	return ((4.0/3.0)*PI*pow(rad,3));
}

int main(int argc, char *argv[])
{
      static char help[] = "Solves a linear system in parallel with KSP.\n";
      
      PetscErrorCode ierr;
      ierr = PetscInitialize(&argc,&argv,(char*)0,help);
      
      double sigma= 5.0e-2; //kg/s^2
      int n,nodes,links,count, i1, i2,inode,nsteps;
      double xt, yt, zt, dt;

      double pinlet ; // Pressure at the inlet  boundary 
      double poutlet; // Pressure at the outlet boundary 
      double mu     ; // viscosity
      double tempe  ;
      int    uniform; // ==1 true 
      double ll     ; //link length 
      double linkr  ; // link radius
      int nx, ny, nz, flowType; 
      string fnode,flink; 

      fstream infile;
      infile.open("input.txt",ios::in); //open a file to perform read operation using file object
      infile >> pinlet ;   // in MPa 
      infile >> poutlet ;  // in MPa 
      infile >> mu ;       // viscosity (Pa. s)
      infile >> dt ;       // time step 
      infile >> nsteps;    // number of time steps
      infile >> tempe;     // temperature in K
      infile >> flowType;  // flow type (1 - incompressible; >1 for compressible)   
      infile.close();      //close the file object.

      double pmin = poutlet, pmax = pinlet;  //MPa
      ifstream rnode;  
      rnode.open(constants::fnode1);
      rnode >> nodes >> xt >> yt >> zt;  
      
      double* nreff = (double*)malloc(nodes * sizeof(double*));
      double* nflow = (double*)malloc(nodes * sizeof(double*));

      //Node Body;  
      OrenNode Bodyadd;
      OrenNode *Body;                // Declare pointer to a class.

     // Save the address of first object
      Body = &Bodyadd;
 
      Body->readdata(); 
      Body->DisplayContents();  
 
      double maxd = 10.0;  // maximum diameter 

      //read link information 
      ifstream rlink;
      rlink.open(constants::flink1); 
      rlink >> links;   // number of links 
     
      //Link Throat;  
      OrenLink Throatadd, *Throat;
      Throat = &Throatadd; 
 
      Throat->readdata(); 
      Throat->DisplayContents();  

      double linkl, density;

      std::vector<double> den;   // density 
      den.reserve(nodes); 
      
      std::vector<int> ni;       // # of connections 
      ni.reserve(nodes);
     
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
      
      double sumcond,tmp,tmpi,tmpj,tmpl,tmps;
      double rho,cnstc, pavg;
      double len = 100.0;         // link length 
      int k1, k2,ic,Id;
      
      int maxIters=90000;    
      double resid;
      int iters, nsec=nx*ny;
      
      double length, area;
      double massin=0.0,massout=0.0,tmass; 
      double vcomp,effr, voli;
      double rid, rjd, qid, qjd;  
      double porevol = 0, klef, ckkb = 1e7;
      double permeability, flowr, flow;  
      double l2norm; 
      double domainvol = Body->xc[nodes]*Body->yc[nodes]*Body->zc[nodes];  

      ofstream tout ("output.txt");


      //find min, max, avg coordination numbers 
      int minBodyRad = 100, maxBodyRad=-100, sumBodyRad=0.0; 
      for (int i=0; i<nodes; ++i) {
          if(minBodyRad > Body->rad[i]) minBodyRad = Body->rad[i]; 
          if(maxBodyRad < Body->rad[i]) maxBodyRad = Body->rad[i]; 
          sumBodyRad = sumBodyRad + Body->rad[i];
      }   

//    Declare memory block of size nodes
      double avgcord = 0; 
      int** conn = new int*[nodes];
      for (int i=0; i<nodes; ++i) { 
  	      conn[i] = new int[Body->conum[i]];
              avgcord = avgcord + Body->conum[i]; 
      } 

      //find min, max, avg coordination numbers 
      int minCord = 100, maxCord=-100;
      for (int i=0; i<nodes; ++i) {
          if(minCord > Body->conum[i]) minCord = Body->conum[i]; 
          if(maxCord < Body->conum[i]) maxCord = Body->conum[i]; 
      } 


      double minThroatLen = 100, maxThroatLen= -100, sumThroatLen=0.0; 
      double minThroatRad = 100, maxThroatRad= -100, sumThroatRad=0.0; 
      for (int i=0; i<links; ++i) {
          if(minThroatRad > Throat->rad[i]) minThroatRad = Throat->rad[i]; 
          if(maxThroatRad < Throat->rad[i]) maxThroatRad = Throat->rad[i];
          sumThroatRad= sumThroatRad + Throat->rad[i]; 

          if(minThroatLen > Throat->len[i]) minThroatLen = Throat->len[i]; 
          if(maxThroatLen < Throat->len[i]) maxThroatLen = Throat->len[i];
          sumThroatLen= sumThroatLen + Throat->len[i]; 
      } 

      // initialize rhs, pressure and density 
      for(int i=0; i<nodes; i++) {
          rhs.push_back(0.0);
          den.push_back(1.0);
          pressures.push_back(poutlet);
      }

      area = (xt/nonle)*(yt/nonle); //PI*(yt/2.0)*(yt/2.0); 
      length = zt/nonle; 
      ofstream myfile ("pressure_out.txt");
      // constant in adsorption term  
      cnstc = 1.0/(pmax-pmin)*(0.1*pmax/(1.0+0.1*pmax) - 0.1*pmin/(1.0+0.1*pmin));

      pold = pressures ;
      //loop through the nodes and make connectivity array 
      for(int i=0; i<nodes; i++)
    	  ni[i] = 0;


      for(int j=0; j<links; j++) {
    	  k1 = Throat->ends.at(j).first;
    	  k2 = Throat->ends.at(j).second;

    	  if(k1 > 0) { //node1 ==0 : inlet or outlet
    		  i1 = ni[k1-1];
    		  conn[k1-1][i1] = j ;
    		  ni[k1-1] = ni[k1-1]+ 1;
    	  }
    	  else if(k2> 0) {
        	   i2 = ni[k2-1];
        	   conn[k2-1][i2] = j ;
        	   ni[k2-1] = ni[k2-1]+ 1;
    	  }

      }

      // calculate initial density 
      for(int i=0; i<nodes; i++) {
    	  den[i] = calrho(pressures[i],pnond,tempe);
      }
      // calculate porosity
      double porosity = 0.0;
      for(int i=0; i<nodes; i++)
    	  porosity = porosity + (4.0/3.0)*PI*pow(Body->rad[i],3);       // add nodes volume
      for(int i=0; i<links; i++)
    	  porosity = porosity + PI*pow(Throat->rad[i],2)*Throat->len[i]; // add links volume

      porosity =  porosity/(xt*yt*zt);

      tout << minBodyRad  <<"\t"<< maxBodyRad <<"\t"<< sumBodyRad/float(nodes) << "\t" ;
      tout << minCord <<"\t"<< maxCord <<"\t"<< avgcord/float(nodes) << "\t"; 
      tout << minThroatLen <<"\t"<< maxThroatLen  <<"\t"<< sumThroatLen/float(nodes) << "\t"; 
      tout << minThroatRad <<"\t"<< maxThroatRad  <<"\t"<< sumThroatRad/float(nodes) << "\t"; 
      tout << porosity << endl;
 
      /*******************************************************************************/
      // loop through time steps 
      for(int k=1; k<nsteps; k++) {
          row.clear();
          col.clear();
          rhs.clear();
          val.clear();

      #pragma omp parallel
          {
          //loop through the nodes
	  std::vector<int> row_th; 
	  std::vector<int> col_th; 
	  std::vector<double> val_th; 
	  std::vector<double> rhs_th; 
          #pragma omp for nowait schedule(static)
    	  for(int i=0; i<nodes; i++) {
    		  sumcond   = 0.0;
    		  rhs_th.push_back(0.0);
    		  nflow[i]  = 0.0;
    		  voli      = calvol(Body->rad[i]);

    		  if(ni[i] ==0) {  //number of connections  = 0
    			  val_th.push_back(1.0);
    			  row_th.push_back(i);
    			  col_th.push_back(i);
    			  rhs_th[i] = poutlet;
    		  }
		     else { 
	   		     tmps = 0.0; 
	   		     for(int j=0; j<ni[i]; j++) {
			    	     Id  = conn[i][j];
				     k1  = Throat->ends.at(Id).first  -1;    // node at one end of the link
				     k2  = Throat->ends.at(Id).second -1;     // node on the other end of the link
				     linkl = Throat->len[Id];  // link length
				     if(Throat->id[Id] == 1) { // connected to inlet
				    	 /*
				    	 //if((k1+1) > 0) tmpl  = conductance(Throat.rad[Id], mu,linkl)*((pinlet+pressures[k1])/(2.0*pressures[k1])) ; // compressible
				    	 //if((k2+1) > 0)  tmpl  = conductance(Throat.rad[Id], mu,linkl)*((pinlet+pressures[k2])/(2.0*pressures[k2])) ; // compressible\
				    	 //klef = (1.0+2.0*(ckkb/pnond)/(pressures[k1]+pressures[k2]));
				    	 tmpl  = conductance(Throat.rad[Id], mu,linkl);
				    	 tmp   = tmpl*dt; //1.0/(1.0/tmpi + 1.0/tmpl + 1.0/tmpj);
				    	 sumcond = sumcond + tmp; //*klef;
				    	 rhs[i] = rhs[i] - pinlet*tmp;
				    	 */ 
		    			  val_th.push_back(1.0);
		    			  row_th.push_back(i);
		    			  col_th.push_back(i);
		    			  rhs_th[i] = pinlet;
				     }
				     else if (Throat->id[Id] == 2) {
				    	 /*
				    	 //if((k1+1) >0 ) tmpl  = conductance(Throat.rad[Id], mu,linkl)*((poutlet+pressures[k1])/(2.0*pressures[k1])) ;
				    	 //if((k2+1) >0 ) tmpl  = conductance(Throat.rad[Id], mu,linkl)*((poutlet+pressures[k2])/(2.0*pressures[k2])) ;
                                         
				    	 tmpl  = conductance(Throat.rad[Id], mu,linkl);
				    	 tmp   = tmpl*dt;
				    	 sumcond = sumcond + tmp;
				    	 rhs[i]  = rhs[i] - poutlet*tmp;
                                         */
		 			  val_th.push_back(1.0);
		    			  row_th.push_back(i);
		    			  col_th.push_back(i);
		    			  rhs_th[i] = poutlet;

				     }
				     else  {
				    	 //only pore throats offer resistence to flow but they do not have volume
				    	 //pore bodies possess volumes but they do not provide resistance to flow

				    	  tmpl  = conductance(Throat->rad[Id], mu,linkl);  
                                          if(flowType >1) tmpl *= ((pressures[k1]+pressures[k2])/(2.0*pressures[i])) ; // compressible
					  tmp   = tmpl*dt;
					  sumcond = sumcond + tmp; //*klef;
  					     if(i==k1) {   // k1 -node number  >0  internal node
					    	 val_th.push_back(tmp);
					         row_th.push_back(i); 
					         col_th.push_back(k2);
					     }
					     else if (i==k2) { // k2 -node number  > 0 internal node
					    	 val_th.push_back(tmp);
						 row_th.push_back(i); 
						 col_th.push_back(k1); 
					     }
				     }
                  
			     }
			     if(voli > 0.0 && flowType >1) vcomp  = voli/pnond;          // compressibility term
			     if(sumcond < 1.e-200  && sumcond > -1.e-200) sumcond = 1.0;
	                     if(isnan(sumcond)) sumcond = 1.0;
                             if(flowType >1) {
			        rhs_th[i] = rhs_th[i]- vcomp ; //modify right hand side vector to consider sorption and compressibility
				val_th.push_back(-sumcond-vcomp/pressures[i]); // add compressibility term
                             }
			     else { 
				val_th.push_back(-sumcond); // add compressibility term
			     }
			     row_th.push_back(i);
			     col_th.push_back(i);
	          }
       	  }

          #pragma omp for schedule(static) ordered
          for(int i=0; i<omp_get_num_threads(); i++) {
          	#pragma omp ordered
           	row.insert(row.end(), row_th.begin(), row_th.end());
           	col.insert(col.end(), col_th.begin(), col_th.end());
           	val.insert(val.end(), val_th.begin(), val_th.end());
           	rhs.insert(rhs.end(), rhs_th.begin(), rhs_th.end());
          }
          }

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

	  // calculate massflow at the inlet 
	  // massin = 0.0; 
    	  flowr = 0.0;
    	  for(int i=0; i<links; i++) {
    		  k1  = Throat->ends.at(i).first  ;
    		  k2  = Throat->ends.at(i).second ;
    		  if(Throat->id[i] == 1) {
    			  if(k1 < 0) i1 = k2-1;
    			  if(k2 < 0) i1 = k1-1;
    			  linkl  = Throat->len[i];
    			  tmpl   = conductance(Throat->rad[i], mu,linkl); 
                          if(flowType>1) tmpl *= ((pinlet+pressures[i1])/(2.0*pressures[i1])) ;
    			  tmp    = tmpl*dt;
    			  density = calrho(pinlet,pnond,tempe);
    			  //massin  = massin + tmp  * density*(pinlet - pressures[i1]);
    			  flowr   = flowr  + tmpl * (pinlet - pressures[i1]);
    		  }
    	  }
    	  //permeability =  flowr*nonle*nonle*mu*length/(area*nonle*(pinlet-poutlet)*pnond);
    	  permeability = nonle*nonle*flowr*mu*length/(area*(pinlet-poutlet));
    	  cout << "permeability "<<flowr<<"\t"<<nonle<<"\t"<< mu <<"\t"<< length <<"\t"<< area <<"\t"<< pinlet <<"\t"<< permeability <<"\t" << permeability/(0.987*1e-12) <<"\t"<< endl;   ///(0.987*1e-18) << endl;
    	  tout << permeability <<"\t" << permeability/(0.987*1e-12) <<"\t"<< endl;   /// <<  endl;
      } // time step loop 
      
     
      for(int i=0; i<nodes; i++) {
	  myfile << pressures[i] <<" "<<Body->xc[i] <<" "<<Body->yc[i]<<" "<<Body->zc[i]<<" "<<  std::endl;
      }
      myfile.close();
     
      writeData out2file; 
      out2file.NodesVTK(nodes, pressures); 
      out2file.LinksVTK(nodes, links, pressures); 

      cout << "print vtp files: completed" << endl; 
      //system("gnuplot -p -e plot_mass.gnu");

      free(nreff);
      free(nflow);
      free(conn);
      
      ierr = PetscFinalize();
      if (ierr) {return ierr;}
}
