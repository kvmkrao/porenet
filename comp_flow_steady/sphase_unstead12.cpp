
#include<cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <sstream>
#include <vector>

// petsc header files 
#include <petscksp.h>
#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscviewer.h>

// include linear solver header file 
#include "linear_solver_petsc.h"

using namespace std;
const double PI  =3.141592653589793238463;


void readnode(double* xn, double* yn, double* zn, double* dia) {
	int nn;
	ifstream nodeloc("mnodes.dat");
	nodeloc >> nn;
	
	ifstream nodedia("mdiameter.dat");
	for(int i=0; i<nn ; i++) {
	     nodeloc >> xn[i] >> yn[i]>> zn[i];	
	     nodedia >> dia[i]; // in  nm 
	    // dia[i] = dia[i]; 
	}
	
	nodeloc.close();
	nodedia.close(); 
	return;
}

void readlink(int* l1, int* l2, double* ld) {
	int nl;
	ifstream linknn("link.dat");
	linknn >> nl;
	
	for(int i=0; i<nl ; i++) {
	     linknn >> l1[i] >> l2[i];
	     ld[i] = 100.0; ////1e-7;
	}
	
	linknn.close();
	return;
}

double conductance( double diameter, double vis, double len) 
{
	double effdia = sqrt(4.0/PI)*diameter; 
	//double effdia = diameter; 
	return (PI*pow(effdia/2.0,4)/(8.0*vis*len));
	//return (PI*pow(diameter/2.0,4)/(8.0*mu*ped))
}

double dv_sorption( double diameter, double cnstc, double rho) 
{
	double effdia = sqrt(4.0/PI)*diameter; 
	//double effdia = diameter; 
	double cmg    = 2.0e-7; 
	return (6.0*pow(effdia,2)*cmg*cnstc/rho);
}

void write_vtp(int nodes, int links, double* xc, double* yc, double* zc, int* l1, int* l2, double* nd, std::vector<double> &pressures) 
{
       ofstream output ("single_pore_out.vtp");
       output <<"<VTKFile byte_order=\"LittleEndian\" type=\"PolyData\" version=\"0.1\">" << endl;
       output <<"        <PolyData>"<< endl;
       output <<"            <Piece NumberOfLines=\"" << links <<"\" NumberOfPoints=\""<< nodes <<"\">" << endl;
       
       
       output <<"                <Points>" << endl;
       output <<"                   <DataArray Name=\"coords\" NumberOfComponents=\"3\" type=\"Float64\">";
       for(int i=0; i<nodes ; i++) {
               output << xc[i]<<" \t "<< yc[i]<<" \t "<< zc[i] <<"\t";
       }
       output <<"</DataArray>" << endl;
       output <<"                </Points>" << endl;
       
       
       output <<"                <Lines>" << endl;
       output <<"                <DataArray Name=\"connectivity\" NumberOfComponents=\"1\" type=\"Int64\">";
       for(int i=0; i<links ; i++) {
       	output  <<  l1[i] <<" \t " << l2[i] << "\t";
       	//cout    <<  l1[i] <<" \t " << l2[i] << endl;  
       }
       
       output <<"</DataArray>" << endl;
       output <<"                </Lines>" << endl;
       
       output <<"                <PointData>" << endl;
       output <<"                   <DataArray Name=\"pore.diameter\" NumberOfComponents=\"1\" type=\"Float64\">";
       for(int i=0; i<nodes ; i++) {
               output << nd[i] <<" \t " ;
       }
       output <<"</DataArray>" << endl;
       
       output <<"                   <DataArray Name=\"pore.pressure\" NumberOfComponents=\"1\" type=\"Float64\">";
       for(int i=0; i<nodes ; i++) {
               output << pressures[i] <<" \t " ;
       }
       output <<"</DataArray>" << endl;
       output <<"                </PointData>"<<endl;
       output <<"            </Piece>" << endl;
       output <<"        </PolyData>"<<endl;
       output <<"</VTKFile>" << endl;
       
       output.close();
       return ; 
}


int main(int argc, char *argv[])
{
      static char help[] = "Solves a linear system in parallel with KSP.\n";
      
      PetscErrorCode ierr;
      ierr = PetscInitialize(&argc,&argv,(char*)0,help);
      
      double mu1, mu2, r, L, sigma;
      double small;
      int n,nodes,links;
      double xt, yt, zt, p1, p2,A,dp,pc;
      double pinlet = 1.0 ; //101325*10;    //           5 atm 
      double poutlet= 0.5; //101325*5;   //101325;  //10 atm 
      double dt     = 0.5e-10; 
      double zleft  = 1.0   ; //0.8e-7;  
      double zright = 960.0  ; // 9.6e-7;  
      
      double pmin = 0.5, pmax = 1.0;  //MPa 
      //const double PI  =3.141592653589793238463;
      
      ifstream rnode("mnodes.dat");
      rnode >> nodes; 
      
      double* xc   = (double*)malloc(nodes * sizeof(double*));
      double* yc   = (double*)malloc(nodes * sizeof(double*));
      double* zc   = (double*)malloc(nodes * sizeof(double*));
      double* nd   = (double*)malloc(nodes * sizeof(double*));
      
      // read node locations and diameters 
      readnode(xc, yc, zc, nd); 
      
      /*
      // check nodes 
      for(int i=0; i<nodes ; i++) {
      //	  rnode>> xc[i] >> yc[i]>> zc[i]; 
      //        rnodedia >> nd[i];
                cout << i<<"\t"<< xc[i] <<"\t" << yc[i] <<"\t"<< zc[i] <<"\t"<< nd[i] << endl;
                if(i >1000) return 0; 	
      }
     */

      //read link information 
      ifstream rlink("link.dat");
      rlink >> links;
      
      int* l1     = (int*)malloc(links * sizeof(int*));   // node at one end of link 
      int* l2     = (int*)malloc(links * sizeof(int*));   // node at the other end of the link 
//      int* lp     = (int*)malloc(links * sizeof(int*)); 
      double*  ld = (double*)malloc(links * sizeof(double*));
//      double*  ll = (double*)malloc(links * sizeof(double*));
//      double* n1l = (double*)malloc(links * sizeof(double*));
//      double* n2l = (double*)malloc(links * sizeof(double*));
//      double*  la = (double*)malloc(links * sizeof(double*));
      
      int*     wet    = (int*)malloc(nodes * sizeof(int*));
      double*  sthres = (double*)malloc(nodes * sizeof(double*));
      double*  disp   = (double*)malloc(nodes * sizeof(double*));
      
      std::vector<double> den;  // density 
      den.reserve(nodes); 
      
      std::vector<double> sat;  // saturation 
      sat.reserve(nodes);
      
      std::vector<int> ni;   // # of connections 
      ni.reserve(nodes);
      
      std::vector<double> rates; // flow rate
      rates.reserve(nodes);
      
      std::vector<double> val;   // non-zero val of matrix 
      val.reserve(links);
      
      std::vector<double> rhs;  // rhs of matrix 
      rhs.reserve(nodes);
      
      std::vector<double> pressures; // unknown vector 
      pressures.reserve(nodes); 
      
      std::vector<double> pold;    // unknown vector 
      pold.reserve(nodes);
      
      std::vector<int> row;       // row index of matrix 
      row.reserve(links); 
      
      std::vector<int> col;       // col index of matrix 
      col.reserve(links); 
      
      double sumcond,tmp,tmpi,tmpj,tmpl,tmps,position,vol,shapeFactor,clay,distance;
      double rho,cnstc; 
      double mu  = 1.5e-5/1e-9;         // viscosity  (scaled nm) 
      double len = 100; //1e-7;           // link length 
      int k1, k2,ic,id;
      
      int maxIters=90000;    
      double resid;
      int iters;
      
      double flow=0.0;
      double length, area;
      double domainvol = xc[nodes]*yc[nodes]*zc[nodes]; 
      double massin=0.0,massout=0.0,tmass; 
      double rpart1, rpart2,effr, voli;
      double ped, pedk1, pedk2; 
      double porevol = 0, klef, ckkb = 1e7; 
      double pnond   =1013250; 
       
      int conn[nodes][10];
      
      //read link information 
      readlink(l1, l2, ld);
      
      /*
      for(int i=0; i<links ; i++) {
      //         rlink >> l1[i] >> l2[i];
      //        rlinkdia >> ld[i];
      //        ld[i] = 1e-7; 
                cout << "link" << i <<"\t "<< l1[i] <<"\t "<< l2[i] <<"\t "<< ld[i] << endl;  
      }
      */

      for(int i=0; i<nodes; i++) {
         rhs.push_back(0);
         pressures.push_back(poutlet); // initial pressure 
         sat.push_back(1e-3);          // initial saturation 
      }
      
      // apply pressure bc 
      for(int i=0; i<nodes; i++) {
	 disp[i] = 0.0; 
         if(zc[i] < zleft)  {
		 pressures[i] = pinlet; 
		 disp[i] = 0.1; //0.1e-9;
	 } 
	 else if(zc[i] > zright) pressures[i] = poutlet;
          // calulcate density at nodes  
         den[i] =  pressures[i]*pnond/(8.31447*298.15)*0.0291;
         //den[i] =  pressures[i]/(8.31447*298.15)*0.0291;
         // initial displacement of the fluid in x-direction; 
	 //cout << "node "<< i <<"\t "<< den[i] << endl; 
      }
      
     
      pold = pressures ;
      //cout << nodes << "\t " << links << endl;   
      //loop through the nodes and make connectivity array 
      for(int i=0; i<nodes; i++) {
        ic    = 0;
        ni[i] = 0;
        // loop through the links
        for(int j=0; j<links; j++) {
            k1 = l1[j];
            k2 = l2[j];
            //printf("%d %d %d\n", k1, k2,i);
            // check the link end nodes with the node
            if(i==k1 || i==k2) {
                conn[i][ic] = j;
                ic = ic + 1;
      	        if(i==1) cout << j <<"\t"<< k1 <<"\t"<< k2 <<endl; 
            }
        }
        ni[i] = ic; 
        //if(i==1) cout << "shape factor i=1 "<< ni[i] <<endl;
        //cout <<" node shape factor "<< i <<"  "<< ni[i] << endl;
      }

  
      area   = 3.120000e+03 ;  // not in use 
      length = 3.120000e+03 ;  // not in use 
      
      ofstream myfile ("pressure_out.txt");
      ofstream outmass ("time_massflow.txt");
      ofstream dispfile ("displacement.txt");
      
      
      // loop to define hydrocarbon/water wet 
      for(int i=0; i<nodes; i++) {
      	wet[i] = 1;                 // water - wet  
      	if(xc[i] >200)   wet[i] =2;  // xc index  = 3- 7 (yz plane) 
      	if(xc[i] <800)   wet[i] =2;  //
        //if(yc[i] == 5e-7) wet[i] =1; 	
        //if(yc[i] == 6e-7) wet[i] =1; 	
      }
      // loop to define hydrocarbon wet 
      
      //calculate saturation threshold 
      for(int i=0; i<nodes; i++) {
              tmp  = 0;
      	
              if(ni[i] ==0 ) {  //no connections - isolated
                 sthres[i] = 0.0; 
              }
              else if(zc[i] <=zleft) {
                 sthres[i] = 0.0; 
              }
              else if(zc[i] >= zright) {
                 sthres[i] = 0.0; 
              }
              else
      	      {
      	         for(int j=0; j<ni[i]; j++) {
      	           id  = conn[i][j];   // link index 
                   tmp = tmp  + 0.4764*(pow((ld[id]-nd[id])/nd[i],3));
      	           //cout << " tmp" << nd[i] << ld[]	
                 }
      	         sthres[i] = 0.4764 ; //abs(tmp)/float(ni[i]+1); 
      	         //cout <<" simb threshold"<< i <<" \t " << sthres[i] << endl; 
      	      }
      }
     
//      exit(0); 

      cnstc = 1.0/(pmax-pmin)*(0.1*pmax/(1.0+0.1*pmax) - 0.1*pmin/(1.0+0.1*pmin));


      for(int i=0; i<nodes; i++) {
	      for(int j=0; j<ni[i]; j++) {
                       id      = conn[i][j];
                       k1      = l1[id];
                       k2      = l2[id];
                       tmp     = (pressures[k1]+ pressures[k2])/2.0; 
                       //cout << "node " << i <<"pressure " << pressures[k1] <<"\t"<< pressures[k2]<<"\t"<< tmp <<"\t"<< k1 <<"\t"<<k2 <<  endl;
	      }
      }
      
      //exit(0); 

      // loop through time steps 
      for(int k=1; k<100; k++) {
      printf("time step loop %i\n", k);

/*
        //loop through nodes 
        for(int i=0; i<nodes; i++) {
              voli      = (4.0/3.0)*PI*pow(nd[i]/2.0,3);
              sumcond   = 0.0;
              rates[i]  = 0.0;
	      den[i]    = pressures[i]*pnond/(8.31447*298.15)*0.0291;
      	      rho       = den[i]; //pressures[i]/(8.31447*298.15)*0.0291; //pressures[i]/(8.31447*298.15)*0.0291;  //kg/mol
      	      tmps      = 0.0; 
      	      ped       = nd[i]; 
              //if(i==1) cout <<" node shape factor " << i  <<"  "<<ni[i] << endl;   
              // loop through the links
              // isolated           left bc        right bc 
              if( ni[i] !=0 || (zc[i] > zleft) || (zc[i] < zright) ) {
                   if(wet[i] != 1) {
                       // gas sorption - modification to effective pore radius 
                       effr    = 2.0*0.4e-9*0.1*pressures[i]*pnond/(1.0+0.1*pnond*pressures[i]);
                       ped     = nd[i] - effr;
                       tmps    = dv_sorption(ped, cnstc, rho);
		      // cout << i << "\t" << ni[i] <<"\t"<< wet[i] <<"\t " << zc[i] << endl; 
	//	       cout << i << "\t" << wet[i] <<"\t " << zc[i] << endl; 
                   }
                 //  cout <<" "<< i <<" "<< ni[i] <<" "<< rho << endl;   
                   for(int j=0; j<ni[i]; j++) {
                       id      = conn[i][j];
                       k1      = l1[id];
                       k2      = l2[id];
                       //if(i==1) cout << i <<" "<<id<<" "<< k1 <<" "<< k2 << endl; 
                       ld[id]= 100.0 - nd[k1]/2.0- nd[k2]/2.0;
                       pedk1   = nd[k1];
                       pedk2   = nd[k2];
      
                       tmpi    = conductance(pedk1, mu, pedk1);
                       tmpj    = conductance(pedk2, mu, pedk2);
                       tmpl    = conductance(ld[id], mu,ld[id])*((pressures[k1]+pressures[k2])/(2.0*pressures[i])) ; // compressible 

		       // conductance 
      		       tmp     = 1.0/(1.0/tmpi  + 1.0/tmpl + 1.0/tmpj);
		      // cout << " node " << i <<"\t "<< ld[id]<<"\t "<< mu << endl; 

		       // klinkenberg factor 
		       klef  = (1.0 + 2.0*(ckkb/pnond)/(pressures[k1]+pressures[k2]));
		       //modify the conductance  
                       sumcond = sumcond  + tmp *klef;

                       massin  = massin   + tmp*(pressures[i]-pinlet)*pnond;
                       massout = massout  + tmp*(pressures[i]-poutlet)*pnond;
                       rates[i]= rates[i] + tmp*(pressures[i]-pressures[k1])*pnond;
                       rates[i]= rates[i] + tmp*(pressures[i]-pressures[k2])*pnond;
               }
               //update saturation 
               rpart1  = dt*(rates[i]) ; // + tmps*(pold[i] - pressures[i]));
               rpart2  = voli*sat[i];
               // update density based on saturation and pressures; 
              // den[i]  = den[i]  * (rpart1 + rpart2)/rpart2 ;
               sat[i]  = sat[i]  + (dt/voli)*abs(rates[i]);
               disp[i] = disp[i] +  dt*abs(rates[i])/(PI*pow(nd[i]/2.0,2.0));
               dispfile << k <<"\t" << i <<"\t "<< disp[i]<<"\t "<< abs(rates[i])/(PI*pow(nd[i]/2.0,2.0))  <<"\t" <<rates[i]<< endl;
               //if(sat[i] > 0.001) cout << "saturation at " << i << " is " << sat[i] << " \t" << sthres[i] <<"\t" <<rates[i] << endl;
              }  // if connection condition  
        }   //for node loop 
     
        //cout << " begin"	 << endl; 
        //exit(0); 

*/
        //loop through the nodes 
        for(int i=0; i<nodes; i++) {
              sumcond = 0;
              rhs[i]  = 0;
      	      voli    = (4.0/3.0)*PI*pow(nd[i]/2.0,3);
	      //cout <<" node i" << i <<"\t"<< pressures[i] << "\t " << poutlet << endl; 
	      den[i]  = pressures[i]/(8.31447*308.15)*0.0291; //pressures[i]/(8.31447*298.15)*0.0291;  //kg/mol
              //if(i==1) cout <<" node shape factor " << i  <<"  "<<ni[i] << endl;   
              // loop through the links
	      tmp = (pressures[k1] + pressures[k2])/2.0; 
              
              cout << rho <<"\t"<< cnstc <<"\t"<<tmp <<"\t"<< pressures[i] << endl; 
              if(ni[i] ==0) {  //shape factor = 0 
                  val.push_back(1.0);
                  row.push_back(i);
                  col.push_back(i);
      	          rhs[i] = poutlet; 
              }
              else if(zc[i] <zleft) { 
        	  rhs[i] = pinlet; 
                  val.push_back(1.0);
                  row.push_back(i);
                  col.push_back(i);
              }
              else if(zc[i] > zright) {
                  rhs[i] = poutlet;
                  val.push_back(1.0);
                  row.push_back(i);
                  col.push_back(i);
              }
	      else {
              //else if(xc[i]<= disp[i] && disp[i] <= xright) {
              //else if(pressures[i] > poutlet) {
              //else if(tmp > poutlet) {
//	      else if(abs(sat[i]) >= sthres[i]) { //    xc[i]<= disp[i] && disp[i] <= xright) {
        	  rho  = den[i];  //pressures[i]/(8.31447*298.15)*0.0291;  //kg/mol
        	  //rho  = pressures[i]/(8.31447*298.15)*0.0291;  //kg/mol
      	          //cout << " density" << rho << endl;  
      	          ped  = nd[i];
      	          tmps = 0.0;  
		  /*
      	          if(wet[i] != 1) {
      		      // gas sorption - modification to effective pore radius 
      		      effr = 2.0*0.4e-9*0.1*pressures[i]*pnond/(1+0.1*pressures[i]*pnond); 
                      ped  = nd[i] - effr; 
      		      //  cout << "effective radius" << effr <<"\t" << nd[i] << endl;  
                      //  cout << i << "radius" << nd[i] <<"\t "<< ped << endl;  
                      tmps = dv_sorption(ped, cnstc, rho);
      	              //tmps  = 6.0*pow(ped,2)*2.0e-7*cnstc/rho;
      	          }
              */	    

                  for(int j=0; j<ni[i]; j++) {
                      id  = conn[i][j];
                      k1  = l1[id];
                      k2  = l2[id];
		      tmp = (pressures[k1] + pressures[k2])/2.0; 
                      //if(tmp > poutlet) 
		      {
			   cout << "node  " <<"\t" << tmp << "\t "<< k1 << "\t" << k2 << endl; 
        	         //if(i==1) cout << i <<" "<<id<<" "<< k1 <<" "<< k2 << endl; 
        	           ld[id] = 100.0 - nd[k1]/2.0- nd[k2]/2.0; 
      		           pedk1 = nd[k1]; 
      		           pedk2 = nd[k2]; 
        	              //ld[id] = 1e-7 - nd[k1]/2.0- nd[k2]/2.0 - effr;
      
                           tmpi = conductance(pedk1, mu, pedk1);
                           tmpj = conductance(pedk2, mu, pedk2);
        
                           // link (ij condutance)
                           //tmpl = ((pressures[k1]+pressures[k2])/(2.0*pressures[i]))* PI*pow((ld[id]/2.0),4)/(8.0*mu*(ld[id])); // compressible 
                           //tmpl = conductance(ld[id], mu,ld[id])*((pressures[k1]+pressures[k2])/(2.0*pressures[i])) ; // compressible 
                           //tmpl = PI*pow((ld[id]/2.0),4)/(8.0*mu*(ld[id]));                                                 //  incompressible 
                            
                           tmpl = conductance(ld[id], mu,ld[id])*((pressures[k1]+pressures[k2])/(2.0*pressures[i])) ; // compressible 
                           tmp  = 1.0/(1.0/tmpi + 1.0/tmpl + 1.0/tmpj);
		           klef = (1.0+2.0*(ckkb/pnond)/(pressures[k1]+pressures[k2])); 

                           sumcond = sumcond + tmp; //*klef;
          	             
        	           // Applying pressure inlet  condition 
                           if(zc[k1] < zleft || zc[k2] < zleft) {
                                   //  if any end connected to inlet 
                                   rhs[i] = rhs[i] - pinlet*tmp;
                                   //cout  << xc[i] <<"\t"<< pinlet << endl; 
                           }
                            // Applying pressure outlet condition 
        	           else if(zc[k1] > zright || zc[k2] > zright) {
                                   // if any end connected to outlet 
                                   rhs[i] = rhs[i] - poutlet*tmp;
                                   //cout  << xc[i] <<"\t"<< poutlet << endl; 
                           }
        	           else {
      		                   if(i==k1) {   // k1 -node number  >0  internal node 
                                      val.push_back(tmp);
                                      row.push_back(i); 
                                      col.push_back(k2);
                              //      if(i==1) cout << i <<"\t"<<id<<"\t"<< k1 <<"\t"<< k2 <<"\t" << tmp << endl;
                                   }
                                   else if (i==k2) { // k2 - node number  > 0 
                                      val.push_back(tmp);
                                      row.push_back(i);
                                      col.push_back(k1);
                                      // if(i==1) cout << i <<"\t"<<id<<"\t"<< k2 <<"\t"<< k1 <<"\t" << tmp << endl;
                                   }
                           }
      		      //rates[i]  = rates[i] + tmp*(pressures[i]-pressures[k1]);
      		      //rates[i]  = rates[i] + tmp*(pressures[i]-pressures[k2]);
		      } // if condition 
                  }  // connection loop

                  // cout <<" sum cond "<< sumcond <<" "<< tmps/dt<<" "<< tmps <<" "<<dt << endl; 
                  //if(i%2 ==0)
      	          rpart1  = voli*sat[i];
      /*
      	          if(wet[i] !=1) {
      		         sumcond  = sumcond + tmps/dt + rpart1/pressures[i]/dt;
      	          }
      	          else {
      		         sumcond  = sumcond + rpart1/pressures[i]/dt;
      	          }
     */

      	          //rates[i] = sumcond ; 
                  if(sumcond < 1.e-200 && sumcond > -1.e-200) sumcond = 1.0;
                  val.push_back(-sumcond);
        //          rhs[i]          = rhs[i]-pressures[i]*(tmps/dt) - rpart1/dt ;
                  row.push_back(i);
                  col.push_back(i);
                  cout << i <<"\t"<<id<<"\t"<< i <<"\t"<< ni[i] <<"\t" << sumcond << endl;
      

              }  // if condition loop end 
      	      /*
      	      else if (disp[i] > xright) {
      		  cout << "flow out of the outlet"<< endl; 
      		  exit(0); 
      	      }
      	      
      	      else {
      	          row.push_back(i); 
      	          col.push_back(i); 
      	          val.push_back(1.0); 
      	          rhs[i] = pressures[i]; 
      	      }
              */ 
              //}  // if connection condition  
        }   //for node loop 
      
	//exit(0); 
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
      	
      //printf("maxiter, res %d %f\n",iters, resid); 
      
        cout << "time step "<< k <<" "<<k*dt<<" "<< massin*dt << " "<<massout*dt << endl; 
      
        cout << "avg flow rate  "<< flow <<" "<< flow*mu*length/area/(pinlet-poutlet) << std::endl;
        outmass   << k <<"\t"<< k*dt <<"\t" << massin*dt <<"\t"<< -massout*dt <<"\t" << flow <<"\t"<< flow*mu*length/area/(pinlet-poutlet) << std::endl;
      } // time step loop 
      //return 0; 
      
      //std::cout << "avg flow rate  "<< flow <<" "<< flow*mu*length/area/(pinlet-poutlet) << std::endl;
      double porosity     =  domainvol/porevol;
      cout << " porosity "<< porosity << std::endl;
      double permeability =  flow*mu*length/area/(pinlet-poutlet);
      cout << "avg flow rate  "<< flow <<" "<< permeability/(0.987*1e-12) << std::endl;
      //1 darcy = 0.987 micro m^2
      
      
      for(int i=0; i<nodes; i++) {
          myfile << pressures[i] <<" "<<xc[i] <<" "<<yc[i]<<" "<<zc[i]<<" "<< sat[i] << std::endl;
      //    cout   << pressures[i] <<" "<<xc[i] <<" "<<yc[i]<<" "<<zc[i]<<" "<<std::endl; 
      }
      myfile.close();
      
      write_vtp(nodes, ic, xc, yc, zc, l1, l2, nd, pressures); 

      ofstream outavgpfile ("z_avg_pressure_out.txt");
      // calculate average pressure
      double zpresavg[nodes/121];
      int indx;
      for(int i=0; i<nodes/121; i++) {
              zpresavg[i] = 0.0;
              for(int j=0; j<121; j++) {
                      indx = i*121+j;
                      zpresavg[i] = zpresavg[i]+ pressures[indx];
              }
              outavgpfile << zc[i*121+1]<<" "<< zpresavg[i]/121.0 << std::endl;
      }
      outavgpfile.close();

      free(xc);
      free(yc);
      free(zc);
      free(nd);
      free(wet); 
      free(sthres); 
      free(disp); 
      
      free(l1);
      free(l2);
      free(ld);
//      free(ll);
//      free(la);
//      free(n1l);
//      free(n2l);
      
      ierr = PetscFinalize();
      if (ierr) {return ierr;}
} 
