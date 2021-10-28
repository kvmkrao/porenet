
#include<cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <sstream>
#include <vector>

#include <petscksp.h>
#include <iostream>

#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscviewer.h>

#include <mpi.h>

#include "linear_solver_petsc.h"

using namespace std;


//int main() { 
int main(int argc, char *argv[])
{
static char help[] = "Solves a linear system in parallel with KSP.\n";

PetscErrorCode ierr;
ierr = PetscInitialize(&argc,&argv,(char*)0,help);
int mpi_rank, comm_size;

MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

int nodes,links;
double dd1, dd2, dd3, dt; 
double pinlet = 101325;    // 1atm 
double poutlet= 101325*5;  // 5atm 

const double PI  =3.141592653589793238463;

ifstream rnode ("node.dat");
rnode >> nodes >>  dd1 >>  dd2 >> dd3; 

//double* xc = new double[nodes];
double* xc = (double*)malloc(nodes * sizeof(double*));
//double* yc = new double[nodes]; 
double* yc = (double*)malloc(nodes * sizeof(double*));
double* zc = (double*)malloc(nodes * sizeof(double*));
//double* zc = new double[nodes]; 
//double* nr = new double[nodes]; 
double* nr = (double*)malloc(nodes * sizeof(double*));
//double* pres = new double[nodes]; 
double* pres = (double*)malloc(nodes * sizeof(double*));

ifstream rlink("link.dat");
rlink >> links;

//int* l1   = new int[links];
int* l1   = (int*)malloc(links * sizeof(int*));
//int* l2   = new int[links];
int* l2   = (int*)malloc(links * sizeof(int*));
int* lp   = (int*)malloc(links * sizeof(int*));
//double* lr  = new double[links];
double* lr = (double*)malloc(links * sizeof(double*));
//double* ll  = new double[links];
double* ll = (double*)malloc(links * sizeof(double*));
//double* n1l = new double[links];
double* n1l = (double*)malloc(links * sizeof(double*));
//double* n2l = new double[links];
double* n2l = (double*)malloc(links * sizeof(double*));

std::vector<double> val; 
val.reserve(links+1);

std::vector<int> shapf; 
shapf.reserve(nodes+1);

//std::vector<int> conn; 
//conn.reserve(nodes*17);

std::vector<double> rhs;
rhs.reserve(nodes+1);

std::vector<double> pressures;
pressures.reserve(nodes+1); 

std::vector<int> row;
row.reserve(links+1); 

std::vector<int> col;
col.reserve(links+1); 

double sumcond,tmp,tmpi,tmpj,tmpl,position,vol,shapeFactor,clay,distance; 
double mu  = 1.5e-5; 
double len = 1e-4; 
int k1, k2,ic,id; 

for(int i=0; i<nodes ; i++) {
	// units are in micro meters 
        rnode >> id >> xc[i] >> yc[i] >> zc[i] >> position >> shapf[i] >> nr[i] >> vol >> shapeFactor >> clay;
  //      cout <<  id <<" "<< shapf[id] << " "<< shapeFactor << endl; 
//        cout <<  id <<" "<< xc[id]<<" "<< yc[id]<<" "<< zc[id]<<" "<<position <<" "<< shapf[id] <<" "<<nr[id] << vol << shapeFactor << clay << endl; 
//	rnode >> id >> x >> y >> z >> position >> coordNumber >> radius >> vol >> shapeFactor >> clay;
}

for(int i=0; i<links ; i++) {
	// units are in micro meters 
        rlink >>    id >> l1[i] >> l2[i]  >>      lp[i]     >> lr[i]  >> vol >> shapeFactor >> clay >> n1l[i]      >>  n2l[i]      >> ll[i]   >> distance;
//	inStream >> id >> nodeInID >> nodeOutID >> position >> radius >> vol >> shapeFactor >> clay >> node1Length >>  node2Length >> length  >> distance;
        //rlink >>  id >> l1[i] >> l2[i]  >>      lp[i]        >> lr[i]  >>     vol >>       shapeFactor >> clay >>     n1l[i]      >>  n2l[i]      >> ll[i]   >> distance;
        cout  <<  id <<" "<< l1[i]<<" " << l2[i] <<"  "<< lp[i] <<"  "<< lr[i] <<" "<< n1l[i] << " "  <<  n2l[i]<<"  "<< ll[i]   << endl;
                //rlink >> nodeInID >> nodeOutID >> position >> radius >>    vol >>       shapeFactor >> clay >>     node1Length >>  node2Length >> length  >> distance;
}


//return 0;

for(int i=0; i<nodes; i++) {
   rhs.push_back(0);
   pressures.push_back(pinlet); 
}   


int conn[nodes+1][20];

//loop through the nodes
for(int i=0; i<nodes; i++) {
  ic = 0; 
  shapf[i] = 0;
  // loop through the links
  for(int j=0; j<links; j++) {
      k1 = l1[j];
      k2 = l2[j];
      //printf("%d %d %d\n", k1, k2,i);
      // check the link end nodes with the node 
      if(i+1==k1 || i+1==k2) {
          conn[i][ic] = j;
	  ic = ic + 1;  
      }
  }
  shapf[i] = ic;
  //cout <<" node shape factor " << i  <<"  "<<shapf[i] << endl;   
}

//return 0; 
//lp position array  (0- internal; 1- inlet; 2- outlet) 

//return 0; 
//loop through the nodes 
for(int i=0; i<nodes; i++) { 
      sumcond = 0;
      rhs[i]  = 0; 
      // loop through the links
      if(shapf[i] ==0 ) {  //shape factor = 0 
        sumcond = 1.0;
        val.push_back(sumcond);
        row.push_back(i);
        col.push_back(i);
      }
      else {
           for(int j=0; j<shapf[i]; j++) {
     	       id = conn[i][j];
               k1 = l1[id];
     	       k2 = l2[id];
               //printf("%d %d %d\n", k1, k2,i); 
               // link one end connected to inlet/outlet  
	       if(k1 == 0) {
		      if(lp[id]==2) { //outlet
			  tmpj = 0.0; 
                          tmpi = (PI*pow(nr[k2-1],4)/(8.0*mu*n1l[id]))*1e-18;
		      }
		      else if (lp[id]==1){
			  tmpi = 0.0; 
			  tmpj = (PI*pow(nr[k2-1],4)/(8.0*mu*n2l[id]))*1e-18;
		      }
	       }      
/*		      
     	       else if(n1l[id] !=0){
		       tmpi = (PI*pow(nr[k1-1],4)/(8.0*mu*n1l[id]))*1e-18;
	       }

               // link other end connected to inlet/outlet 
               if(k2 == 0) {
		       tmpj = 0;
	       }
     	       else if(n2l[id] !=0){
		       tmpj = (PI*pow(nr[k2-1],4)/(8.0*mu*n2l[id]))*1e-18;
	       }
*/

               // link (ij condutance)
               tmpl = (PI*pow(lr[id],4)/(8.0*mu*ll[id]))*1e-18;
     
              //  cout << i <<" "<<tmpi<<" "<<tmpj<<" "<<tmpl<< endl; 
              //  printf("%d %d %d %d \n", i, id, k1, k2); 
     	      //  printf("%f %f %f \n", nr[k1], nr[k2], lr[id]); 
             //	  printf("%f  %f  %f  %f  %f  %f \n", nr[k1],nr[k2],lr[id],n1l[id], n2l[id], ll[id]); 
             //	  printf("%lf  %f \n", PI, mu);
             //   printf("%d %d %d %f %f %f %f %f %lf\n", k1, k2, id, lr[id],mu,ll[id],tmpi, tmpj, tmpl);
             //   printf("%lf %lf %lf \n", tmpi, tmpj, tmpl);
     	          tmp = tmpi + tmpl + tmpj; 
                  sumcond = sumcond + tmp;
		  //if(tmpi >1000.0 ) std::cout <<" nodei "<< i <<  tmpi <<" "<< n1l[id]<<" "<<nr[k1-1]<<" "<< sumcond << std::endl; 
		  //if(tmpj >1000.0 ) std::cout <<" nodej "<< i <<  tmpj <<" "<< n2l[id]<<" "<<nr[k2-1]<<" "<< sumcond << std::endl; 
		  //if(tmpl >1000.0 ) std::cout <<" nodel "<< i <<  tmpi <<" "<<lr[id]<<" "<<ll[id]<<" "<< sumcond << std::endl; 
     	          if(i+1==k1 && lp[id]==0) {   // k1 -node number  >0  internal node 
                    val.push_back(-tmp);
                    row.push_back(k1-1);
                    col.push_back(k2-1);
                  }
     	          else if (i+1==k2 && lp[id]==0) { // k2 - node number  > 0 
                    val.push_back(-tmp);
                    row.push_back(k2-1);
                    col.push_back(k1-1);
                  }
                  // pressure bc at the inlet 
                  if(lp[id] == 1) { 
     		       //  if any end connected to inlet 
     		      if (k1 == 0 || k2 == 0) rhs[i] = rhs[i] + pinlet*tmp;  
     	          }		    
     
                  // pressure bc at the outlet 
                  if(lp[id] == 2) { 
     		    // if any end connected to outlet 
     		    if (k1 == 0 || k2 == 0) rhs[i] = rhs[i] + poutlet*tmp; 
     	          }		
           }  // connection loop 
      if(sumcond < 1.e-200 && sumcond > -1.e-200) sumcond = 1.0;
      val.push_back(sumcond);
      row.push_back(i);
      col.push_back(i);
      }  // if connection condition  

      //if(i>=12348) { 
       //std::cout << i <<" "<< k1 <<" "<< k2<<" "<<sumcond << std::endl;
      //	      //return 0; 
      //}
//      printf("%d %d %f\n", i, i, sumcond);
}   //for node loop 

printf("print made matrix \n");
//return 0; 

cout << " nodes" << nodes <<" links" << links <<" "<< row.size()<<" "<<col.size()<< val.size()<< endl;
//return 0; 
/*   
   if (link->inlet())
       rhs.back() =(pressureIn * link->avgConductivity());

   else if (link->outlet())
        rhs.back() =(pressureOut * link->avgConductivity());
*/


/*
for(int i=0; i<nodes ; i++) 
    pres[i] = 5.0; 

for(int i=1211; i<nodes ; i++) {
    pres[i] = 10.0; 
}
*/

int maxIters=10000;
double resid; 
int iters;

linear_solver_petsc(val,
                    rhs,
                    pressures,
                    row,
                    col,
                    nodes,
                    maxIters,
                    resid,
                    iters); 
	
printf("maxiter, res %d %lf\n",iters, resid); 

ofstream myfile ("single_pore_out.txt");

for(int i=0; i<nodes-1; i++) {
    myfile << pressures[i] <<" "<<xc[i] <<" "<<yc[i]<<" "<<zc[i]<<" "<<std::endl; 
//    cout   << pressures[i] <<" "<<xc[i] <<" "<<yc[i]<<" "<<zc[i]<<" "<<std::endl; 
}

myfile.close();

double flow=0.0; 

//loop through the nodes 
for(int i=0; i<nodes; i++) { 
      sumcond = 0;
      rhs[i]  = 0; 
      // loop through the links
      if(shapf[i] ==0 ) {  //shape factor = 0 
        sumcond = 1.0;
      }
      else {
           for(int j=0; j<shapf[i]; j++) {
     	       id = conn[i][j];
               k1 = l1[id];
     	       k2 = l2[id];
               //printf("%d %d %d\n", k1, k2,i); 
               // link one end connected to inlet/outlet  
	       if(k1 == 0) {
		      if(lp[id]==2) { //outlet
			  tmpj = 0.0; 
                          tmpi = (PI*pow(nr[k2-1],4)/(8.0*mu*n1l[id]))*1e-18;
		      }
		      else if (lp[id]==1){
			  tmpi = 0.0; 
			  tmpj = (PI*pow(nr[k2-1],4)/(8.0*mu*n2l[id]))*1e-18;
		      }
	       }      

               // link (ij condutance)
               tmpl = (PI*pow(lr[id],4)/(8.0*mu*ll[id]))*1e-18;
     
     	          tmp = tmpi + tmpl + tmpj; 

           }  // connection loop 
      flow = flow + tmp*(pressure[i] - poutlet); 
      }  // if connection condition  

}   //for node loop 

std::cout << "avg flow rate" << flow << std::endl; 

//return 0; 
free(xc); 
free(yc); 
free(zc); 
free(nr); 
free(pres); 

free(l1); 
free(l2); 
free(lr); 
free(ll); 
free(n1l); 
free(n2l); 


ierr = PetscFinalize();
if (ierr) {return ierr;}

} 

