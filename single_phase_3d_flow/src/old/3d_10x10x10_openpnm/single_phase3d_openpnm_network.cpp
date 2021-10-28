
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


#include "linear_solver_petsc.h"

using namespace std;


//int main() { 
int main(int argc, char *argv[])
{
static char help[] = "Solves a linear system in parallel with KSP.\n";

PetscErrorCode ierr;
ierr = PetscInitialize(&argc,&argv,(char*)0,help);


double mu1, mu2, r, L, sigma;
double small;
int n,nodes,links;
double x1, t1, t2, p1, p2,A,dp,pc,dt;
double pinlet = 1000.0;    // 1atm 
double poutlet= 50.0;  //101325;  // 5atm 

const double PI  =3.141592653589793238463;

ifstream rnode("node.dat");
rnode >> nodes; 

ifstream rnodedia("pore_diameter");
double* xc   = (double*)malloc(nodes * sizeof(double*));
double* yc   = (double*)malloc(nodes * sizeof(double*));
double* zc   = (double*)malloc(nodes * sizeof(double*));
double* nd   = (double*)malloc(nodes * sizeof(double*));
double* pres = (double*)malloc(nodes * sizeof(double*));

for(int i=0; i<nodes ; i++) {
	rnode>> xc[i] >> yc[i] >> zc[i]; 
        rnodedia >> nd[i];
//        cout << i<<"\t"<< xc[i] <<"\t" << yc[i] <<"\t"<< zc[i] <<"\t"<< nd[i] << endl; 
}
ifstream rlink("link.dat");
rlink >> links;

ifstream rlinkdia("throat_diameter");
ifstream condl1("throat.conduit_lengths.pore1"); 
ifstream condl2("throat.conduit_lengths.pore2"); 
ifstream condll("throat.conduit_lengths.throat"); 
ifstream condla("throat_area"); 

int* l1     = (int*)malloc(links * sizeof(int*));
int* l2     = (int*)malloc(links * sizeof(int*));
int* lp     = (int*)malloc(links * sizeof(int*));
double*  ld = (double*)malloc(links * sizeof(double*));
double*  ll = (double*)malloc(links * sizeof(double*));
double* n1l = (double*)malloc(links * sizeof(double*));
double* n2l = (double*)malloc(links * sizeof(double*));
double*  la = (double*)malloc(links * sizeof(double*));

std::vector<double> val; 
val.reserve(links);

std::vector<double> rhs;
rhs.reserve(nodes);

std::vector<double> pressures;
pressures.reserve(nodes); 

std::vector<int> row;
row.reserve(links); 

std::vector<int> col;
col.reserve(links); 

std::vector<int> shapf;
shapf.reserve(nodes);

double sumcond,tmp,tmpi,tmpj,tmpl,position,vol,shapeFactor,clay,distance;
double mu  = 1.002e-3;
double len = 1e-4;
int k1, k2,ic,id;

double porevol = 0; 
for(int i=0; i<links ; i++) {
        rlink >> l1[i] >> l2[i];
        rlinkdia >> ld[i];
        condl1   >> n1l[i]; 
        condl2   >> n2l[i]; 
        condll   >> ll[i]; 
        condla   >> la[i]; 
        //cout     << i<<"\t"<< l1[i] <<"\t" << l2[i] <<"\t"<< ld[i]<<"\t"<<n1l[i]<<"\t"<<n2l[i]<<"\t"<<ll[i]<<"\t"<<la[i] << endl;
}
//return 0; 

for(int i=0; i<nodes; i++) {
   rhs.push_back(0);
   pres[i] = pinlet;
   pressures.push_back(pinlet);
}

int conn[nodes][10];

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
      if(i==k1 || i==k2) {
          conn[i][ic] = j;
          ic = ic + 1;
      }
  }
  shapf[i] = ic;
  //cout <<" node shape factor " << i  <<"  "<<shapf[i] << endl;   
}

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
               //tmpi = (PI*pow(nd[k1]/2.0,4)/(8.0*mu*n1l[id]));
               tmpi = (PI*pow(nd[k1]/2.0,3)/(8.0*mu));
               tmpj = (PI*pow(nd[k2]/2.0,3)/(8.0*mu));
               //tmpj = (PI*pow(nd[k2]/2.0,4)/(8.0*mu*n2l[id]));

               // link (ij condutance)
               //tmpl = (PI*pow(ld[id]/2.0,4)/(8.0*mu*ld[id]));
               tmpl = (la[id]*pow(ld[id]/2.0,2)/(8.0*mu*ll[id]));

                  tmp = tmpi + tmpl + tmpj;
                  sumcond = sumcond + tmp;
                  if(i==k1) {   // k1 -node number  >0  internal node 
                      val.push_back(-tmp);
                      row.push_back(k1);
                      col.push_back(k2);
                  }
                  else if (i==k2) { // k2 - node number  > 0 
                      val.push_back(-tmp);
                      row.push_back(k2);
                      col.push_back(k1);
                  }
                  // pressure bc at the inlet 
                  if(xc[i] <= 5.0e-05) {
                      //  if any end connected to inlet 
                      rhs[i] = rhs[i] + pinlet*tmp;
                      //cout  << xc[i] <<"\t"<< pinlet << endl; 
                  }

                  // pressure bc at the outlet 
                  if(xc[i] >=1.0e-03) {
                      // if any end connected to outlet 
                      rhs[i] = rhs[i] + poutlet*tmp;
                      //cout  << xc[i] <<"\t"<< poutlet << endl; 
                  }
           }  // connection loop 
      if(sumcond < 1.e-200 && sumcond > -1.e-200) sumcond = 1.0;
      val.push_back(sumcond);
      row.push_back(i);
      col.push_back(i);
      }  // if connection condition  

}   //for node loop 

//return 0; 
/*
for(int i=0; i<nodes ; i++) 
    pres[i] = 5.0; 

for(int i=1211; i<nodes ; i++) {
    pres[i] = 10.0; 
}
*/

int maxIters=50000;
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
	

printf("maxiter, res %d %f\n",iters, resid); 

ofstream myfile ("single_pore_out.txt");

for(int i=0; i<nodes; i++) {
    myfile << pressures[i] <<" "<<xc[i] <<" "<<yc[i]<<" "<<zc[i]<<" "<<std::endl;
   // cout   << pressures[i] <<" "<<xc[i] <<" "<<yc[i]<<" "<<zc[i]<<" "<<std::endl; 
}

myfile.close();
free(xc);
free(yc);
free(zc);
free(nd);
free(pres);

free(l1);
free(l2);
free(ld);
free(ll);
free(la);
free(n1l);
free(n2l);

ierr = PetscFinalize();
if (ierr) {return ierr;}
} 

