
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

double mu1, mu2, r, L, sigma;
double small;
int n,nodes,links;
double x1, t1, t2, p1, p2,A,dp,pc,dt;
double a,b,v,q,ve,R ; 

const double PI  =3.141592653589793238463;

fstream readrad;
readrad.open("inputfile",fstream::in); 
readrad >> r ;
readrad >> mu1 ;
readrad >> sigma ;

fstream rnode;
rnode.open("node.dat",fstream::in);
rnode >> nodes; 

fstream rnodedia;
rnodedia.open("pore_diameter",fstream::in);
double* xc = new double[nodes]; 
double* yc = new double[nodes]; 
double* zc = new double[nodes]; 
double* nd = new double[nodes]; 
double* pres = new double[nodes]; 

for(int i=0; i<nodes ; i++) {
	rnode>> xc[i] >> yc[i] >> zc[i]; 
        rnodedia >> nd[i];
//        cout << i<<"\t"<< xc[i] <<"\t" << yc[i] <<"\t"<< zc[i] <<"\t"<< nd[i] << endl; 
}

fstream rlink;
rlink.open("link.dat",fstream::in);
rlink >> links;

fstream rlinkdia;
rlinkdia.open("throat_diameter",fstream::in);
int* l1 = new int[links];
int* l2 = new int[links];
double* ld = new double[links];

std::vector<double> val; 
val.reserve(nodes*7);

std::vector<double> rhs;
rhs.reserve(nodes);

std::vector<double> pressures;
pressures.reserve(nodes); 

std::vector<int> row;
row.reserve(nodes*7); 

std::vector<int> col;
col.reserve(nodes*7); 

double sumcond,tmp; 
double mu  = 1.5e-5; 
double len = 1e-4; 
int k1, k2; 

for(int i=0; i<links ; i++) {
        rlink >> l1[i] >> l2[i];
        rlinkdia >> ld[i];
        //cout << i<<"\t"<< l1[i] <<"\t" << l2[i] <<"\t"<< ld[i]  << endl;
}

for(int i=0; i<nodes; i++)
   rhs.push_back(0);

//loop through the nodes 
for(int i=0; i<nodes; i++) { 
      sumcond = 0; 
      // loop through the links 
      for(int j=0; j<links; j++) { 
          k1 = l1[j]; 
	  k2 = l2[j]; 
          //printf("%d %d %d\n", k1, k2,i); 
	  if(i==k1 || i==k2) {
            tmp = (PI*pow(100.*ld[j]/2.0,4)/(8.0*mu*len));
            printf("%d %d %f %f %f %lf\n", k1, k2, ld[j],mu,len,tmp); 
	    if(i==k1) {
              val.push_back(-tmp);
              row.push_back(k1);
              col.push_back(k2);
            }
	    else if (i==k2) {
              val.push_back(-tmp);
              row.push_back(k2);
              col.push_back(k1);
            }
            sumcond += tmp;
          }
// pressure bc at the outlet 
      if(xc[k1]>1.0e-3) rhs[k1] = 10.0* tmp;  
      if(xc[k2]>1.0e-3) rhs[k2] = 10.0* tmp;  
      }
      if(sumcond < 1.e-200 && sumcond > -1.e-200) sumcond = 1.0;
      val.push_back(sumcond);
      row.push_back(i);
      col.push_back(i);
      printf("%d %d %f\n", i, i, sumcond);
}

printf("print made matrix \n");

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

int maxIters;
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

return 0; 
//mu1 = 6.69e-03; 
mu2 = 0*1e-6;  L=10e-6; 
//sigma=1.60;
//r=1E-6;
small=1e-2;
n=10; x1=small; t1=small; t2=10;  p1=0; p2=0;

A=PI*pow(r,2);
pc= 2*PI*r*sigma/A;   dp=p1-p2+pc;

dt = 0.1 ; 
n  = (50 - t1)/dt ; 

double t[n+1];
double xe[n+1];
double  x[n+1];

x[0]=x1; t[0]=t1; xe[0]=x1;

//nv=50;
//for(int  n; n <51 ; n=n+10) { 
//dt=(t2-t1)/n;

ofstream myfile;
myfile.open ("single_pore_out.txt");

for(int i=0; i <=n; i=i+1) {
    a=8*(mu1*x[i] + mu2*(L-x[i])) ; b=PI*pow(r,4);
    q = dp * b/a;  v=q/A; // q=pi*r^4 * dp / (8  (x*mu1 + (L-x)*mu2) )
    x[i+1] = x[i] + v*dt;
    t[i+1] = t[i] + dt;
    xe[i+1] = sqrt(x1*x1 + r*sigma*(t[i+1]-t1) / (2*mu1));
    ve = (xe[i+1]-xe[i])/dt;
    std::cout << t[i+1] << " "<< x[i+1] <<" "<< q <<" "<< v <<" "<< ve << std::endl; 
    myfile <<  t[i+1]   << " "<< x[i+1] <<" "<< q <<" "<< v <<" "<< ve << std::endl; 
}

myfile.close();

ierr = PetscFinalize();
if (ierr) {return ierr;}

} 

