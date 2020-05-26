
#include<cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <sstream>
using namespace std;

int main() { 
double mu1, mu2, r, L, sigma;
double small;
int n;
double x1, t1, t2, p1, p2,A,dp,pc,dt;
double a,b,v,q,ve,R ; 

const double PI  =3.141592653589793238463;

fstream readrad;
readrad.open("inputfile",fstream::in); 
readrad >> r ;
readrad >> mu1 ;
readrad >> sigma ;


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
} 

