
#ifndef WRITE_VTK_H
#define WRITE_VTK_H

#include <vector>

void printNodesVTK(int nodes,double* xc, double* yc, double* zc, double* nd, std::vector<double> &pressures); 
void printLinksVTK(int nodes, int links, double* xc, double* yc, double* zc, int* l1, int* l2, double* ld, std::vector<double> &pressures); 

#endif // WRITE_VTK_H
