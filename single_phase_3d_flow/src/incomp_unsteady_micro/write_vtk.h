
#ifndef WRITE_VTK_H
#define WRITE_VTK_H

#include <vector>
#include <string>
#include <iostream>
void printNodesVTK(std::string fnode, int nodes, double* nd, std::vector<double> &pressures); 
void printLinksVTK(std::string fnode, int nodes, int links, int* l1, int* l2, double* ld, std::vector<double> &pressures); 

#endif // WRITE_VTK_H
