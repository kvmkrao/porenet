
#ifndef WRITE_VTK_H
#define WRITE_VTK_H

#include <vector>

void printNodesVTK(int nodes, std::vector<double> &pressures); 
void printLinksVTK(int nodes, int links, std::vector<double> &pressures); 

#endif 
