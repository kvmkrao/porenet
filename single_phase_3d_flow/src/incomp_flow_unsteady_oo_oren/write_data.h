
#ifndef WRITE_DATA_H
#define WRITE_DATA_H

#include <vector>

class writeData 
{
    public: 
        writeData(); 
	~writeData(); 
	void NodesVTK(int nodes, std::vector<double> &pressures); 
	void LinksVTK(int nodes, int links, std::vector<double> &pressures); 
}; 

#endif
