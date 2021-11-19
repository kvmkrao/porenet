
#include<cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <sstream>
#include <vector>
#include <iomanip>      // std::setiosflags
#include <stdlib.h>     /* srand, rand */

#include "sphase_unsteady_comp.hpp"

using namespace std;
void printNodesVTK(int nodes, std::vector<double> &pressures) 
{
    Node Body;  
    Body.readdata();  
    ofstream output ("single_pore_nodes.vtp");
    output << "<VTKFile byte_order=\"LittleEndian\" type=\"PolyData\" version=\"0.1\">" << endl;
    output << " <PolyData>" << endl;
    output << "     <Piece   NumberOfPoints=\"" << nodes << "\">" << endl;
    output << "         <Points>" << endl;
    output << "             <DataArray Name=\"coords\" NumberOfComponents=\"3\" type=\"Float64\">" << endl;

    //Node Positions
    //output << setprecision(4) << setiosflags(ios::scientific);
    for(int k=0; k<nodes; ++k)
    {
        output << Body.xc[k] << "\t" << Body.yc[k] << "\t" << Body.zc[k] << "\t";
    }
    output << endl;
    output << "             </DataArray>" << endl;
    output << "         </Points>" << endl;
    //Node Radii
    output << "         <PointData Scalars=\"Point Data\">" << endl;
    output << "             <DataArray type=\"Float32\" Name=\"diameter\" format=\"ascii\">" << endl;
    for(int k=0; k<nodes; ++k)
    {
        output << Body.rad[k] << "\t";
    }
    output << endl;
    output << "             </DataArray>" << endl;

    output << "             <DataArray type=\"Float32\" Name=\"pressure\" format=\"ascii\">" << endl;
    for(int k=0; k<nodes; ++k)
    {
        output << pressures[k] << "\t";
    }
    output << endl;
    output << "             </DataArray>" << endl;
    output << "         </PointData>" << endl;
    output << "     </Piece>" << endl;
    output << " </PolyData>" << endl;
    output << "</VTKFile>" << endl;
    output << endl;
    output.close();
    return; 
}


void printLinksVTK(int nodes, int links, std::vector<double> &pressures) 
{

    Node Body; 
    Body.readdata(); 
    Link Throat; 
    Throat.readdata(); 
    ofstream output ("single_pore_links.vtp");
    output << "<VTKFile byte_order=\"LittleEndian\" type=\"PolyData\" version=\"0.1\">" << endl;
    output << " <PolyData>" << endl;
    output << "     <Piece  NumberOfLines=\"" << links << "\" NumberOfPoints=\"" << nodes << "\">" << endl;
    output << "         <Points>" << endl;
    output << "             <DataArray Name=\"nodes\" NumberOfComponents=\"3\" type=\"Float64\">" << endl;
  //  output << setprecision(4) << setiosflags(ios::scientific);

    for(int k=0; k<nodes; ++k)
    {
        output << Body.xc[k] << "\t" << Body.yc[k] << "\t" << Body.zc[k] << "\t";
    }
    output << endl;
    output << "             </DataArray>" << endl;
    output << "         </Points>" << endl;

    output << "         <Lines>" << endl;
    output << "             <DataArray Name=\"connectivity\" NumberOfComponents=\"1\" type=\"Int64\">" << endl;

    for(int k=0; k<links; k++)
    {
        output << Throat.ends.at(k).first  << "\t" << Throat.ends.at(k).second << "\t";
    }
    output << endl;
    output << "             </DataArray>" << endl;
    output << "             <DataArray Name=\"offsets\" NumberOfComponents=\"1\" type=\"Int64\">" << endl;

    cout << " read links.dat \n" << endl;

    for(int k=0; k<links; k++)
    {
        output << 2*(k+1) << "\t";
    }
    output << endl;
    output << "             </DataArray>" << endl;
    output << "         </Lines>" << endl;

    output << "         <CellData Scalars=\"Cell Data\">" << endl;
    output << "             <DataArray type=\"Float32\" Name=\"radii\" format=\"ascii\">" << endl;

    cout << " write information \n" << endl;

    //double tmp = 100.0; 
    for(int k=0; k<links; ++k)
    {
        //output << tmp << "\t";
        output << Throat.rad[k] << "\t";
    }
    output << endl;
    output << "             </DataArray>" << endl;
/*
    output << "             <DataArray type=\"Float32\" Name=\"test\" format=\"ascii\">" << endl;

    cout << " write information 2 \n" << endl;

    for(int k=0; links; ++k)
    {
        output << rand() << "\t";
    }

    cout << " write information 2 \n" << endl;
    output << endl;
    output << "             </DataArray>" << endl;

*/
    output << "         </CellData>" << endl;
    output << "     </Piece>" << endl;
    output << " </PolyData>" << endl;
    output << "</VTKFile>" << endl;
//    exit(0);  
    cout << " complete \n" << endl; 
    output << endl;   
    output.close();
    return; 
}

