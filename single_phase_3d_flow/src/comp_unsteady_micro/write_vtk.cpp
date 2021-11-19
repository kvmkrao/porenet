
#include<cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <sstream>
#include <vector>
#include <iomanip>      // std::setiosflags
#include <stdlib.h>     /* srand, rand */

/*
instructions to load the vtp files in paraview can be found at
https://openpnm.readthedocs.io/en/latest/getting_started/visualize_in_paraview.html

To visualize the pore data, we need to add some Glyphs. First click on the Glyph button in the tool bar. Then, you can plot the pore data as spheres, where their size and/or their color can be mapped to some variables. In the images below spherical Glyphs are assigned to the pores, where the diameter is linked to the pore diameters and the color to the concentration. Clicking on the Apply button renders these settings.

To visualize the throat data (like diameters or flow rate) we need to set up the following filters in Paraview. First, use the Shrink filter and set it up to 1. Then, the cell data needs to be transposed to the data point with CellDatatoPointdata filter. Then extract the surface with the filter ExtractSurface. Finally, the data may be plotted as tube by using the Tube filter. As previously for the pore data, either the Tube radius or color can be linked to throat data.

Throat data plotted with tubes radius proportional to throat diameter and tubes color linked to the throat flow rate:
*/

using namespace std;
void printNodesVTK(string fnode, int nodes, double* nd, std::vector<double> &pressures) 
{
	size_t id,nn;
	size_t position, shapf; 
	double nr, vol, shapeFactor, clay, xc, yc, zc; 
	double dd1, dd2, dd3,dd4; 
	ifstream nodeloc;
	nodeloc.open(fnode.c_str()); 
        nodeloc >> nn >>  dd1 >>  dd2 >> dd3 >> dd4;
	if(nn =! nodes) {
		cout << "check the input file" << endl; 
		exit(0); 
	}
        cout.precision(3);
        string str1, str2, str3, str4;

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
        nodeloc >> id >> str1 >> str2 >> str3 >> position >> shapf >> str4 >> vol >> shapeFactor >> clay;
        xc = std::stod(str1); 
	yc = std::stod(str2); 
	zc = std::stod(str3); 
        output << xc << "\t" << yc << "\t" << zc << "\t";
    }
    nodeloc.close();
    output << endl;
    output << "             </DataArray>" << endl;
    output << "         </Points>" << endl;
    //Node Radii
    output << "         <PointData Scalars=\"Point Data\">" << endl;
    output << "             <DataArray type=\"Float32\" Name=\"diameter\" format=\"ascii\">" << endl;
    for(int k=0; k<nodes; ++k)
    {
        output << nd[k] << "\t";
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


void printLinksVTK(string fnode, int nodes, int links, int* l1, int* l2, double* ld, std::vector<double> &pressures) 
{
	size_t id,nn;
	size_t position, shapf; 
	double nr, vol, shapeFactor, clay, xc, yc, zc; 
	double dd1, dd2, dd3,dd4; 
	ifstream nodeloc;
	nodeloc.open(fnode.c_str()); 
        nodeloc >> nn >>  dd1 >>  dd2 >> dd3 >> dd4;
	if(nn =! nodes) {
		cout << "check the input file" << endl; 
		exit(0); 
	}
        cout.precision(3);
        string str1, str2, str3, str4;

    int icount = 0;
    for(int k=0; k<links; k++)
            if (l1[k] > 0 && l2[k]>0) icount = icount + 1;

    ofstream output ("single_pore_links.vtp");
    output << "<VTKFile byte_order=\"LittleEndian\" type=\"PolyData\" version=\"0.1\">" << endl;
    output << " <PolyData>" << endl;
    output << "     <Piece  NumberOfLines=\"" << icount << "\" NumberOfPoints=\"" << nodes << "\">" << endl;
    output << "         <Points>" << endl;
    output << "             <DataArray Name=\"nodes\" NumberOfComponents=\"3\" type=\"Float64\">" << endl;
  //  output << setprecision(4) << setiosflags(ios::scientific);

    for(int k=0; k<nodes; ++k)
    {
	nodeloc >> id >> str1 >> str2 >> str3 >> position >> shapf >> str4 >> vol >> shapeFactor >> clay;
        xc = std::stod(str1); 
	yc = std::stod(str2); 
	zc = std::stod(str3); 
        output << xc << "\t" << yc << "\t" << zc << "\t";
    }
    nodeloc.close(); 

    output << endl;
    output << "             </DataArray>" << endl;
    output << "         </Points>" << endl;

    output << "         <Lines>" << endl;
    output << "             <DataArray Name=\"connectivity\" NumberOfComponents=\"1\" type=\"Int64\">" << endl;

    for(int k=0; k<links; k++)
	    if(l1[k] > 0 && l2[k] > 0) output << l1[k]-1 << "\t" << l2[k]-1 << "\t";
    output << endl;
    output << "             </DataArray>" << endl;
    output << "             <DataArray Name=\"offsets\" NumberOfComponents=\"1\" type=\"Int64\">" << endl;

    cout << " read links.dat \n" << endl;

    for(int k=0; k<icount; k++)
    {
        output << 2*(k+1) << "\t";
    }
    output << endl;
    output << "             </DataArray>" << endl;
    output << "         </Lines>" << endl;

    output << "         <CellData Scalars=\"Cell Data\">" << endl;
    output << "             <DataArray type=\"Float32\" Name=\"radii\" format=\"ascii\">" << endl;

    cout << " write information \n" << endl;

    for(int k=0; k<links; ++k)
    {
        if(l1[k] > 0 && l2[k] > 0) output << ld[k] << "\t";
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
    cout << " complete \n" << endl; 
    output << endl;   
    output.close();
    return; 
}

