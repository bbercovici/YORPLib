
/** MIT License

Copyright (c) 2014 Jay McMahon

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/



//
//  Body.h
//  Spacecraft_SRP
//
//  Created by Jay McMahon on 5/19/14.
//  Copyright (c) 2014 Jay McMahon. All rights reserved.
//

#ifndef __Spacecraft_SRP__Body__
#define __Spacecraft_SRP__Body__

//#include <iostream>
//#include <vector>
#include <string>
//#include "Facet.h"
#include "VoxelGrid.h"

//using namespace std;

class Body
{
public:
    Body();
    Body(string bodyFile);
    Body(string bodyFile, double rho, double spec);

    // Added by Benjamin Bercovici
    Body(std::vector<std::vector<double > > vertices, 
        std::vector<std::vector<int> > facets, 
        double rho, double spec);

    Body(string bodyFile, string optFile);
    int getNumFacets();
    int getNumVerts();
//    void getFacetInfo(int fnum, Facet* fac);
    double* getFacetNormal(int fnum);
    double* getFacetRc(int fnum);
    double getFacetArea(int fnum);
    int* getNeighbors(int fnum);
    vector<int>* getInView(int fnum);
    vector<int>* getInViewRc(int fnum);
    vector<int>* getF2V(int fnum);
    int getVinVox(int voxnum);
    int getFinVox(int voxnum);
    Facet getFacet( int fnum);
    void setNeighbors();
    void setView();
    void setViewRc();
    void setViewGrid();
    void setFacets(double rho, double s);
    void setFacetsVec(vector<double>* rho, vector<double>* s);
    void setF2V();
    void setVoxelGrid(double xmax_in, double ymax_in, double zmax_in, int N);
    int ray_intersection(double* pt, double* uhat, double* hit_pnt, vector<int>* extraFacets);
    double getMaxDimVoxGrid();
//    void getOpticalFourier(int fnum, double& rho, double& s, double& a2, double* nn, double* Ar1, double* Ar2, double* Ar3, double* Ar1_2, double* Ar2_2, double* Ar3_2);
    void getOpticalFourier(int fnum, double& rho, double& s, double& a2, double nn[3][3], double Ar1[3], double Ar2[3][3], double Ar3[3][3], double Ar1_2[3], double Ar2_2[3][3], double Ar3_2[3][3]);
    double getMaxDim(int axis);
private:
    int numFacets, numVerts;
    bool viewCheck; // true if setView has been run, false otherwise
    vector <int> neighbors;
    vector < vector<int>* > inView;
    vector < vector<int>* > inViewRc;
    vector < vector<int>* > inViewGrid;
    vector <int> facetList; // contains the facet numbers as referenced IN THE FILE, so the first vertex = 1...
    vector <double> vertices;
    vector <Facet> facetData;
    vector < vector<int> > f2vlist; // list of facets associate with each vertex (index 1 is associated with vertex 
    VoxelGrid bodyVox;
    double dot ( double* a, double* b);
};

#endif /* defined(__Spacecraft_SRP__Body__) */
