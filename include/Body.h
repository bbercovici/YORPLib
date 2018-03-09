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
