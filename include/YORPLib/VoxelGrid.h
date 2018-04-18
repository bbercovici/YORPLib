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
//  VoxelGrid.h
//  Spacecraft_SRP
//
//  Created by Jay McMahon on 5/19/14.
//  Copyright (c) 2014 Jay McMahon. All rights reserved.
//

#ifndef __Spacecraft_SRP__VoxelGrid__
#define __Spacecraft_SRP__VoxelGrid__

//#include <iostream>
//#include <vector>
#include "Facet.h"

//using namespace std;
namespace YORPLib{

    class VoxelGrid
    {
    public:
        VoxelGrid();
        VoxelGrid(double xmax_in, double ymax_in, double zmax_in, int N);
        void fillGrid(vector<int>* facets, vector<double>* vertices, int numVerts, vector< vector<int> >* f2vlist, vector<Facet>* fdata);
        int ray_intersect_voxel(double* pt, double* uhat, vector<int>* facets, vector<double>* vertices, vector <Facet>* facetData, double* hit_pnt, vector<int>* extraFacets);
        void find_pt_in_voxel(double* pt, vector<int>* voxel_list);
        int test_intersection(double* pt, double* uhat, vector<int>* facets, vector<double>* vertices, vector <Facet>* facetData, double*hit_pnt, int current_voxel, vector<int>* extraFacets);
        int getNumFinVox(int VoxNum);
        int getNumVinVox(int VoxNum);
        double getMaxDim();
        void getVoxInd(int voxnum, int* x, int* y, int* z);
        bool intersectTriVox(int voxNum, double* v1, double* v2, double* v3, double* nhat);
    private:
        vector <double> xmax;
        vector <double> xmin;
        vector <double> ymax;
        vector <double> ymin;
        vector <double> zmax;
        vector <double> zmin;
        vector <int> posx;
        vector <int> negx;
        vector <int> posy;
        vector <int> negy;
        vector <int> posz;
        vector <int> negz;
        vector <int> fnum;
        vector <int> vnum;
        vector < vector <int> > flist;
        vector < vector <int> > vlist;
        double dx, dy, dz;
        int N;
        void inverse ( double* answer, double* x, double* y, double* n);
        bool in_triangle_3D( double *A, double *B, double *C, double *P, double *nhat);
        double dot ( double* a, double* b);
        bool satAABBTri(double rBox, double* axis, double* v1, double* v2, double* v3);
        void cross(double* a, double* b, double* c);
    };
}
#endif /* defined(__Spacecraft_SRP__VoxelGrid__) */
