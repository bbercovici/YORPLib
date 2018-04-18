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
//  SRPModel.h
//  Spacecraft_SRP
//
//  Created by Jay McMahon on 5/27/14.
//  Copyright (c) 2014 Jay McMahon. All rights reserved.
//

#ifndef __Spacecraft_SRP__SRPModel__
#define __Spacecraft_SRP__SRPModel__

//#include <iostream>
//#include <vector>
//#include <string>
#include <fstream>
#include "FCoeffs.h"
//#include "Facet.h"
//#include "VoxelGrid.h"
//#include "Body.h"
//#include "Shadow.h"

//using namespace std;

class SRPModel
{
public:
    SRPModel();
    SRPModel(double lambdaDel, double deltaDel, int MaxFourier, Body* bodyExist, int bounces, int numRefine);
    void writeSRPCoeffs(int deltaWrite);
    void writeSRPCoeffsFile(string outputFileBaseName, int numDelta);
private:
    Body Spacecraft;
    vector <double> deltaSunList;
    vector <double> lambdaSunList;
    int numBounces, numGrid, numShadRefine; // if these are 0, then not using this logic (do this to start)
    int FourierOrder;
    vector <Shadow> shadowing; // first index corresponds to lambda_s, each shadow has the value for all facets
    vector < vector <double> > shadowBoundsLow; // 1st facet, 2nd bound case
    vector < vector <double> > shadowBoundsHigh;
    vector < vector<double> > riseLambda; // outer index corresponds to delta_s, inner to each facet
    vector < vector<double> > setLambda;  // outer index corresponds to delta_s, inner to each facet 
    vector < vector <FCoeffs> > Coefficients; // outer index corresponds to delta_s, inner to Fourier order
    void computeMulti();
    vector < vector < double > > fMult;
    vector < vector < double > > tMult;
    vector < vector < vector < double > > > fSpec; // facet, specular case num, force vector
    vector < vector < vector < double > > > tSpec; // facet, specular case num, torque vector
    vector < vector < double > > latSpec; // facet, solar latitude for each specular case num
    vector < vector < double > > longSpec; // facet, solar longitude for each specular case num
    void computeRiseSet(const double* normal, const double latitude, double* rise, double* set);
    double dot ( double* a, double* b);
    void cross(double* a, double* b, double* c);
    void specularMulti(const int dsnum, const int lsnum);
    void compute_fi(const double* uhat, Facet hitF, double* fi);
    void computeShadowBounds(const int dsindex, const int numf, const int lsnum, const double lambdaDel);
//    double I0_1[3]; // = zeros(3,1);
//    double I0_2[3][3]; // = zeros(3);
//    vector <double[3]> Ic_1; // = cell(fourier_max,1);
//    vector <double[3]> Is_1; // = cell(fourier_max,1);
//    vector <double[3][3]> Ic_2; // = cell(fourier_max,1);
//    vector <double[3][3]> Is_2; // = cell(fourier_max,1);
/*    for k = 1:fourier_max
        Ic_1{k} = zeros(3,1);
        Is_1{k} = zeros(3,1);
        Ic_2{k} = zeros(3);
        Is_2{k} = zeros(3);
      end */
};

#endif /* defined(__Spacecraft_SRP__SRPModel__) */
