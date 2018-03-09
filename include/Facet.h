//
//  Facet.h
//  Spacecraft_SRP
//
//  Created by Jay McMahon on 5/19/14.
//  Copyright (c) 2014 Jay McMahon. All rights reserved.
//

#ifndef __Spacecraft_SRP__Facet__
#define __Spacecraft_SRP__Facet__

#include <iostream>
#include <vector>

using namespace std;

class Facet
{
public:
    Facet();
    Facet(int* verts, vector<double>* vertices, double rho, double s);
    void updateRhoS(double rhoSet, double sSet);
    void computeOptical(); 
    int getVert(int vnum);
    double getRho();
    double getS();
    double getArea();
    double* getRc();
    double* getNormal();
    double geta2();
    void getOpticalFourier(double nnOut[3][3], double Ar1Out[3], double Ar2Out[3][3], double Ar3Out[3][3], double Ar1_2Out[3], double Ar2_2Out[3][3], double Ar3_2Out[3][3]);
//    void getOpticalFourier(double* nnOut, double* Ar1Out, double* Ar2Out, double* Ar3Out, double* Ar1_2Out, double* Ar2_2Out, double* Ar3_2Out);
private:
    int vert_nums[3]; // index number of first element of each vertex (vertex 1 spans indices 0-2 in the vector)
    //double v1[3], v2[3], v3[3];
    double r_c[3];
    double normal[3];
    double rho;
    double s;
    double area;
// Below are the properties needed for the Lambertian BRDF; consider packaging these into some sort of BRDF class so that it will be easier in the future to switch out BRDFs
//    double B; // Assumed = 2/3 for now; just doing Lambertian to start
    double a2;
    double nn[3][3]; // = normals(:,j)*normals(:,j)';
    double Ar1[3]; // = -areas(j)./pi.*normals(:,j);
    double Ar2[3][3]; // = -areas(j)./pi.*nn;
    double Ar3[3][3]; // = -areas(j)./pi.*(2.*nn - eye(3));
    double Ar1_2[3]; // = Ar1./2;
    double Ar2_2[3][3]; // = Ar2./2;
    double Ar3_2[3][3]; // = Ar3./2;
    // Private member functions
    void computeGeometry(double* v1, double* v2, double* v3, double* normal, double* r_c, double* area);
};

#endif /* defined(__Spacecraft_SRP__Facet__) */
