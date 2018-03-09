//
//  FCoeffs.h
//  Spacecraft_SRP
//
//  Created by Jay McMahon on 5/27/14.
//  Copyright (c) 2014 Jay McMahon. All rights reserved.
//

#ifndef __Spacecraft_SRP__FourierCoeffs__
#define __Spacecraft_SRP__FourierCoeffs__

//#include <iostream>
//#include <vector>
//#include <string>
//#include "Facet.h"
//#include "VoxelGrid.h"
//#include "Body.h"
#include "Shadow.h"

class FCoeffs
{
public:
    FCoeffs();
    FCoeffs(int order);
    void computeCoeffs(const double delta_s, Body* Target, const vector < vector <double> >* low, const vector < vector <double> >* high, const double lamDel, const int numbounces, const vector < vector <double> >* fmult, const vector < vector <double> >* tmult, const vector < vector < double > >* latSpec, const vector < vector < double > >* longSpec, const vector < vector < vector < double > > >* fSpec, const vector < vector < vector < double > > >* tSpec);
    void computeYear();
    void setOrder(int order);
    double* getA();
    double* getB();
    double* getC();
    double* getD();

private:
    int order;
    double A[3];
    double B[3];
    double Ayear[3];
    double Byear[3];
    double C[3];
    double D[3];
    double Cyear[3];
    double Dyear[3];
    double compute_uc1n(double lbnd, double ubnd, int nint);
    double compute_uc2n(double lbnd, double ubnd, int nint);
    double compute_uc3n(double lbnd, double ubnd, int nint);
    double compute_uuc11n(double lbnd, double ubnd, int nint);
    double compute_uuc12n(double lbnd, double ubnd, int nint);
    double compute_uuc22n(double lbnd, double ubnd, int nint);
    double compute_us1n(double lbnd, double ubnd, int nint);
    double compute_us2n(double lbnd, double ubnd, int nint);
    double compute_us3n(double lbnd, double ubnd, int nint);
    double compute_uus11n(double lbnd, double ubnd, int nint);
    double compute_uus12n(double lbnd, double ubnd, int nint);
    double compute_uus22n(double lbnd, double ubnd, int nint);
    double compute_uc1n_shad(double a, double b, int nint, double Ha, double dH);
    double compute_uc2n_shad(double a, double b, int nint, double Ha, double dH);
    double compute_uc3n_shad(double a, double b, int nint, double Ha, double dH);
    double compute_uuc11n_shad(double a, double b, int nint, double Ha, double dH);
    double compute_uuc12n_shad(double a, double b, int nint, double Ha, double dH);
    double compute_uuc22n_shad(double a, double b, int nint, double Ha, double dH);
    double compute_us1n_shad(double a, double b, int nint, double Ha, double dH);
    double compute_us2n_shad(double a, double b, int nint, double Ha, double dH);
    double compute_us3n_shad(double a, double b, int nint, double Ha, double dH);
    double compute_uus11n_shad(double a, double b, int nint, double Ha, double dH);
    double compute_uus12n_shad(double a, double b, int nint, double Ha, double dH);
    double compute_uus22n_shad(double a, double b, int nint, double Ha, double dH);
    double compute_u1_shad(double a, double b, double Ha, double dH);
    double compute_u2_shad(double a, double b, double Ha, double dH);
    double compute_u3_shad(double a, double b, double Ha, double dH);
    double compute_uu11_shad(double a, double b, double Ha, double dH);
    double compute_uu12_shad(double a, double b, double Ha, double dH);
    double compute_uu22_shad(double a, double b, double Ha, double dH);
};

#endif /* defined(__Spacecraft_SRP__FourierCoeffs__) */
