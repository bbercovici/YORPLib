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
//  FCoeffs.cpp
//  Spacecraft_SRP
//
//  Created by Jay McMahon on 5/27/14.
//  Copyright (c) 2014 Jay McMahon. All rights reserved.
//

#include "FCoeffs.h"
#include "Body.h"
#include <math.h>
#include <algorithm>

FCoeffs::FCoeffs()
{
    FCoeffs(1);
}

FCoeffs::FCoeffs(int orderIn)
{
    order = orderIn;
    for (int ii=0; ii<3; ii++) {
        A[ii] = 0.0;
        B[ii] = 0.0;
        Ayear[ii] = 0.0;
        Byear[ii] = 0.0;
        C[ii] = 0.0;
        D[ii] = 0.0;
        Cyear[ii] = 0.0;
        Dyear[ii] = 0.0;
    }
}

void FCoeffs::computeCoeffs(const double delta_s, Body* Target, const vector < vector <double> >* low, const vector < vector <double> >* high, const double lamDel, const int numbounces, const vector < vector <double> >* fmult, const vector < vector <double> >* tmult, const vector < vector < double > >* latSpec, const vector < vector < double > >* longSpec, const vector < vector < vector < double > > >* fSpec, const vector < vector < vector < double > > >* tSpec)
{
    int numf = Target->getNumFacets();
    double cdels = cos(delta_s);
    double sdels = sin(delta_s);
    double cdels2 = cdels*cdels;
    double sdels2 = sdels*sdels;
    double csdels = cdels*sdels;
    double u1, u2, u3, uu11, uu12, uu22;
    double I0_1[3], I0_2[3][3];
    int jj, kk, ll, bound_count;
    double rho, s, a2;
    double nn[3][3], Ar1[3], Ar2[3][3], Ar3[3][3];
    double Ar1_2[3], Ar2_2[3][3], Ar3_2[3][3];
    double nhat[3], rc[3];
    double* nhat_p;
    double* rc_p;
//    double* rho_p;
    double A0[3], A1[3], A2[3], temp[3];
    double B0[3], B1[3], B2[3], temp2[3];
    double uc1n, uc2n, uc3n, us1n, us2n, us3n;
    double uuc11n, uuc12n, uuc22n, uus11n, uus12n, uus22n;
    double Ic_1[3], Is_1[3], Ic_2[3][3], Is_2[3][3];
    double A01m, A1m, B1m;
    vector <double>::const_iterator iter_d;
    int specIndex;
    
    for (jj=0; jj<3; jj++) {
        A[jj] = 0.0;
        B[jj] = 0.0;
        C[jj] = 0.0;
        D[jj] = 0.0;
    }
    
    for (int ii=0; ii<numf; ii++) {
        // test to see if size of third index for this facet in bounds is zero - if so skip this facet b/c not lit
        if (!(*low)[ii].empty()) {
            // get optical data for this facet
            Target->getOpticalFourier(ii, rho, s, a2, nn, Ar1, Ar2, Ar3, Ar1_2, Ar2_2, Ar3_2);
            nhat_p = Target->getFacetNormal(ii);
            nhat[0] = *nhat_p;
            nhat[1] = *(nhat_p + 1);
            nhat[2] = *(nhat_p + 2);
            rc_p = Target->getFacetRc(ii);
            rc[0] = *rc_p;
            rc[1] = *(rc_p+1);
            rc[2] = *(rc_p+2);
            // get number of segments
            bound_count = (*low)[ii].size();
            double dHall = 1.0/lamDel;
            if (order==0) {
                // Compute integrals analytically over each unshadowed segment
                u1 = 0.0;
                u2 = 0.0;
                u3 = 0.0;
                uu11 = 0.0;
                uu12 = 0.0;
                uu22 = 0.0;
                for (kk = 0; kk<bound_count; kk++) {
                    u1 += sin((*high)[ii][kk]) - sin((*low)[ii][kk]);
                    u2 += -1.0*cos((*high)[ii][kk]) + cos((*low)[ii][kk]);
                    u3 += (*high)[ii][kk] - (*low)[ii][kk];
                    uu11 += ((*high)[ii][kk]/2) + (sin(2*(*high)[ii][kk])/4) - ((*low)[ii][kk]/2) - (sin(2*(*low)[ii][kk])/4);
                    uu12 += (-0.5*cos((*high)[ii][kk])*cos((*high)[ii][kk])) + (.5*cos((*low)[ii][kk])*cos((*low)[ii][kk]));
                    uu22 += ((*high)[ii][kk]/2) - (sin(2*(*high)[ii][kk])/4) - ((*low)[ii][kk]/2) + (sin(2*(*low)[ii][kk])/4);
                    // Compute the partial shadowing before the low boundary
                    if ((*low)[ii][kk] > 0.0) {
                        u1 += compute_u1_shad((*low)[ii][kk] - lamDel, (*low)[ii][kk], 0.0, dHall);
                        u2 += compute_u2_shad((*low)[ii][kk] - lamDel, (*low)[ii][kk], 0.0, dHall);
                        u3 += compute_u3_shad((*low)[ii][kk] - lamDel, (*low)[ii][kk], 0.0, dHall);
                        uu11 += compute_uu11_shad((*low)[ii][kk] - lamDel, (*low)[ii][kk], 0.0, dHall);
                        uu12 += compute_uu12_shad((*low)[ii][kk] - lamDel, (*low)[ii][kk], 0.0, dHall);
                        uu22 += compute_uu22_shad((*low)[ii][kk] - lamDel, (*low)[ii][kk], 0.0, dHall);
                    }
                    // Compute the partial shadowing after the high boundary
                    if ((*high)[ii][kk] < 2.0*M_PI) {
                        u1 += compute_u1_shad((*high)[ii][kk], (*high)[ii][kk]+lamDel, 1.0, -dHall);
                        u2 += compute_u2_shad((*high)[ii][kk], (*high)[ii][kk]+lamDel, 1.0, -dHall);
                        u3 += compute_u3_shad((*high)[ii][kk], (*high)[ii][kk]+lamDel, 1.0, -dHall);
                        uu11 += compute_uu11_shad((*high)[ii][kk], (*high)[ii][kk]+lamDel, 1.0, -dHall);
                        uu12 += compute_uu12_shad((*high)[ii][kk], (*high)[ii][kk]+lamDel, 1.0, -dHall);
                        uu22 += compute_uu22_shad((*high)[ii][kk], (*high)[ii][kk]+lamDel, 1.0, -dHall);
                    }
                }
                
                // Compute zero order integrals
                I0_1[0] = cdels*u1;
                I0_1[1] = cdels*u2;
                I0_1[2] = sdels*u3;
                
                I0_2[0][0] = cdels2*uu11;
                I0_2[1][0] = cdels2*uu12;
                I0_2[2][0] = csdels*u1;
                I0_2[0][1] = I0_2[1][0];
                I0_2[1][1] = cdels2*uu22;
                I0_2[2][1] = csdels*u2;
                I0_2[0][2] = I0_2[2][0];
                I0_2[1][2] = I0_2[2][1];
                I0_2[2][2] = sdels2*u3;
                
                // Compute Total Coefficients
                A01m = 0.0;
                for (jj=0; jj<3; jj++) {
                    A0[jj] = 0.0;
                    A1[jj] = 0.0;
                    A2[jj] = 0.0;
                    for (kk=0; kk<3; kk++) {
                        A0[jj] += I0_2[jj][kk]*Ar1_2[kk];
                        A1[jj] += Ar2_2[jj][kk]*I0_1[kk];
                        for (ll=0; ll<3; ll++) {
                            A2[jj] += Ar3_2[jj][kk]*I0_2[kk][ll]*nhat[ll];
                        }
                        if ((numbounces>0) && (jj==0)) {
                            A01m += Ar1_2[kk]*I0_1[kk];
                        }
                    }
                    temp[jj] = A0[jj] + a2*A1[jj] + rho*s*A2[jj];
                    A[jj] += temp[jj];
                    if (numbounces>0) {
                        A[jj] += A01m*(*fmult)[ii][jj];
                    }
                }
                C[0] += rc[1]*temp[2] - rc[2]*temp[1];
                C[1] += rc[2]*temp[0] - rc[0]*temp[2];
                C[2] += rc[0]*temp[1] - rc[1]*temp[0];
                if (numbounces>0) {
                    for (jj=0; jj<3; jj++) {
                        C[jj] += A01m*(*tmult)[ii][jj];
                    }
                }
            } else {
                // Compute integrals analytically over each unshadowed segment
                if (order == 1) {
                    u1 = 0.0;
                    u2 = 0.0;
                    uu11 = 0.0;
                    uu12 = 0.0;
                    uu22 = 0.0;
                    for (kk = 0; kk<bound_count; kk++) {
                        u1 += sin((*high)[ii][kk]) - sin((*low)[ii][kk]);
                        u2 += -1.0*cos((*high)[ii][kk]) + cos((*low)[ii][kk]);
                        uu11 += ((*high)[ii][kk]/2) + (sin(2*(*high)[ii][kk])/4) - ((*low)[ii][kk]/2) - (sin(2*(*low)[ii][kk])/4);
                        uu12 += (-0.5*cos((*high)[ii][kk])*cos((*high)[ii][kk])) + (.5*cos((*low)[ii][kk])*cos((*low)[ii][kk]));
                        uu22 += ((*high)[ii][kk]/2) - (sin(2*(*high)[ii][kk])/4) - ((*low)[ii][kk]/2) + (sin(2*(*low)[ii][kk])/4);
                        // Compute the partial shadowing before the low boundary
                        if ((*low)[ii][kk] > 0.0) {
                            u1 += compute_u1_shad((*low)[ii][kk] - lamDel, (*low)[ii][kk], 0.0, dHall);
                            u2 += compute_u2_shad((*low)[ii][kk] - lamDel, (*low)[ii][kk], 0.0, dHall);
                            uu11 += compute_uu11_shad((*low)[ii][kk] - lamDel, (*low)[ii][kk], 0.0, dHall);
                            uu12 += compute_uu12_shad((*low)[ii][kk] - lamDel, (*low)[ii][kk], 0.0, dHall);
                            uu22 += compute_uu22_shad((*low)[ii][kk] - lamDel, (*low)[ii][kk], 0.0, dHall);
                        }
                        // Compute the partial shadowing after the high boundary
                        if ((*high)[ii][kk] < 2.0*M_PI) {
                            u1 += compute_u1_shad((*high)[ii][kk], (*high)[ii][kk]+lamDel, 1.0, -dHall);
                            u2 += compute_u2_shad((*high)[ii][kk], (*high)[ii][kk]+lamDel, 1.0, -dHall);
                            uu11 += compute_uu11_shad((*high)[ii][kk], (*high)[ii][kk]+lamDel, 1.0, -dHall);
                            uu12 += compute_uu12_shad((*high)[ii][kk], (*high)[ii][kk]+lamDel, 1.0, -dHall);
                            uu22 += compute_uu22_shad((*high)[ii][kk], (*high)[ii][kk]+lamDel, 1.0, -dHall);
                        }
                    }
                    uc1n = uu11;
                    uc2n = uu12;
                    uc3n = u1;
                    us1n = uu12;
                    us2n = uu22;
                    us3n = u2;
                } else {
                    uc1n = 0.0;
                    uc2n = 0.0;
                    uc3n = 0.0;
                    us1n = 0.0;
                    us2n = 0.0;
                    us3n = 0.0;
                    for (kk = 0; kk<bound_count; kk++) {
                        uc1n += compute_uc1n((*low)[ii][kk],(*high)[ii][kk],order);
                        uc2n += compute_uc2n((*low)[ii][kk],(*high)[ii][kk],order);
                        uc3n += compute_uc3n((*low)[ii][kk],(*high)[ii][kk],order);
                        us1n += compute_us1n((*low)[ii][kk],(*high)[ii][kk],order);
                        us2n += compute_us2n((*low)[ii][kk],(*high)[ii][kk],order);
                        us3n += compute_us3n((*low)[ii][kk],(*high)[ii][kk],order);
                        // Compute the partial shadowing before the low boundary
                        if ((*low)[ii][kk] > 0.0) {
                            uc1n += compute_uc1n_shad((*low)[ii][kk] - lamDel, (*low)[ii][kk], order, 0.0, dHall);
                            uc2n += compute_uc2n_shad((*low)[ii][kk] - lamDel, (*low)[ii][kk], order, 0.0, dHall);
                            uc3n += compute_uc3n_shad((*low)[ii][kk] - lamDel, (*low)[ii][kk], order, 0.0, dHall);
                            us1n += compute_us1n_shad((*low)[ii][kk] - lamDel, (*low)[ii][kk], order, 0.0, dHall);
                            us2n += compute_us2n_shad((*low)[ii][kk] - lamDel, (*low)[ii][kk], order, 0.0, dHall);
                            us3n += compute_us3n_shad((*low)[ii][kk] - lamDel, (*low)[ii][kk], order, 0.0, dHall);
                        }
                        // Compute the partial shadowing after the high boundary
                        if ((*high)[ii][kk] < 2.0*M_PI) {
                            uc1n += compute_uc1n_shad((*high)[ii][kk], (*high)[ii][kk]+lamDel, order, 1.0, -dHall);
                            uc2n += compute_uc2n_shad((*high)[ii][kk], (*high)[ii][kk]+lamDel, order, 1.0, -dHall);
                            uc3n += compute_uc3n_shad((*high)[ii][kk], (*high)[ii][kk]+lamDel, order, 1.0, -dHall);
                            us1n += compute_us1n_shad((*high)[ii][kk], (*high)[ii][kk]+lamDel, order, 1.0, -dHall);
                            us2n += compute_us2n_shad((*high)[ii][kk], (*high)[ii][kk]+lamDel, order, 1.0, -dHall);
                            us3n += compute_us3n_shad((*high)[ii][kk], (*high)[ii][kk]+lamDel, order, 1.0, -dHall);
                        }
                    }
                }
                uuc11n = 0.0;
                uuc12n = 0.0;
                uuc22n = 0.0;
                uus11n = 0.0;
                uus12n = 0.0;
                uus22n = 0.0;
                for (kk = 0; kk<bound_count; kk++) {
                    uuc11n += compute_uuc11n((*low)[ii][kk],(*high)[ii][kk],order);
                    uuc12n += compute_uuc12n((*low)[ii][kk],(*high)[ii][kk],order);
                    uuc22n += compute_uuc22n((*low)[ii][kk],(*high)[ii][kk],order);
                    uus11n += compute_uus11n((*low)[ii][kk],(*high)[ii][kk],order);
                    uus12n += compute_uus12n((*low)[ii][kk],(*high)[ii][kk],order);
                    uus22n += compute_uus22n((*low)[ii][kk],(*high)[ii][kk],order);
                    // Compute the partial shadowing before the low boundary
                    if ((*low)[ii][kk] > 0.0) {
                        uuc11n += compute_uuc11n_shad((*low)[ii][kk] - lamDel, (*low)[ii][kk], order, 0.0, dHall);
                        uuc12n += compute_uuc12n_shad((*low)[ii][kk] - lamDel, (*low)[ii][kk], order, 0.0, dHall);
                        uuc22n += compute_uuc22n_shad((*low)[ii][kk] - lamDel, (*low)[ii][kk], order, 0.0, dHall);
                        uus11n += compute_uus11n_shad((*low)[ii][kk] - lamDel, (*low)[ii][kk], order, 0.0, dHall);
                        uus12n += compute_uus12n_shad((*low)[ii][kk] - lamDel, (*low)[ii][kk], order, 0.0, dHall);
                        uus22n += compute_uus22n_shad((*low)[ii][kk] - lamDel, (*low)[ii][kk], order, 0.0, dHall);
                    }
                    // Compute the partial shadowing after the high boundary
                    if ((*high)[ii][kk] < 2.0*M_PI) {
                        uuc11n += compute_uuc11n_shad((*high)[ii][kk], (*high)[ii][kk]+lamDel, order, 1.0, -dHall);
                        uuc12n += compute_uuc12n_shad((*high)[ii][kk], (*high)[ii][kk]+lamDel, order, 1.0, -dHall);
                        uuc22n += compute_uuc22n_shad((*high)[ii][kk], (*high)[ii][kk]+lamDel, order, 1.0, -dHall);
                        uus11n += compute_uus11n_shad((*high)[ii][kk], (*high)[ii][kk]+lamDel, order, 1.0, -dHall);
                        uus12n += compute_uus12n_shad((*high)[ii][kk], (*high)[ii][kk]+lamDel, order, 1.0, -dHall);
                        uus22n += compute_uus22n_shad((*high)[ii][kk], (*high)[ii][kk]+lamDel, order, 1.0, -dHall);
                    }
                }
                
                // Compute nth order integrals
                // Cosine 1 Integral
                Ic_1[0] = cdels*uc1n;
                Ic_1[1] = cdels*uc2n;
                Ic_1[2] = sdels*uc3n;
                
                // Sine 1 Integral
                Is_1[0] = cdels*us1n;
                Is_1[1] = cdels*us2n;
                Is_1[2] = sdels*us3n;
                
                // Cosine 2 Integral
                Ic_2[0][0] = cdels2*uuc11n;
                Ic_2[1][0] = cdels2*uuc12n;
                Ic_2[2][0] = csdels*uc1n;
                Ic_2[0][1] = Ic_2[1][0];
                Ic_2[1][1] = cdels2*uuc22n;
                Ic_2[2][1] = csdels*uc2n;
                Ic_2[0][2] = Ic_2[2][0];
                Ic_2[1][2] = Ic_2[2][1];
                Ic_2[2][2] = sdels2*uc3n;
                
                // Sine 2 Integral
                Is_2[0][0] = cdels2*uus11n;
                Is_2[1][0] = cdels2*uus12n;
                Is_2[2][0] = csdels*us1n;
                Is_2[0][1] = Is_2[1][0];
                Is_2[1][1] = cdels2*uus22n;
                Is_2[2][1] = csdels*us2n;
                Is_2[0][2] = Is_2[2][0];
                Is_2[1][2] = Is_2[2][1];
                Is_2[2][2] = sdels2*us3n;
                
                // Compute Total Coefficients
                A1m = 0.0;
                B1m = 0.0;
                for (jj=0; jj<3; jj++) {
                    A0[jj] = 0.0;
                    A1[jj] = 0.0;
                    A2[jj] = 0.0;
                    B0[jj] = 0.0;
                    B1[jj] = 0.0;
                    B2[jj] = 0.0;
                    for (kk=0; kk<3; kk++) {
                        A0[jj] += Ic_2[jj][kk]*Ar1[kk];
                        A1[jj] += Ar2[jj][kk]*Ic_1[kk];
                        B0[jj] += Is_2[jj][kk]*Ar1[kk];
                        B1[jj] += Ar2[jj][kk]*Is_1[kk];
                        for (ll=0; ll<3; ll++) {
                            A2[jj] += Ar3[jj][kk]*Ic_2[kk][ll]*nhat[ll];
                            B2[jj] += Ar3[jj][kk]*Is_2[kk][ll]*nhat[ll];
                        }
                        if ((numbounces>0) && (jj==0)) {
                            A1m += Ar1[kk]*Ic_1[kk];
                            B1m += Ar1[kk]*Is_1[kk];
                        }
                    }
                    temp[jj] = A0[jj] + a2*A1[jj] + rho*s*A2[jj];
                    A[jj] += temp[jj];
                    
                    temp2[jj] = B0[jj] + a2*B1[jj] + rho*s*B2[jj];
                    B[jj] += temp2[jj];
                    
                    if (numbounces>0) {
                        A[jj] += A1m*(*fmult)[ii][jj];
                        B[jj] += B1m*(*fmult)[ii][jj];
                    }
                }
                C[0] += rc[1]*temp[2] - rc[2]*temp[1];
                C[1] += rc[2]*temp[0] - rc[0]*temp[2];
                C[2] += rc[0]*temp[1] - rc[1]*temp[0];
                
                D[0] += rc[1]*temp2[2] - rc[2]*temp2[1];
                D[1] += rc[2]*temp2[0] - rc[0]*temp2[2];
                D[2] += rc[0]*temp2[1] - rc[1]*temp2[0];
                
                if (numbounces>0) {
                    for (jj=0; jj<3; jj++) {
                        C[jj] += A1m*(*tmult)[ii][jj];
                        D[jj] += B1m*(*tmult)[ii][jj];
                    }
                }

            }
        }
        
        // add specular reflections
        if ((numbounces > 0) && (!(*latSpec)[ii].empty())) {

            iter_d = find((*latSpec)[ii].begin(),(*latSpec)[ii].end(), delta_s);

            while (iter_d != (*latSpec)[ii].end()) {
                // get associated index
                specIndex = iter_d - (*latSpec)[ii].begin();
                
                // Add contribution from matching latitutde specular reflections to current order coefficients
                if (order == 0) {
                    for (jj=0; jj<3; jj++) {
                        A[jj] += (*fSpec)[ii][specIndex][jj]/(2*M_PI);
                        C[jj] += (*tSpec)[ii][specIndex][jj]/(2*M_PI);
                    }
                } else {
                    for (jj=0; jj<3; jj++) {
                        A[jj] += (*fSpec)[ii][specIndex][jj] * cos(order * (*longSpec)[ii][specIndex])/M_PI;
                        B[jj] += (*fSpec)[ii][specIndex][jj] * sin(order * (*longSpec)[ii][specIndex])/M_PI;
                        C[jj] += (*tSpec)[ii][specIndex][jj] * cos(order * (*longSpec)[ii][specIndex])/M_PI;
                        D[jj] += (*tSpec)[ii][specIndex][jj] * sin(order * (*longSpec)[ii][specIndex])/M_PI;
                    }
                }
                
                // update iterator for next search
                iter_d = find(iter_d+1,(*latSpec)[ii].end(), delta_s);
            }
        }
    }
}

double FCoeffs::compute_uc1n(double lbnd, double ubnd, int nint)
{
    double n = (double)nint;
    double fu = (1.0/(n*n-1.0))*(-cos(n*ubnd)*sin(ubnd) + n*cos(ubnd)*sin(n*ubnd));
    double fl = (1.0/(n*n-1.0))*(-cos(n*lbnd)*sin(lbnd) + n*cos(lbnd)*sin(n*lbnd));

    return fu - fl;

}

double FCoeffs::compute_uc2n(double lbnd, double ubnd, int nint)
{
    double n = (double)nint;
    double fu = (1.0/(n*n-1.0))*(cos(n*ubnd)*cos(ubnd) + n*sin(ubnd)*sin(n*ubnd));
    double fl = (1.0/(n*n-1.0))*(cos(n*lbnd)*cos(lbnd) + n*sin(lbnd)*sin(n*lbnd));

    return fu - fl;

}

double FCoeffs::compute_uc3n(double lbnd, double ubnd, int nint)
{
    double n = (double)nint;
    double fu = (1.0/n)*(sin(n*ubnd));
    double fl = (1.0/n)*(sin(n*lbnd));

    return fu - fl;

}

double FCoeffs::compute_uuc11n(double lbnd, double ubnd, int nint)
{
    double n = (double)nint;
    double fu, fl;
    if (nint != 2) {
        fu = (1.0/4.0)*((sin((n-2.0)*ubnd)/(n-2.0)) + (2.0*sin(n*ubnd)/n) + (sin((n+2.0)*ubnd)/(n+2.0)));
        fl = (1.0/4.0)*((sin((n-2.0)*lbnd)/(n-2.0)) + (2.0*sin(n*lbnd)/n) + (sin((n+2.0)*lbnd)/(n+2.0)));
    } else {
        fu = (ubnd/4.0) + (sin(2.0*ubnd)/4.0) + (sin(4.0*ubnd)/16.0);
        fl = (lbnd/4.0) + (sin(2.0*lbnd)/4.0) + (sin(4.0*lbnd)/16.0);
    }

    return fu - fl;

}

double FCoeffs::compute_uuc12n(double lbnd, double ubnd, int nint)
{
    double n = (double)nint;
    double fu, fl;
    if (nint != 2) {
        fu = (1.0/4.0)*((cos((n-2.0)*ubnd)/(n-2.0)) - (cos((n+2.0)*ubnd)/(n+2.0)));
        fl = (1.0/4.0)*((cos((n-2.0)*lbnd)/(n-2.0)) - (cos((n+2.0)*lbnd)/(n+2.0)));
    } else {
        fu = -cos(4.0*ubnd)/16.0;
        fl = -cos(4.0*lbnd)/16.0;
    }

    return fu - fl;

}

double FCoeffs::compute_uuc22n(double lbnd, double ubnd, int nint)
{
    double n = (double)nint;
    double fu, fl;

    if (nint != 2) {
        fu = (1.0/4.0)*((-sin((n-2.0)*ubnd)/(n-2.0)) + (2.0*sin(n*ubnd)/n) - (sin((n+2.0)*ubnd)/(n+2.0)));
        fl = (1.0/4.0)*((-sin((n-2.0)*lbnd)/(n-2.0)) + (2.0*sin(n*lbnd)/n) - (sin((n+2.0)*lbnd)/(n+2.0)));
    } else {
        fu = (-ubnd/4.0) + (sin(2.0*ubnd)/4.0) - (sin(4.0*ubnd)/16.0);
        fl = (-lbnd/4.0) + (sin(2.0*lbnd)/4.0) - (sin(4.0*lbnd)/16.0);
    }

    return fu - fl;

}

double FCoeffs::compute_us1n(double lbnd, double ubnd, int nint)
{
    double n = (double)nint;
    double fu = (1.0/(1.0-n*n))*(sin(n*ubnd)*sin(ubnd) + n*cos(ubnd)*cos(n*ubnd));
    double fl = (1.0/(1.0-n*n))*(sin(n*lbnd)*sin(lbnd) + n*cos(lbnd)*cos(n*lbnd));

    return fu - fl;

}

double FCoeffs::compute_us2n(double lbnd, double ubnd, int nint)
{
    double n = (double)nint;
    double fu = (1.0/(n*n-1.0))*(-n*cos(n*ubnd)*sin(ubnd) + cos(ubnd)*sin(n*ubnd));
    double fl = (1.0/(n*n-1.0))*(-n*cos(n*lbnd)*sin(lbnd) + cos(lbnd)*sin(n*lbnd));

    return fu - fl;

}

double FCoeffs::compute_us3n(double lbnd, double ubnd, int nint)
{
    double n = (double)nint;
    double fu = (1.0/n)*(-cos(n*ubnd));
    double fl = (1.0/n)*(-cos(n*lbnd));

    return fu - fl;

}

double FCoeffs::compute_uus11n(double lbnd, double ubnd, int nint)
{
    double n = (double)nint;
    double fu, fl;

    if (nint != 2) {
        fu = (1.0/4.0)*((-cos((n-2.0)*ubnd)/(n-2.0)) - (2.0*cos(n*ubnd)/n) - (cos((n+2.0)*ubnd)/(n+2.0)));
        fl = (1.0/4.0)*((-cos((n-2.0)*lbnd)/(n-2.0)) - (2.0*cos(n*lbnd)/n) - (cos((n+2.0)*lbnd)/(n+2.0)));
    } else {
        fu = -cos(ubnd)*cos(ubnd)*cos(ubnd)*cos(ubnd)/2.0;
        fl = -cos(lbnd)*cos(lbnd)*cos(lbnd)*cos(lbnd)/2.0;
    }

    return fu - fl;

}

double FCoeffs::compute_uus12n(double lbnd, double ubnd, int nint)
{
    double n = (double)nint;
    double fu, fl;

    if (n != 2) {
        fu = (1.0/4.0)*((sin((n-2.0)*ubnd)/(n-2.0)) - (sin((n+2.0)*ubnd)/(n+2.0)));
        fl = (1.0/4.0)*((sin((n-2.0)*lbnd)/(n-2.0)) - (sin((n+2.0)*lbnd)/(n+2.0)));
    } else {
        fu = (ubnd/4.0) - (sin(4.0*ubnd)/16.0);
        fl = (lbnd/4.0) - (sin(4.0*lbnd)/16.0);
    }

    return fu - fl;

}

double FCoeffs::compute_uus22n(double lbnd, double ubnd, int nint)
{
    double n = (double)nint;
    double fu, fl;

    if (n != 2) {
        fu = (1.0/4.0)*((cos((n-2.0)*ubnd)/(n-2.0)) - (2.0*cos(n*ubnd)/n) + (cos((n+2.0)*ubnd)/(n+2.0)));
        fl = (1.0/4.0)*((cos((n-2.0)*lbnd)/(n-2.0)) - (2.0*cos(n*lbnd)/n) + (cos((n+2.0)*lbnd)/(n+2.0)));
    } else {
        fu = sin(ubnd)*sin(ubnd)*sin(ubnd)*sin(ubnd)/2.0;
        fl = sin(lbnd)*sin(lbnd)*sin(lbnd)*sin(lbnd)/2.0;
    }

    return fu - fl;

}

double FCoeffs::compute_uc1n_shad(double a, double b, int nint, double Ha, double dH)
{
    double n = (double)nint;
    double ans;
    
    if (nint == 1) {
        ans = ((-(dH*cos(2*a)) + dH*cos(2*b) + 2*((a - b)*(a*dH - b*dH - 2*Ha) - Ha*sin(2*a) + (-(a*dH) + b*dH + Ha)*sin(2*b))))/8.0;
    } else {
        ans = ((-(dH*cos(a*(-1 + n))) + dH*cos(b*(-1 + n)) + ((-1 + n)*(-(dH*(-1 + n)*cos(a*(1 + n))) + 2*Ha*(1 + n)*(cos(a*n)*sin(a) - n*cos(a)*sin(a*n))))/pow(1 + n,2) + ((-1 + n)*(-(sin(b)*(2*(-(a*dH) + b*dH + Ha)*(1 + n)*cos(b*n) + dH*(-1 + n)*sin(b*n))) + cos(b)*(dH*(-1 + n)*cos(b*n) + 2*(-(a*dH) + b*dH + Ha)*n*(1 + n)*sin(b*n))))/pow(1 + n,2)))/(2.*pow(-1 + n,2));
    }
    
    return ans;
}

double FCoeffs::compute_uc2n_shad(double a, double b, int nint, double Ha, double dH)
{
    double n = (double)nint;
    double ans;
    
    if (nint == 1) {
        ans = ((2*Ha*cos(2*a) + 2*(a*dH - b*dH - Ha)*cos(2*b) + dH*(-sin(2*a) + sin(2*b))))/8.0;
    } else {
        ans = ((-((Ha*cos(a*(-1 + n)))/(-1 + n)) + ((-(a*dH) + b*dH + Ha)*cos(b*(-1 + n)))/(-1 + n) + (Ha*cos(a*(1 + n)))/(1 + n) - ((-(a*dH) + b*dH + Ha)*cos(b*(1 + n)))/(1 + n) + (dH*sin(a*(-1 + n)))/pow(-1 + n,2) - (dH*sin(b*(-1 + n)))/pow(-1 + n,2) - (dH*sin(a*(1 + n)))/pow(1 + n,2) + (dH*sin(b*(1 + n)))/pow(1 + n,2)))/2.0;
    }
    
    return ans;
}

double FCoeffs::compute_uc3n_shad(double a, double b, int nint, double Ha, double dH)
{
    double n = (double)nint;
    double ans;
    
    ans = ((-(dH*cos(a*n)) + dH*cos(b*n) - Ha*n*sin(a*n) + (-(a*dH) + b*dH + Ha)*n*sin(b*n)))/pow(n,2);
    
    return ans;
}

double FCoeffs::compute_us1n_shad(double a, double b, int nint, double Ha, double dH)
{
    double n = (double)nint;
    double ans;
    
    if (nint == 1) {
        ans = ((2*Ha*cos(2*a) + 2*(a*dH - b*dH - Ha)*cos(2*b) + dH*(-sin(2*a) + sin(2*b))))/8.0;
    } else {
        ans = (((Ha*cos(a*(-1 + n)))/(-1 + n) - ((-(a*dH) + b*dH + Ha)*cos(b*(-1 + n)))/(-1 + n) + (Ha*cos(a*(1 + n)))/(1 + n) - ((-(a*dH) + b*dH + Ha)*cos(b*(1 + n)))/(1 + n) - (dH*sin(a*(-1 + n)))/pow(-1 + n,2) + (dH*sin(b*(-1 + n)))/pow(-1 + n,2) - (dH*sin(a*(1 + n)))/pow(1 + n,2) + (dH*sin(b*(1 + n)))/pow(1 + n,2)))/2.0;
    }
    
    return ans;
}

double FCoeffs::compute_us2n_shad(double a, double b, int nint, double Ha, double dH)
{
    double n = (double)nint;
    double ans;
    
    if (nint == 1) {
        ans = ((dH*cos(2*a) - dH*cos(2*b) + 2*(Ha*(sin(2*a) - sin(2*b)) + (a - b)*(a*dH - b*dH - 2*Ha + dH*sin(2*b)))))/8.0;
    } else {
        ans = ((-(dH*cos(a*(-1 + n))) + dH*cos(b*(-1 + n)) + ((-1 + n)*(dH*(-1 + n)*cos(a*(1 + n)) + 2*Ha*(1 + n)*(n*cos(a*n)*sin(a) - cos(a)*sin(a*n))))/pow(1 + n,2) + ((-1 + n)*(-(sin(b)*(2*(-(a*dH) + b*dH + Ha)*n*(1 + n)*cos(b*n) - dH*(-1 + n)*sin(b*n))) + cos(b)*(-(dH*(-1 + n)*cos(b*n)) + 2*(-(a*dH) + b*dH + Ha)*(1 + n)*sin(b*n))))/pow(1 + n,2)))/(2.*pow(-1 + n,2));
    }
    
    return ans;
}

double FCoeffs::compute_us3n_shad(double a, double b, int nint, double Ha, double dH)
{
    double n = (double)nint;
    double ans;
    
    ans = ((Ha*n*cos(a*n) + (a*dH - b*dH - Ha)*n*cos(b*n) + dH*(-sin(a*n) + sin(b*n))))/pow(n,2);
    
    return ans;
}

double FCoeffs::compute_uuc11n_shad(double a, double b, int nint, double Ha, double dH)
{
    double n = (double)nint;
    double ans;
    
    if (nint == 2) {
        ans = ((8*a*(a*dH - 2*Ha) + 8*b*(-2*a*dH + b*dH + 2*Ha) - 8*dH*cos(2*a) - dH*cos(4*a) + 8*dH*cos(2*b) + dH*cos(4*b) - 4*Ha*(4*sin(2*a) + sin(4*a)) + 4*(-(a*dH) + b*dH + Ha)*(4*sin(2*b) + sin(4*b))))/64.0;
    } else {
        ans = ((-((dH*cos(a*(-2 + n)))/pow(-2 + n,2)) + (dH*cos(b*(-2 + n)))/pow(-2 + n,2) - (2*dH*cos(a*n))/pow(n,2) + (2*dH*cos(b*n))/pow(n,2) - (dH*cos(a*(2 + n)))/pow(2 + n,2) + (dH*cos(b*(2 + n)))/pow(2 + n,2) - (Ha*sin(a*(-2 + n)))/(-2 + n) - (a*dH*sin(b*(-2 + n)))/(-2 + n) + (b*dH*sin(b*(-2 + n)))/(-2 + n) + (Ha*sin(b*(-2 + n)))/(-2 + n) - (2*Ha*sin(a*n))/n - (2*a*dH*sin(b*n))/n + (2*b*dH*sin(b*n))/n + (2*Ha*sin(b*n))/n - (Ha*sin(a*(2 + n)))/(2 + n) - (a*dH*sin(b*(2 + n)))/(2 + n) + (b*dH*sin(b*(2 + n)))/(2 + n) + (Ha*sin(b*(2 + n)))/(2 + n)))/4.0;
    }
    
    return ans;
}

double FCoeffs::compute_uuc12n_shad(double a, double b, int nint, double Ha, double dH)
{
    double n = (double)nint;
    double ans;
    
    if (nint == 2) {
        ans = ((4*Ha*cos(4*a) - 4*(-(a*dH) + b*dH + Ha)*cos(4*b) - dH*sin(4*a) + dH*sin(4*b)))/64.0;
    } else {
        ans = ((-((Ha*cos(a*(-2 + n)))/(-2 + n)) + ((-(a*dH) + b*dH + Ha)*cos(b*(-2 + n)))/(-2 + n) + (Ha*cos(a*(2 + n)))/(2 + n) - ((-(a*dH) + b*dH + Ha)*cos(b*(2 + n)))/(2 + n) + (dH*sin(a*(-2 + n)))/pow(-2 + n,2) - (dH*sin(b*(-2 + n)))/pow(-2 + n,2) - (dH*sin(a*(2 + n)))/pow(2 + n,2) + (dH*sin(b*(2 + n)))/pow(2 + n,2)))/4.0;
    }
    
    return ans;
}

double FCoeffs::compute_uuc22n_shad(double a, double b, int nint, double Ha, double dH)
{
    double n = (double)nint;
    double ans;
    
    if (nint == 2) {
        ans = ((-8*a*(a*dH - 2*Ha) - 8*b*(-2*a*dH + b*dH + 2*Ha) - 8*dH*cos(2*a) + dH*cos(4*a) + 8*dH*cos(2*b) - dH*cos(4*b) + 4*Ha*(-4*sin(2*a) + sin(4*a)) - 4*(-(a*dH) + b*dH + Ha)*(-4*sin(2*b) + sin(4*b))))/64.0;
    } else {
        ans = ((-((a*dH*sin(a*(-2 + n)))/(-2 + n)) + (Ha*sin(a*(-2 + n)))/(-2 + n) + (dH*(cos(a*(-2 + n)) + a*(-2 + n)*sin(a*(-2 + n))))/pow(-2 + n,2) + (a*dH*sin(b*(-2 + n)))/(-2 + n) - (Ha*sin(b*(-2 + n)))/(-2 + n) - (dH*(cos(b*(-2 + n)) + b*(-2 + n)*sin(b*(-2 + n))))/pow(-2 + n,2) + (2*a*dH*sin(a*n))/n - (2*Ha*sin(a*n))/n - (2*dH*(cos(a*n) + a*n*sin(a*n)))/pow(n,2) - (2*a*dH*sin(b*n))/n + (2*Ha*sin(b*n))/n + (2*dH*(cos(b*n) + b*n*sin(b*n)))/pow(n,2) - (a*dH*sin(a*(2 + n)))/(2 + n) + (Ha*sin(a*(2 + n)))/(2 + n) + (dH*(cos(a*(2 + n)) + a*(2 + n)*sin(a*(2 + n))))/pow(2 + n,2) + (a*dH*sin(b*(2 + n)))/(2 + n) - (Ha*sin(b*(2 + n)))/(2 + n) - (dH*(cos(b*(2 + n)) + b*(2 + n)*sin(b*(2 + n))))/pow(2 + n,2)))/4.0;
    }
    
    return ans;
}

double FCoeffs::compute_uus11n_shad(double a, double b, int nint, double Ha, double dH)
{
    double n = (double)nint;
    double ans;
    
    if (nint == 2) {
        ans = ((4*Ha*(4*cos(2*a) + cos(4*a) - 4*cos(2*b) - cos(4*b)) + dH*(4*(a - b)*(4*cos(2*b) + cos(4*b)) - 8*sin(2*a) - sin(4*a) + 8*sin(2*b) + sin(4*b))))/64.0;
    } else {
        ans = (((Ha*cos(a*(-2 + n)))/(-2 + n) - ((-(a*dH) + b*dH + Ha)*cos(b*(-2 + n)))/(-2 + n) + (2*Ha*cos(a*n))/n - (2*(-(a*dH) + b*dH + Ha)*cos(b*n))/n + (Ha*cos(a*(2 + n)))/(2 + n) + (a*dH*cos(b*(2 + n)))/(2 + n) - (b*dH*cos(b*(2 + n)))/(2 + n) - (Ha*cos(b*(2 + n)))/(2 + n) - (dH*sin(a*(-2 + n)))/pow(-2 + n,2) + (dH*sin(b*(-2 + n)))/pow(-2 + n,2) - (2*dH*sin(a*n))/pow(n,2) + (2*dH*sin(b*n))/pow(n,2) - (dH*sin(a*(2 + n)))/pow(2 + n,2) + (dH*sin(b*(2 + n)))/pow(2 + n,2)))/4.0;
    }
    
    return ans;
}

double FCoeffs::compute_uus12n_shad(double a, double b, int nint, double Ha, double dH)
{
    double n = (double)nint;
    double ans;
    
    if (nint == 2) {
        ans = ((dH*cos(4*a) - dH*cos(4*b) + 4*(Ha*(sin(4*a) - sin(4*b)) + (a - b)*(2*a*dH - 2*b*dH - 4*Ha + dH*sin(4*b)))))/64.0;
    } else {
        ans = ((-(dH*cos(a*(-2 + n))) + dH*cos(b*(-2 + n)) + ((-2 + n)*(dH*(-2 + n)*cos(a*(2 + n)) + 4*Ha*(2 + n)*(n*cos(a)*cos(a*n)*sin(a) - pow(cos(a),2)*sin(a*n) + pow(sin(a),2)*sin(a*n))))/pow(2 + n,2) - ((-2 + n)*(dH*(-2 + n)*cos(b*(2 + n)) + 4*(-(a*dH) + b*dH + Ha)*(2 + n)*(n*cos(b)*cos(b*n)*sin(b) - pow(cos(b),2)*sin(b*n) + pow(sin(b),2)*sin(b*n))))/pow(2 + n,2)))/(4.0*pow(-2 + n,2));
    }
    
    return ans;
}

double FCoeffs::compute_uus22n_shad(double a, double b, int nint, double Ha, double dH)
{
    double n = (double)nint;
    double ans;
    
    if (nint == 2) {
        ans = ((4*Ha*(4*cos(2*a) - cos(4*a) - 4*cos(2*b) + cos(4*b)) + dH*(-4*(a - b)*(-4*cos(2*b) + cos(4*b)) - 8*sin(2*a) + sin(4*a) + 8*sin(2*b) - sin(4*b))))/64.0;
    } else {
        ans = ((-((Ha*cos(a*(-2 + n)))/(-2 + n)) + ((-(a*dH) + b*dH + Ha)*cos(b*(-2 + n)))/(-2 + n) + (2*Ha*cos(a*n))/n - (2*(-(a*dH) + b*dH + Ha)*cos(b*n))/n - (Ha*cos(a*(2 + n)))/(2 + n) - (a*dH*cos(b*(2 + n)))/(2 + n) + (b*dH*cos(b*(2 + n)))/(2 + n) + (Ha*cos(b*(2 + n)))/(2 + n) + (dH*sin(a*(-2 + n)))/pow(-2 + n,2) - (dH*sin(b*(-2 + n)))/pow(-2 + n,2) - (2*dH*sin(a*n))/pow(n,2) + (2*dH*sin(b*n))/pow(n,2) + (dH*sin(a*(2 + n)))/pow(2 + n,2) - (dH*sin(b*(2 + n)))/pow(2 + n,2)))/4.0;
    }
    
    return ans;
}

double FCoeffs::compute_u1_shad(double a, double b, double Ha, double dH)
{
    double ans = -(dH*cos(a)) + dH*cos(b) - Ha*sin(a) + (-(a*dH) + b*dH + Ha)*sin(b);
    return ans;
    
}

double FCoeffs::compute_u2_shad(double a, double b, double Ha, double dH)
{
    double ans = Ha*cos(a) + (a*dH - b*dH - Ha)*cos(b) + dH*(-sin(a) + sin(b));
    return ans;
}

double FCoeffs::compute_u3_shad(double a, double b, double Ha, double dH)
{
    double ans = (a - b)*(a*dH - b*dH - 2*Ha)/2.0;
    return ans;
}

double FCoeffs::compute_uu11_shad(double a, double b, double Ha, double dH)
{
    double ans = (-(dH*cos(2*a)) + dH*cos(2*b) + 2*((a - b)*(a*dH - b*dH - 2*Ha) - Ha*sin(2*a) + (-(a*dH) + b*dH + Ha)*sin(2*b)))/8.0;
    return ans;
    
}

double FCoeffs::compute_uu12_shad(double a, double b, double Ha, double dH)
{
    double ans = (2*Ha*cos(2*a) + 2*(a*dH - b*dH - Ha)*cos(2*b) + dH*(-sin(2*a) + sin(2*b)))/8.0;
    return ans;
}

double FCoeffs::compute_uu22_shad(double a, double b, double Ha, double dH)
{
    double ans = (dH*cos(2*a) - dH*cos(2*b) + 2*(Ha*(sin(2*a) - sin(2*b)) + (a - b)*(a*dH - b*dH - 2*Ha + dH*sin(2*b))))/8.0;
    return ans;
}

// Set the order
void FCoeffs::setOrder(int orderIn) {
    order = orderIn;
}

// Get A
double* FCoeffs::getA()
{
    return &A[0];
}

// Get B
double* FCoeffs::getB()
{
    return &B[0];
}

// Get C
double* FCoeffs::getC()
{
    return &C[0];
}

// Get D
double* FCoeffs::getD()
{
    return &D[0];
}