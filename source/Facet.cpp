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
//  Facet.cpp
//  Spacecraft_SRP
//
//  Created by Jay McMahon on 5/19/14.
//  Copyright (c) 2014 Jay McMahon. All rights reserved.
//

#include "YORPLib/Facet.h"
#include <math.h>
#include <iostream>

namespace YORPLib{


// Overridden default constructor
    Facet::Facet()
    {
        vert_nums[0] = 0;
        vert_nums[1] = 0;
        vert_nums[2] = 0;
        rho = 0.0;
        s = 0.0;
    }

// Preferred constructor
    Facet::Facet(int* verts, vector<double>* vertices, double rhoSet, double sSet)
    {
    // Note that verts numbers start at 1, not 0! may want to change this on input...
        vert_nums[0] = 3*(verts[0]-1);
        vert_nums[1] = 3*(verts[1]-1);
        vert_nums[2] = 3*(verts[2]-1);

        rho = rhoSet;
        
        s = sSet;
        
    // Get the vertex coordinates
        double* v1;
        double* v2;
        double* v3;
        
    // std::cout << "before computeGeometry" << std::endl;


    // std::cout << verts[0] << " " << verts[1] <<" "<< verts[2] << std::endl;
        
    // std::cout << vert_nums[0] << " " << vert_nums[1] <<" "<< vert_nums[2] << std::endl;



        v1 = &vertices->at(vert_nums[0]);
        v2 = &vertices->at(vert_nums[1]);
        v3 = &vertices->at(vert_nums[2]);

    // std::cout << "computeGeometry" << std::endl;
    // Call member functions to fill remaining values
        computeGeometry(&v1[0], &v2[0], &v3[0], &normal[0], &r_c[0], &area);
        
    }

// Private function to compute geometry
    void Facet::computeGeometry(double* v1, double* v2, double* v3, double* normal, double* r_c, double* area)
    {
        int ii;
        
    // Compute r_c
        for (ii=0; ii<3; ii++) {
            r_c[ii] = (v1[ii]+v2[ii]+v3[ii])/3;
        }
        
    // Compute edge vectors
        double e1[3];
        double e2[3];
        double e3[3];
        for (ii=0; ii<3; ii++) {
            e1[ii] = v2[ii] - v1[ii];
            e2[ii] = v3[ii] - v2[ii];
            e3[ii] = v1[ii] - v3[ii];
        }
        
        
    // Compute normal
    // This could benefit from a cross product and norm function
        normal[0] = e1[1]*e2[2] - e1[2]*e2[1];
        normal[1] = e1[2]*e2[0] - e1[0]*e2[2];
        normal[2] = e1[0]*e2[1] - e1[1]*e2[0];
        double nmag;
        nmag = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
        for (ii=0; ii<3; ii++) {
            normal[ii] = normal[ii]/nmag;
        }
        
    // Compute area
        double a, b, c, perim_2;
        a = sqrt(e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2]);
        b = sqrt(e2[0]*e2[0] + e2[1]*e2[1] + e2[2]*e2[2]);
        c = sqrt(e3[0]*e3[0] + e3[1]*e3[1] + e3[2]*e3[2]);
        perim_2 = (a+b+c)/2;
        
        *area = sqrt(perim_2*(perim_2-a)*(perim_2-b)*(perim_2-c));
    }

// Get one of the vertex numbers
    int Facet::getVert(int vnum)
    {
        return vert_nums[vnum];
    }

// Get reflectivity
    double Facet::getRho()
    {
        return rho;
    }

// Get specularity
    double Facet::getS()
    {
        return s;
    }

// Get area
    double Facet::getArea()
    {
        return area;
    }

// Get facet center
    double* Facet::getRc()
    {
        return &r_c[0];
    }

// Get normal
    double* Facet::getNormal()
    {
        return &normal[0];
    }

// Get a2 parameter
    double Facet::geta2()
    {
        return a2;
    }

// Set Lambertian model optical properties needed to compute coefficients
    void Facet::computeOptical()
    {
    // Lambertian
        double B = 2.0/3.0;
        
        a2 = B*rho*(1.0-s) + B*(1.0-rho);
        
        int ii, jj;
        
        for (ii=0; ii<3; ii++) {
            Ar1[ii] = -area/M_PI*normal[ii];
            Ar1_2[ii] = Ar1[ii]/2.0;
            for (jj=0; jj<3; jj++) {
                nn[ii][jj] = normal[ii]*normal[jj];
                Ar2[ii][jj] = -area/M_PI*nn[ii][jj];
                Ar2_2[ii][jj] = Ar2[ii][jj]/2.0;
                if (ii == jj) {
                    Ar3[ii][jj] = -area/M_PI*(2.0*nn[ii][jj] - 1.0);
                } else {
                    Ar3[ii][jj] = -area/M_PI*(2.0*nn[ii][jj]);
                }
                Ar3_2[ii][jj] = Ar3[ii][jj]/2.0;
            }
        }
    }

// update Rho and S, generally for the case when they weren't originally input
    void Facet::updateRhoS(double rhoSet, double sSet)
    {
        rho = rhoSet;
        
        s = sSet;
    }

/*void Facet::getOpticalFourier(double* nnOut, double* Ar1Out, double* Ar2Out, double* Ar3Out, double* Ar1_2Out, double* Ar2_2Out, double* Ar3_2Out)
{
    nnOut = &nn[0][0];
    Ar1Out = &Ar1[0];
    Ar2Out = &Ar2[0][0];
    Ar3Out = &Ar3[0][0];
    Ar1_2Out = &Ar1_2[0];
    Ar2_2Out = &Ar2_2[0][0];
    Ar3_2Out = &Ar3_2[0][0];
}*/

    void Facet::getOpticalFourier(double nnOut[3][3], double Ar1Out[3], double Ar2Out[3][3], double Ar3Out[3][3], double Ar1_2Out[3], double Ar2_2Out[3][3], double Ar3_2Out[3][3])
    {
        int ii, jj;
        
        for (ii=0; ii<3; ii++) {
            Ar1Out[ii] = Ar1[ii];
            Ar1_2Out[ii] = Ar1_2[ii];
            for (jj=0; jj<3; jj++) {
                nnOut[ii][jj] = nn[ii][jj];
                Ar2Out[ii][jj] = Ar2[ii][jj];
                Ar3Out[ii][jj] = Ar3[ii][jj];
                Ar2_2Out[ii][jj] = Ar2_2[ii][jj];
                Ar3_2Out[ii][jj] = Ar3_2[ii][jj];
            }
        }
    }
}