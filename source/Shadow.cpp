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
//  Shadow.cpp
//  Spacecraft_SRP
//
//  Created by Jay McMahon on 5/27/14.
//  Copyright (c) 2014 Jay McMahon. All rights reserved.
//

#include "YORPLib/Shadow.h"
#include <math.h>
namespace YORPLib{

    Shadow::Shadow()
    {
        
    }

    Shadow::Shadow(int numfacets, double delta, double lambda)
    {
        deltaSun = delta;
        lambdaSun = lambda;
        
        percentUnShadowed.resize(numfacets);
    }

    void Shadow::ComputeShadowing(Body* Target, vector<double>* riseLam, vector<double>* setLam)
    {
    // this is where we can use the inView information to simplify what facets to look through when testing to see if shadowed... will have to do comparison tests to see if it makes a significant difference (voxels might be fast enough)
        
    // to start, just ignore inView and test every case
        int numf = Target->getNumFacets();
        int ii, jj;
        int testHit;
        vector<int> dummy;
        
    // Define uhat
    double uhat[3], negu[3], pt[3];//, rmag;
    uhat[0] = cos(lambdaSun)*cos(deltaSun);
    uhat[1] = sin(lambdaSun)*cos(deltaSun);
    uhat[2] = sin(deltaSun);
    for (jj=0; jj<3; jj++) {
        negu[jj] = -1.0*uhat[jj];
    }
    
    double maxDim = Target->getMaxDimVoxGrid();
    
    double* rc;
//    double* rc2;
    double hit_pnt[3];
//    double degBlock = M_PI/180.0; // would be better to input this on a per facet basis, giving a radius of a circle around rc that is within the circle, and then this can be computed to an angle based on facet-to-facet range.
    
    for (ii=0; ii<numf; ii++) {
        // Check to see if Sun is in view based on rise/set information
        percentUnShadowed[ii] = 1;
        if ((*riseLam)[ii]==10.0) {
            // Never sees the sun
            percentUnShadowed[ii] = 0;
        } else if (((*riseLam)[ii]==2*M_PI) && ((*setLam)[ii]==0.0)) {
            // Facet sees the sun at all longitudes
            // Test to see if ray from Sun direction intersects this facet
            // To add inView, do so here; start ray at the current facet and test only for intersection of the inView facets
            if ((*Target->getInView(ii)).size()>0) {
                rc = Target->getFacetRc(ii);
                // This logic starts "at the Sun" and looks toward the facet in quesiton
                // Doesn't need to check for multiple facets as it should hit the center
                for (jj=0; jj<3; jj++) {
                    pt[jj] = rc[jj] + 2.0*maxDim*uhat[jj];
                }
                testHit = Target->ray_intersection(&pt[0], &negu[0], &hit_pnt[0], &dummy);
                if (testHit != ii) {
                    // Facet was shadowed by another facet
                    percentUnShadowed[ii] = 0;
                }
            }
            // Alternative is to use already computed inViewRc data; check each possible facet in view, and if the relative position is within X degrees of uhat, this facet is shadowed
//            if ((*Target->getInViewRc(ii)).size()>0) {
//                rc = Target->getFacetRc(ii);
//                for (int kk=0; kk<(*Target->getInViewRc(ii)).size(); kk++ ){
//                    rc2 = Target->getFacetRc((*Target->getInViewRc(ii))[kk]);
//                    for (jj=0; jj<3; jj++) {
//                        pt[jj] = rc2[jj] - rc[jj];
//                    }
//                    rmag = sqrt(dot(pt,pt));
//                    if (acos(dot(uhat, pt)/rmag) < degBlock) {
//                        // Facet was shadowed by another facet
//                        percentUnShadowed[ii] = 0;
//                        break;
//                    }
//                }
//            }
        } else {
            if (((*riseLam)[ii] < (*setLam)[ii]) && (lambdaSun > (*riseLam)[ii]) && (lambdaSun < (*setLam)[ii])) {
                // Visible region wraps around 2PI, and Sun is outside of it so facet doesn't see Sun
                percentUnShadowed[ii] = 0;
            } else if (((*riseLam)[ii] > (*setLam)[ii]) && ((lambdaSun > (*riseLam)[ii]) || (lambdaSun < (*setLam)[ii]))) {
                // Visible region is continuous, and Sun is outside of it so facet doesn't see Sun
                percentUnShadowed[ii] = 0;
            } else {
                // Facet currently sees the sun
                // Test to see if ray from Sun direction intersects this facet
                // To add inView, do so here; start ray at the current facet and test only for intersection of the inView facets
                if ((*Target->getInView(ii)).size()>0) {
                    rc = Target->getFacetRc(ii);
                    // This logic starts "at the Sun" and looks toward the facet in quesiton
                    // Doesn't need to check for multiple facets as it should hit the center
                    for (jj=0; jj<3; jj++) {
                        pt[jj] = rc[jj] + 2.0*maxDim*uhat[jj];
                    }
                    testHit = Target->ray_intersection(&pt[0], &negu[0], &hit_pnt[0], &dummy);
                    if (testHit != ii) {
                        // Facet was shadowed by another facet
                        percentUnShadowed[ii] = 0;
                    }
                }
                // Alternative is to use already computed inViewRc data; check each possible facet in view, and if the relative position is within X degrees of uhat, this facet is shadowed
//                if ((*Target->getInViewRc(ii)).size()>0) {
//                    rc = Target->getFacetRc(ii);
//                    for (int kk=0; kk<(*Target->getInViewRc(ii)).size(); kk++ ){
//                        rc2 = Target->getFacetRc((*Target->getInViewRc(ii))[kk]);
//                        for (jj=0; jj<3; jj++) {
//                            pt[jj] = rc2[jj] - rc[jj];
//                        }
//                        rmag = sqrt(dot(pt,pt));
//                        if (acos(dot(uhat, pt)/rmag) < degBlock) {
//                            // Facet was shadowed by another facet
//                            percentUnShadowed[ii] = 0;
//                            break;
//                        }
//                    }
//                }
            }
        }
    }
}

int Shadow::getPcntLit(int fnum)
{
    return percentUnShadowed[fnum];
}

double Shadow::dot ( double* a, double* b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    
}
}