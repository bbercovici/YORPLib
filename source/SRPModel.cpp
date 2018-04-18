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
//  SRPModel.cpp
//  Spacecraft_SRP
//
//  Created by Jay McMahon on 5/27/14.
//  Copyright (c) 2014 Jay McMahon. All rights reserved.
//

#include "SRPModel.h"
#include <math.h>
namespace YORPLib{

    SRPModel::SRPModel(){
        return;
    }

    SRPModel::SRPModel(double lambdaDel, double deltaDel, int MaxFourier, Body* bodyExist, int bounces, int numRefine){

    // Assign body info
        Spacecraft = *bodyExist;
        
    // Assign max order for Fourier coefficient
        FourierOrder = MaxFourier;
        
    // Assign number of bounces
        numBounces = bounces;
        
    // Assign how many refinement steps to take
        numShadRefine = numRefine;
        
    // Set delta_s and lambda_s lists
        int ii, jj;
        int dsnum = 180/deltaDel + 1;
        int lsnum = 360/lambdaDel + 1;
        deltaSunList.resize(dsnum);
        lambdaSunList.resize(lsnum);
        riseLambda.resize(dsnum);
        setLambda.resize(dsnum);
        
        for (ii=0; ii<dsnum; ii++) {
            deltaSunList[ii] = M_PI_2 - ii*deltaDel*M_PI/180.0;
        }
    // If deltaDel wasn't a perfect divisor of 180, will make the last entry off - correct it here
        if (deltaSunList.back() != -1.0*M_PI_2) {
            deltaSunList.back() = -1.0*M_PI_2;
        }
        
        lambdaDel = lambdaDel*M_PI/180.0;
        for (ii=0; ii < lsnum; ii++) {
            lambdaSunList[ii] = ii*lambdaDel;
        }
    // If lambdaDel wasn't a perfect divisor of 360, will make the last entry off - correct it here by adjusting lambdaDel to perfectly divide 360 by the same number of divisions as before
        if (lambdaSunList.back() != 2.0*M_PI) {
//        lambdaSunList.back() = 2.0*M_PI;
            lambdaDel = 360.0/((double)lsnum -1.0);

            lambdaDel = lambdaDel*M_PI/180.0;
            for (ii=0; ii < lsnum; ii++) {
                lambdaSunList[ii] = ii*lambdaDel;
            }
        }
        
    // Fill out the rise/set information
        int numf = Spacecraft.getNumFacets();
        vector<double> row(numf);
        
        for (ii=0; ii<dsnum; ii++) {
            riseLambda[ii] = row;
            setLambda[ii] = row;
            for (jj=0; jj<numf; jj++) {
                computeRiseSet(Spacecraft.getFacetNormal(jj), deltaSunList[ii], &riseLambda[ii][jj], &setLambda[ii][jj]);
            }
        }
        
    // Compute effect of secondary bounces
        if (numBounces > 0) {
            Spacecraft.setViewRc();
            computeMulti();
        // Specular reflections
        // if flag is set? or always
            specularMulti(dsnum, lsnum);
        }
        
    // Create coefficients array and fill
        vector<FCoeffs> Frow(FourierOrder);
        Coefficients.resize(dsnum);
        shadowing.resize(lsnum);
        shadowBoundsLow.resize(numf);
        shadowBoundsHigh.resize(numf);

        for (ii=0; ii<dsnum; ii++) {
            
        // Compute shadowing information for this delta_s
        #pragma omp parallel for
            for (jj=0; jj<lsnum; ++jj) {
                Shadow initShad(numf, deltaSunList[ii], lambdaSunList[jj]);
                shadowing[jj] = initShad;
                shadowing[jj].ComputeShadowing(&Spacecraft, &riseLambda[ii], &setLambda[ii]);
            }
            
        // Compute shadowing bounds and refine if desired
            computeShadowBounds(ii, numf, lsnum, lambdaDel);
            
            Coefficients[ii] = Frow;
        #pragma omp parallel for
            for (jj=0; jj<FourierOrder; ++jj) {

                Coefficients[ii][jj].setOrder(jj);
                Coefficients[ii][jj].computeCoeffs(deltaSunList[ii], &Spacecraft, &shadowBoundsLow, &shadowBoundsHigh, lambdaDel/pow(2,numShadRefine), numBounces, &fMult, &tMult, &latSpec, &longSpec, &fSpec, &tSpec);
            }
        }
        
    }

    void SRPModel::computeRiseSet(const double* normal, const double latitude, double* rise, double* set)
    {
    // Originally programmed by Dan Scheeres
    // Jay changed so the second input is latitude instead of zenith, and so the
    // output is modded by 2pi.
    // If facet is never illuminated, both entries = 10
        
        double zenith, czs, szs, czo, szo, nLong, arg, dellon;
        
    // Input is latitude [-PI/2 PI/2], zenith is needed [0 PI] measured from North Pole
        zenith = M_PI_2 - latitude;
        
    czs = cos(zenith);  // Zb-component of uhat
    szs = sin(zenith);  // Xb/Yb-component of uhat
    
    czo = normal[2];        // Zb-component of nhat
    szo = sqrt(1-czo*czo);    // Xb/Yb-component of nhat
    
    // Test dot(uhat,nhat) if uhat is at the same longitude as nhat
    if ((czs*czo + szs*szo) <= 1e-15){
        // Never illuminated
        *rise = 10.0;
        *set = 10.0;
    }
    // Test dot(uhat,nhat) if uhat is at the longitude 180 degrees from nhat
    else if ((czs*czo - szs*szo) >= 0){
        // Always illuminated
        *rise = 2*M_PI;
        *set = 0.0;
    }
    else{
        // Sometimes illuminated
        nLong = atan2(normal[1],normal[0]);
        arg = -(czo*czs)/(szo*szs);
        dellon = acos(arg);
        *rise = fmod((nLong + dellon), (2*M_PI));
        if (*rise < 0.0) *rise += 2*M_PI;
        *set = fmod((nLong - dellon), (2*M_PI));
        if (*set < 0.0) *set += 2*M_PI;
    }
    return;
}

void SRPModel::writeSRPCoeffs(int deltaWrite)
{
    int ii;
    double* Ap;
    double* Bp;
    
    // cout << deltaSunList[deltaWrite]*180.0/M_PI << " ";
    Ap = Coefficients[deltaWrite][0].getA();
    // cout << *Ap << " " << *(Ap+1) << " " << *(Ap+2) << " ";
    for (ii=1; ii<FourierOrder; ii++) {
        Ap = Coefficients[deltaWrite][ii].getA();
        Bp = Coefficients[deltaWrite][ii].getB();
        // cout << *Ap << " " << *(Ap+1) << " " << *(Ap+2) << " ";
        // cout << *Bp << " " << *(Bp+1) << " " << *(Bp+2) << " ";
    }
    // cout << "\n";
    
    return;
}

void SRPModel::writeSRPCoeffsFile(string outputFileBaseName, int numDelta)
{
    string FfileName = outputFileBaseName + "/F.txt" ;
    string MfileName = outputFileBaseName + "/M.txt" ;
    
    ofstream Ffile(FfileName, std::ofstream::out);
    ofstream Mfile(MfileName, std::ofstream::out);
    
    int dd, ii;
    double* Ap;
    double* Bp;
    double* Cp;
    double* Dp;
    
    for (dd = 0; dd<numDelta; dd++) {
        // Print force coefficients
        Ffile << deltaSunList[dd]*180.0/M_PI << " ";
        Ap = Coefficients[dd][0].getA();
        Ffile << *Ap << " " << *(Ap+1) << " " << *(Ap+2) << " ";
        for (ii=1; ii<FourierOrder; ii++) {
            Ap = Coefficients[dd][ii].getA();
            Bp = Coefficients[dd][ii].getB();
            Ffile << *Ap << " " << *(Ap+1) << " " << *(Ap+2) << " ";
            Ffile << *Bp << " " << *(Bp+1) << " " << *(Bp+2) << " ";
        }
        Ffile << "\n";
        
        // Print moment coefficients
        Mfile << deltaSunList[dd]*180.0/M_PI << " ";
        Cp = Coefficients[dd][0].getC();
        Mfile << *Cp << " " << *(Cp+1) << " " << *(Cp+2) << " ";
        for (ii=1; ii<FourierOrder; ii++) {
            Cp = Coefficients[dd][ii].getC();
            Dp = Coefficients[dd][ii].getD();
            Mfile << *Cp << " " << *(Cp+1) << " " << *(Cp+2) << " ";
            Mfile << *Dp << " " << *(Dp+1) << " " << *(Dp+2) << " ";
        }
        Mfile << "\n";
    }
    
    Ffile.close();
    Mfile.close();
    
}

void SRPModel::computeMulti()
{
    // run through each facet and just use inView information?
    // or compute adjacency matrix like in previous code from inView?
    // Easiest way is to re-use, so just convert inView into adj matrix here
    
    int ii, jj, kk, multiCount;
    int numF = Spacecraft.getNumFacets();
    vector<int>* view;
    vector<int> adjRows;
    vector<int> adjCols;
    vector < vector < double > > relVec;
    vector < double > relVecR;
    vector < vector < double > > relVecHat;
    vector<double> row(3);
    double* rc1;
    double* rc2;
    double temp;
    
    // Get the count and required data for facets that see eachother
    multiCount = 0;
    for (ii=0; ii<numF; ii++) {
        view = Spacecraft.getInViewRc(ii); // need to use getInViewRc here since this logic only uses rc information
        for (jj=0; jj < view->size(); jj++) {
            if (view->at(jj) > ii) {
                multiCount++;
            }
        }
    }
    adjRows.resize(multiCount);
    adjCols.resize(multiCount);
    relVec.resize(multiCount);
    relVecR.resize(multiCount);
    relVecHat.resize(multiCount);
    int currentCount = 0;
    for (ii=0; ii<numF; ii++) {
        view = Spacecraft.getInViewRc(ii); // need to use getInViewRc here since this logic only uses rc information
        rc1 = Spacecraft.getFacetRc(ii);
        for (jj=0; jj < view->size(); jj++) {
            if (view->at(jj) > ii) {
                adjRows[currentCount] = ii;
                adjCols[currentCount] = view->at(jj);
                rc2 = Spacecraft.getFacetRc(view->at(jj));
                for (kk=0; kk<3; kk++) {
                    row[kk] = rc2[kk] - rc1[kk];
                }
                relVec[currentCount] = row;
                temp = sqrt(row[0]*row[0]+row[1]*row[1]+row[2]*row[2]);
                relVecR[currentCount] = temp;
                for (kk=0; kk<3; kk++) {
                    row[kk] /= temp;
                }
                relVecHat[currentCount] = row;
                currentCount++;
            }
        }
    }
    
    // Compute reflectivity matrix entries
    vector<int> reflMatRows(2*multiCount);
    vector<int> reflMatCols(2*multiCount);
    vector<double> reflMatVals(2*multiCount);
    double theta_mid, deltaAng, pwrBounces;
    for (ii=0; ii<multiCount; ii++) {
        
        // Compute fraction from facet1 to facet2
        theta_mid = acos(dot(&relVecHat[ii][0], Spacecraft.getFacetNormal(adjRows[ii])));
        deltaAng = sqrt(Spacecraft.getFacetArea(adjCols[ii]) * (-1.0) * dot(&relVecHat[ii][0], Spacecraft.getFacetNormal(adjCols[ii]))/(relVecR[ii]*relVecR[ii]));
        reflMatRows[2*ii] = adjRows[ii];
        reflMatCols[2*ii] = adjCols[ii];
        reflMatVals[2*ii] = (deltaAng/(2*M_PI))*sin(2.0*theta_mid)*sin(deltaAng);
        
        /*        if (count > 4154){
         printf("i=%d\n",i);
         printf("theta = %.5g\n",theta_mid);
         printf("area col = %.5g\n",areas[adjCols[i]]);
         printf("r = %.5g\n",relVecR[i]);
         printf("rhat = %.5g %.5g %.5g\n",relVecHat[3*i],relVecHat[3*i+1],relVecHat[3*i+2]);
         printf("norm col = %.5g %.5g %.5g\n",normals[3*adjCols[i]],normals[3*adjCols[i]+1],normals[3*adjCols[i]+2]);
         printf("dAng = %.5g\n",deltaAng);
         }
         */
        
        // Compute fraction from facet2 to facet1
        theta_mid = acos(-1.0 * dot(&relVecHat[ii][0], Spacecraft.getFacetNormal(adjCols[ii])));
        deltaAng = sqrt(Spacecraft.getFacetArea(adjRows[ii]) * dot(&relVecHat[ii][0], Spacecraft.getFacetNormal(adjRows[ii]))/(relVecR[ii]*relVecR[ii]));
        reflMatRows[2*ii + 1] = adjCols[ii];
        reflMatCols[2*ii + 1] = adjRows[ii];
        reflMatVals[2*ii + 1] = (deltaAng/(2*M_PI))*sin(2.0*theta_mid)*sin(deltaAng);
        
        
        // if only modeling 1 or 2 bounces, use the series power bounce to increase fidelity
        if (numBounces < 3) {
            // Compute power series for repeated reflections
            pwrBounces = 1/(1-reflMatVals[2*ii]*reflMatVals[2*ii + 1]);
            
            // Put power series into reflectivity matrix
            reflMatVals[2*ii] = pwrBounces*reflMatVals[2*ii];
            reflMatVals[2*ii + 1] = pwrBounces*reflMatVals[2*ii + 1];
        }
        
    }
    
    //    printf("Refl[2 to 19] = %.5g\n",reflMatVals[0]);
    //    printf("Refl[19 to 2] = %.5g\n",reflMatVals[1]);
    
    // Compute all paths up to N bounces
    // Currently, for the sake of getting it working, I'm only going to do 1 bounce
    double tempvec[3], tempvec2[3];
    double Lambertian = 2.0/3.0;
    double *tempn;
    double *temprc;
    
    // initialize all entries with zeros
    fMult.resize(numF);
    tMult.resize(numF);
    row[0] = 0.0;
    row[1] = 0.0;
    row[2] = 0.0;
    for (ii=0; ii<numF; ii++) {
        fMult[ii] = row;
        tMult[ii] = row;
    }
    
    for (ii=0; ii<multiCount; ii++) {
        // Upper triangular component
        // Compute the Lambertian component out of target facet
        //        printf("count %d facets: %d %d\n",i,reflMatRows[2*i],reflMatRows[2*i+1]);
        tempn = Spacecraft.getFacetNormal(reflMatCols[2*ii]);
        for (jj=0; jj<3; jj++) {
            tempvec[jj] = Lambertian*tempn[jj]*reflMatVals[2*ii];
            fMult[reflMatRows[2*ii]][jj] += tempvec[jj];
        }
        
        // Compute the torque component
//        cross(&rc[3*reflMatCols[2*i]], &temp[0], &temp2[0]);
        temprc = Spacecraft.getFacetRc(reflMatCols[2*ii]);
        cross(temprc, &tempvec[0], &tempvec2[0]);
        for (jj=0; jj<3; jj++) {
            tMult[reflMatRows[2*ii]][jj] += tempvec2[jj];
        }
        
        // Compute the directed component along relative vector line
        // Note the +/- signs when adding tempvecs are correct b/c in FCoeffs, there is a negative sign
        // from using Ar1...
        if (reflMatRows[2*ii]<reflMatCols[2*ii]) { // check the inequality sign here...
            // Compute the force component
            for (jj=0; jj<3; jj++) {
                tempvec[jj] = reflMatVals[2*ii]*relVecHat[ii][jj];
                fMult[reflMatRows[2*ii]][jj] -= tempvec[jj];
            }
            
            // Compute the torque component
            cross(temprc, &tempvec[0], &tempvec2[0]);
            for (jj=0; jj<3; jj++) {
                tMult[reflMatRows[2*ii]][jj] -= tempvec2[jj];
            }
        } else {
            // Compute the force component
            for (jj=0; jj<3; jj++) {
                tempvec[jj] = reflMatVals[2*ii]*relVecHat[ii][jj];
                fMult[reflMatRows[2*ii]][jj] += tempvec[jj];
            }
            
            // Compute the torque component
            cross(temprc, &tempvec[0], &tempvec2[0]);
            for (jj=0; jj<3; jj++) {
                tMult[reflMatRows[2*ii]][jj] += tempvec2[jj];
            }
        }
        
        //        printf("fmult[%d] = %.5g %.5g %.5g\n",reflMatRows[2*i],fMult[3*reflMatRows[2*i]],fMult[3*reflMatRows[2*i]+1],fMult[3*reflMatRows[2*i]+2]);
        //        printf("tmult[%d] = %.5g %.5g %.5g\n",reflMatRows[2*i],tMult[3*reflMatRows[2*i]],tMult[3*reflMatRows[2*i]+1],tMult[3*reflMatRows[2*i]+2]);
        
        // Lower triangular component
        // Compute the Lambertian component out of target facet
        tempn = Spacecraft.getFacetNormal(reflMatCols[2*ii+1]);
        for (jj=0; jj<3; jj++) {
            tempvec[jj] = Lambertian*tempn[jj]*reflMatVals[2*ii+1];
            fMult[reflMatRows[2*ii+1]][jj] += tempvec[jj];
        }
        
        // Compute the torque component
        temprc = Spacecraft.getFacetRc(reflMatCols[2*ii+1]);
        cross(temprc, &tempvec[0], &tempvec2[0]);
        for (jj=0; jj<3; jj++) {
            tMult[reflMatRows[2*ii+1]][jj] += tempvec2[jj];
        }
        
        // Compute the directed component along relative vector line
        if (reflMatRows[2*ii+1]<reflMatCols[2*ii+1]) { // check the inequality sign here...
            // Compute the force component
            for (jj=0; jj<3; jj++) {
                tempvec[jj] = reflMatVals[2*ii+1]*relVecHat[ii][jj];
                fMult[reflMatRows[2*ii+1]][jj] -= tempvec[jj];
            }
            
            // Compute the torque component
            cross(temprc, &tempvec[0], &tempvec2[0]);
            for (jj=0; jj<3; jj++) {
                tMult[reflMatRows[2*ii+1]][jj] -= tempvec2[jj];
            }
        } else {
            // Compute the force component
            for (jj=0; jj<3; jj++) {
                tempvec[jj] = reflMatVals[2*ii+1]*relVecHat[ii][jj];
                fMult[reflMatRows[2*ii+1]][jj] += tempvec[jj];
            }
            
            // Compute the torque component
            cross(temprc, &tempvec[0], &tempvec2[0]);
            for (jj=0; jj<3; jj++) {
                tMult[reflMatRows[2*ii+1]][jj] += tempvec2[jj];
            }
        }
        
        //        printf("fmult[%d] = %.5g %.5g %.5g\n",reflMatRows[2*i+1],fMult[3*reflMatRows[2*i+1]],fMult[3*reflMatRows[2*i+1]+1],fMult[3*reflMatRows[2*i+1]+2]);
        //        printf("tmult[%d] = %.5g %.5g %.5g\n",reflMatRows[2*i+1],tMult[3*reflMatRows[2*i+1]],tMult[3*reflMatRows[2*i+1]+1],tMult[3*reflMatRows[2*i+1]+2]);
    }
    
    if (numBounces>1) {
        
        // Make a copy of the single bounce fmult and tmult
        vector < vector < double > > fMultLast(fMult);
        vector < vector < double > > tMultLast(tMult);
        
        // Loop over number of bounces requested
        int bounce = 2;
        vector < vector < double > > fMultNew;
        vector < vector < double > > tMultNew;
        // initialize all entries with zeros
        fMultNew.resize(numF);
        tMultNew.resize(numF);
        row[0] = 0.0;
        row[1] = 0.0;
        row[2] = 0.0;
        for (ii=0; ii<numF; ii++) {
            fMultNew[ii] = row;
            tMultNew[ii] = row;
        }
        
        
        while (bounce <= numBounces) {
            // compute the contribution of the new bounce
            for (ii=0; ii<2*multiCount; ii++) {
                for (jj=0; jj<3; jj++) {
                    fMultNew[reflMatRows[ii]][jj] += reflMatVals[ii]*fMultLast[reflMatCols[ii]][jj];
                    tMultNew[reflMatRows[ii]][jj] += reflMatVals[ii]*tMultLast[reflMatCols[ii]][jj];
                }
            }
            
            // add new contribution to totals
            // DO I NEED THE LOOP OVER JJ HERE? OR CAN I JUST SAY fMult[ii] += fMultNew[ii]; ??
            for (ii=0; ii<numF; ii++) {
                for (jj=0; jj<3; jj++) {
                    fMult[ii][jj] += fMultNew[ii][jj];
                    tMult[ii][jj] += tMultNew[ii][jj];
                }
            }
            
            // copy new contribution to last contribution and reset new values to zero
            for (ii=0; ii<numF; ii++) {
                fMultLast[ii] = fMultNew[ii];
                fMultNew[ii] = row;
                tMultLast[ii] = tMultNew[ii];
                tMultNew[ii] = row;
            }
            
            bounce++;
        }
        
    }
    
    // Print results
//    fprintf(fp,"Multi Reflection facetnum fmult tmult\n");
//    for (i=0;i<numfac;i++){
//        fprintf(fp,"%d %.5g %.5g %.5g %.5g %.5g %.5g\n",i+1,fMult[3*i],fMult[3*i + 1],fMult[3*i+2],tMult[3*i],tMult[3*i + 1],tMult[3*i+2]);
//    }
    
    return;
}

double SRPModel::dot ( double* a, double* b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    
}

void SRPModel::cross(double* a, double* b, double* c)
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
    return;
}

void SRPModel::specularMulti(const int dsnum, const int lsnum)
{
    
    int ii, kk, numF, test, bounce, mm, nn;
    int lstemp;
    vector<int>* view;
    double uhat[3], vmag, pt[3];
    double hit_pnt[3], lightSource[3];
    double ftemp[3], ttemp[3], fi[3], mi[3];
    double uhat_new[3], scale_light;
    double* rc1;
    double *tempn;
    int currF, currSpecCase;
    double currDelS, currLamS;
    Facet tempF, tempF2;
    vector<double> rowd;
    vector < vector < double > > rowvd;
    vector<int> extraFacets;
    
    numF = Spacecraft.getNumFacets();
    
    fSpec.resize(numF);
    tSpec.resize(numF);
    latSpec.resize(numF);
    longSpec.resize(numF);
    
    for (ii=0; ii<numF; ii++) {
        view = Spacecraft.getInViewRc(ii);
        // If running into memory issues with specular, use reserve here on a loop through each component of the following multidimensional vectors
        // Actually this is probably a place to use vectors of pointers?
        latSpec[ii] = rowd;
        longSpec[ii] = rowd;
        fSpec[ii] = rowvd;
        tSpec[ii] = rowvd;
        currSpecCase = 0;
        currF = ii;
        tempF = Spacecraft.getFacet(currF);
        
        if (!view->empty() && (tempF.getS()>0.0)) {
            
            // lat/long loop
            for (mm=0; mm<dsnum; mm++) {
                currDelS = deltaSunList[mm];
                if ((mm==0) || (mm == dsnum-1)) {
                    lstemp = 1;
                } else {
                    lstemp = lsnum;
                }
                for (nn=0; nn<lstemp; nn++) {
                    currLamS = lambdaSunList[nn];
                    // Check rise/set conditions
                    if ((riseLambda[mm][ii]==10.0) || ((riseLambda[mm][ii] < setLambda[mm][ii]) && (currLamS > riseLambda[mm][ii]) && (currLamS < setLambda[mm][ii])) || ((riseLambda[mm][ii] > setLambda[mm][ii]) && ((currLamS > riseLambda[mm][ii]) || (currLamS < setLambda[mm][ii]))) ) {
                        continue;
                    }
                    lightSource[0] = cos(currDelS)*cos(currLamS);
                    lightSource[1] = cos(currDelS)*sin(currLamS);
                    lightSource[2] = sin(currDelS);
                    for (kk=0; kk<3; kk++) {
                        ftemp[kk] = 0.0;
                        ttemp[kk] = 0.0;
                    }
                    
                    bounce = 1;
                    scale_light = tempF.getRho() * tempF.getS();
                    rc1 = Spacecraft.getFacetRc(currF);
                    while (bounce <= numBounces) {
                        // determine direction of ray to search
                        tempn = Spacecraft.getFacetNormal(currF);
                        vmag = dot(&lightSource[0], tempn);
                        // skip if dot product is too close to zero - neglible bouncing power and gives bad numerical results
                        if (vmag < 1e-15) {
                            test = -1;
                        } else {
                            for (kk=0; kk<3; kk++) {
                                uhat[kk] = -lightSource[kk] + 2*vmag*tempn[kk];
                            }
                            vmag = sqrt(dot(&uhat[0], &uhat[0]));
                            for (kk=0; kk<3; kk++) {
                                uhat[kk] = uhat[kk]/vmag;
                                pt[kk] = rc1[kk] + 0.001*uhat[kk]; // move starting point 0.1% of the distance along the connecting line so the ray tracing check doesn't just return the starting facet
                            }
                            test = Spacecraft.ray_intersection(&pt[0], &uhat[0], &hit_pnt[0], &extraFacets);
                        }
                        // If there is a hit, and it hits from the front of the facet
                        if ((test != -1) && (dot(&uhat[0], Spacecraft.getFacetNormal(test))<-1e-15)) {
                            // Figure out the force/moment updates
                            // should just be a multiplier/factor I can plug into my coefficient calculations
                            // associate with each facet that receives a subsequent bounce
                            // computing f_i and m_i here; units of [area^2]
                            for (kk=0; kk<3; kk++) {
                                uhat_new[kk] = -uhat[kk];
                            }
                            tempF = Spacecraft.getFacet(test);
                            compute_fi(&uhat_new[0], tempF, &fi[0]);
                            rc1 = Spacecraft.getFacetRc(test);
                            cross(rc1, &fi[0], &mi[0]);
                            for (kk=0; kk<3; kk++) {
                                ftemp[kk] += scale_light*fi[kk];
                                ttemp[kk] += scale_light*mi[kk];
                            }
                            
                            if (!extraFacets.empty()) {
                                // Go through any other facets hit at the same time
                                // at this point, we only compute the force from these extra facets, but bounce is only computed from the first hit facet
                                // I think what I want to do is an average of the slopes of the facets hit to compute the next bounce... also compute an average for updating lightSource? not sure I can do this with the dependence on currF though...
                                for (int qq=0; qq<extraFacets.size(); qq++) {
                                    
                                    tempF2 = Spacecraft.getFacet(extraFacets[qq]);
                                    compute_fi(&uhat_new[0], tempF2, &fi[0]);
                                    rc1 = Spacecraft.getFacetRc(extraFacets[qq]);
                                    cross(rc1, &fi[0], &mi[0]);
                                    for (kk=0; kk<3; kk++) {
                                        ftemp[kk] += scale_light*fi[kk];
                                        ttemp[kk] += scale_light*mi[kk];
                                    }
                                    
                                }
                                
                                extraFacets.erase(extraFacets.begin(), extraFacets.end());
                            }
                            
                            // update for next bounce
                            bounce++;
                            currF = test;
                            for (kk=0; kk<3; kk++) {
                                lightSource[kk] = uhat_new[kk];
                            }
                            scale_light *= tempF.getRho() * tempF.getS();
                            
                        } else {
                            // Didn't hit anything so no more bounces
                            bounce = numBounces + 1;
                        }
                    }
                    
                    // if there was a non-zero force computed, store in class
                    if ((ftemp[0] != 0.0) || (ftemp[1] != 0.0) || (ftemp[2] != 0.0)) {
                        latSpec[ii].push_back(currDelS);
                        longSpec[ii].push_back(currLamS);
                        fSpec[ii].push_back(rowd);
                        tSpec[ii].push_back(rowd);
                        for (kk=0; kk<3; kk++) {
                            fSpec[ii][currSpecCase].push_back(ftemp[kk]);
                            tSpec[ii][currSpecCase].push_back(ttemp[kk]);
                        }
                        
                        currSpecCase++;
                    }
                }
            }
            
        }
        
    }
    
    return;
}

void SRPModel::compute_fi(const double* uhat, Facet hitF, double* fi)
{
    double rho = hitF.getRho();
    double s = hitF.getS();
    double a2 = hitF.geta2();
    double nn[3][3], Ar1[3], Ar2[3][3], Ar3[3][3];
    double Ar1_2[3], Ar2_2[3][3], Ar3_2[3][3];
    double uu[3][3], uun[3];
    double nhat[3];
    double* nhat_p;
    int ii, jj, kk;
    
    hitF.getOpticalFourier(nn, Ar1, Ar2, Ar3, Ar1_2, Ar2_2, Ar3_2);
    nhat_p = hitF.getNormal();
    nhat[0] = *nhat_p;
    nhat[1] = *(nhat_p+1);
    nhat[2] = *(nhat_p+2);
    
    for (ii=0; ii<3; ii++) {
        for (jj=0; jj<3; jj++) {
            uu[ii][jj] = uhat[ii]*uhat[jj];
        }
    }
    
    for (jj=0; jj<3; jj++) {
        uun[jj] = 0.0;
        for (kk=0; kk<3; kk++) {
            uun[jj] += uu[jj][kk]*nhat[kk];
        }
    }
    
    for (ii=0; ii<3; ii++) {
        fi[ii] =0.0;
        for (jj=0 ; jj<3; jj++) {
            fi[ii] += M_PI*rho*s*Ar3[ii][jj]*uun[jj] + M_PI*uu[ii][jj]*Ar1[jj] + M_PI*a2*Ar2[ii][jj]*uhat[jj];
        }
    }
    
}

void SRPModel::computeShadowBounds(const int dsindex, const int numf, const int lsnum, const double lambdaDel)
{
    vector<double> rowd;
    int kk, lowhigh;
    
    // Compute shadowing information
    for (int jj=0; jj<numf; jj++) {
        shadowBoundsLow[jj] = rowd;
        shadowBoundsHigh[jj] = rowd;
        // Know this facet is always shadowed if indicated in the rise/set information
        if (riseLambda[dsindex][jj] != 10.0) {
            // Step through shadowing results for each facet and find the lambda bounds
            lowhigh = 0;
            kk=0;
            while (kk<lsnum) {
                if ((lowhigh==0) && (shadowing[kk].getPcntLit(jj)==1)) {
                    // searching for a 1.0
                    shadowBoundsLow[jj].push_back(lambdaSunList[kk]);
                    lowhigh = 1;
                } else if ((lowhigh==1) && (shadowing[kk].getPcntLit(jj)==0)) {
                    shadowBoundsHigh[jj].push_back(lambdaSunList[kk-1]); // use previous index as it is the last one with getPcntLit == 1
                    lowhigh = 0;
                }
                kk++;
            }
            
            // If no end was found for the last segment, it should just be 2*PI
            if (shadowBoundsLow[jj].size() > shadowBoundsHigh[jj].size()) {
                shadowBoundsHigh[jj].push_back(lambdaSunList[lsnum-1]);
            }
        }
    }
    
    // Refine shadowing longitude boundaries if desired
    if (numShadRefine > 0) {
        int bound_count, shadFlag, mm, nn, testHit;
        double newLow, newHigh, midp, maxDim;
        double uhat[3], negu[3], pt[3];
        double* rc;
        double hit_pnt[3];
        vector<int> dummy;
        maxDim = Spacecraft.getMaxDimVoxGrid();
        for (int jj=0; jj<numf; jj++) {
            bound_count = shadowBoundsLow[jj].size();
            for (kk = 0; kk<bound_count; kk++) {
                if (shadowBoundsLow[jj][kk] > 0.0) {
                    shadFlag = 0;
                    newLow = shadowBoundsLow[jj][kk] - lambdaDel;
                    newHigh = shadowBoundsLow[jj][kk];
                    for (mm=0; mm<numShadRefine; mm++) {
                        midp = (newLow+newHigh)/2.0;
                        uhat[0] = cos(midp)*cos(deltaSunList[dsindex]);
                        uhat[1] = sin(midp)*cos(deltaSunList[dsindex]);
                        uhat[2] = sin(deltaSunList[dsindex]);
                        for (nn=0; nn<3; nn++) {
                            negu[nn] = -1.0*uhat[nn];
                        }
                        if ((Spacecraft.getInView(jj))->size()>0) {
                            rc = Spacecraft.getFacetRc(jj);
                            for (nn=0; nn<3; nn++) {
                                pt[nn] = rc[nn] + 2.0*maxDim*uhat[nn];
                            }
                            testHit = Spacecraft.ray_intersection(&pt[0], &negu[0], &hit_pnt[0], &dummy);
                            if (testHit != jj) {
                                // Facet was shadowed by another facet
                                if (shadFlag==0) {
                                    newLow = midp;
                                } else if (shadFlag==1) {
                                    newHigh = midp;
                                }
                            } else {
                                if (shadFlag==0) {
                                    newHigh = midp;
                                } else if (shadFlag==1) {
                                    newLow = midp;
                                }
                            }
                        }
                    }
                    if (shadFlag==0) {
                        shadowBoundsLow[jj][kk] = newHigh;
                    } else if (shadFlag==1) {
                        shadowBoundsHigh[jj][kk] = newLow;
                    }
                }
                if (shadowBoundsHigh[jj][kk] < 2.0*M_PI) {
                    shadFlag = 1;
                    newLow = shadowBoundsHigh[jj][kk];
                    newHigh = shadowBoundsHigh[jj][kk] + lambdaDel;
                    for (mm=0; mm<numShadRefine; mm++) {
                        midp = (newLow+newHigh)/2.0;
                        uhat[0] = cos(midp)*cos(deltaSunList[dsindex]);
                        uhat[1] = sin(midp)*cos(deltaSunList[dsindex]);
                        uhat[2] = sin(deltaSunList[dsindex]);
                        for (nn=0; nn<3; nn++) {
                            negu[nn] = -1.0*uhat[nn];
                        }
                        if ((Spacecraft.getInView(jj))->size()>0) {
                            rc = Spacecraft.getFacetRc(jj);
                            for (nn=0; nn<3; nn++) {
                                pt[nn] = rc[nn] + 2.0*maxDim*uhat[nn];
                            }
                            testHit = Spacecraft.ray_intersection(&pt[0], &negu[0], &hit_pnt[0], &dummy);
                            if (testHit != jj) {
                                // Facet was shadowed by another facet
                                if (shadFlag==0) {
                                    newLow = midp;
                                } else if (shadFlag==1) {
                                    newHigh = midp;
                                }
                            } else {
                                if (shadFlag==0) {
                                    newHigh = midp;
                                } else if (shadFlag==1) {
                                    newLow = midp;
                                }
                            }
                        }
                    }
                    if (shadFlag==0) {
                        shadowBoundsLow[jj][kk] = newHigh;
                    } else if (shadFlag==1) {
                        shadowBoundsHigh[jj][kk] = newLow;
                    }
                }
            }
        }
    }
}
}