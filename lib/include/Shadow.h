//
//  Shadow.h
//  Spacecraft_SRP
//
//  Created by Jay McMahon on 5/27/14.
//  Copyright (c) 2014 Jay McMahon. All rights reserved.
//

#ifndef __Spacecraft_SRP__Shadow__
#define __Spacecraft_SRP__Shadow__

//#include <iostream>
//#include <vector>
#include "Body.h"

//using namespace std;

class Shadow
{
public:
    Shadow();
    Shadow(int numfacets, double delta, double lambda);
    void ComputeShadowing(Body* Target, vector<double>* riseLam, vector<double>* setLam);
    void ComputeShadowingGrid(Body* Target, vector<double>* riseLam, vector<double>* setLam); // add grid defining inputs
    int getPcntLit(int fnum);
private:
    double deltaSun;
    double lambdaSun;
    vector <int> percentUnShadowed; // Need percent of area that sees the Sun for coeff calculations
    double dot ( double* a, double* b);
};

#endif /* defined(__Spacecraft_SRP__Shadow__) */
