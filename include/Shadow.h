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
