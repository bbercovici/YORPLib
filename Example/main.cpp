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
//  main.cpp
//  Spacecraft_SRP
//
//  Created by Jay McMahon on 5/19/14.
//  Copyright (c) 2014 Jay McMahon. All rights reserved.
//

#include <sstream>

#include <SRPModel.h>

#include <chrono>
using namespace std;

int main(int argc, const char * argv[])
{

    auto start = std::chrono::system_clock::now();
    
    // Inputs have two choices, either a file was input on the command line, or there will be prompts on the screen
    // to start, just do the file option
    // When implment screen option, write the inputs to a file for re-use
    string objFileName, opticalFileName, outputFileBaseName;
    string inFileName = "../input_use1.txt";
    int numVox, opticalFlag;
    double rho, spec;
    double lambdaDel;
    double deltaDel;
    double MaxFourier;
    int howManyBounces;
    int numrefine;
    
    ifstream inFile(inFileName);
    string linestr;
    
    // First line is .obj file name
    getline(inFile, objFileName);
    
    // Second line is number of voxels per axis
    getline(inFile, linestr);
    numVox = stoi(linestr);
    
    // Third line contains either rho and s, or a file name to read from
    getline(inFile, linestr);
    
    // split line based on white space
    stringstream ss(linestr);
    string item;
    vector<string> elems;
    while (getline(ss, item, ' ')) {
        elems.push_back(item);
    }
    if (elems.size()==2) {
        // Two elements = rho then s
        rho = stod(elems[0]);
        spec = stod(elems[1]);
        opticalFlag = 0;
    } else {
        // One element = optical property file name
        opticalFileName = elems[0];
        opticalFlag = 1;
    }
    
    // Fourth line contains lambdaDel
    getline(inFile, linestr);

    lambdaDel = stod(linestr);
    
    // Fifth line contains deltaDel
    getline(inFile, linestr);
    deltaDel = stod(linestr);
    
    // Sixth line contains MaxFourier
    getline(inFile, linestr);
    MaxFourier = stoi(linestr);
    
    // Seventh line contains howManyBounces
    getline(inFile, linestr);
    howManyBounces = stoi(linestr);
    
    // Eigth line contains numrefine
    getline(inFile, linestr);
    numrefine = stoi(linestr);
    
    // Ninth line contains output file base name
    getline(inFile, outputFileBaseName);
    
    
    inFile.close();
    
    // Construct the Body object with the appropriate method based on opticalFlag
    YORPLib::Body targetObj;
    if (opticalFlag == 0) {
        YORPLib::Body inputObj(objFileName,rho,spec);
        targetObj = inputObj;
    } else if (opticalFlag == 1) {
        YORPLib::Body inputObj(objFileName,opticalFileName);
        targetObj = inputObj;

    }

    
    // Determine limits for voxel grid
    double xmax = 1.01 * targetObj.getMaxDim(1);
    double ymax = 1.01 * targetObj.getMaxDim(2);
    double zmax = 1.01 * targetObj.getMaxDim(3);
    
    // Fill out Voxel Grid
    targetObj.setVoxelGrid(xmax, ymax, zmax, numVox);
    
    // Compute and print SRP Fourier coefficients
    YORPLib::SRPModel targetSRP(lambdaDel, deltaDel, MaxFourier, &targetObj, howManyBounces, numrefine);
    targetSRP.writeSRPCoeffsFile(outputFileBaseName, 2*(90.0/deltaDel) + 1);
    
    // stop timer and report
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_secs = end - start;
    cout << "Time elapsed = " << elapsed_secs.count() << " s" << "\n";
    
    return 0;
}

