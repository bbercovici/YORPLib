//
//  Body.cpp
//  Spacecraft_SRP
//
//  Created by Jay McMahon on 5/19/14.
//  Copyright (c) 2014 Jay McMahon. All rights reserved.
//

#include "Body.h"
#include <algorithm>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

Body::Body(){
    // cerr << "Need inputs to define Body\n"; // change this to throw an error
    
    return;
}

Body::Body(string bodyFile, double rho, double spec)
{
    ifstream inFile ( bodyFile );
    numFacets = 0;
    numVerts = 0;
    
    if ( !inFile.is_open() ) {
        // The file could not be opened
        cerr << "Requested body file doesn't exist\n"; // change this to throw an error
    }
    else {


        // Safely use the file stream
        string linestr, item;
        vector<string> elems;
        
        // Assume that the input file is a .obj file, so the format has lines that start with "v" or "f"
        // get a line
        while (getline(inFile, linestr)) {
            // split line based on white space
            stringstream ss(linestr);
            while (getline(ss, item, ' ')) {
                elems.push_back(item);
            }

            if (elems[0]=="v") {
                vertices.push_back(stod(elems[1]));
                vertices.push_back(stod(elems[2]));
                vertices.push_back(stod(elems[3]));
                numVerts++;
            } else if (elems[0]=="f") {
                facetList.push_back(stoi(elems[1]));
                facetList.push_back(stoi(elems[2]));
                facetList.push_back(stoi(elems[3]));
                numFacets++;
            }
            elems.clear();
        }
    }
    inFile.close();

    // Call other member functions to fill out the rest of the object
    facetData.resize(numFacets);

    setFacets(rho, spec);
    neighbors.resize(3*numFacets);

    setNeighbors();
    inView.resize(numFacets);

    setView();
    f2vlist.resize(numVerts);
    setF2V();

}


Body::Body(std::vector<std::vector<double > > vertices, 
    std::vector<std::vector<int> > facets, double rho, double spec){
   
    numFacets = facets.size();
    numVerts = vertices.size();


    for (int i = 0; i < vertices.size(); ++i){
        this -> vertices.push_back(vertices[i][0]);
        this -> vertices.push_back(vertices[i][1]);
        this -> vertices.push_back(vertices[i][2]);
    }

    for (int i = 0; i < facets.size(); ++i){
        this -> facetList.push_back(facets[i][0] + 1);
        this -> facetList.push_back(facets[i][1] + 1);
        this -> facetList.push_back(facets[i][2] + 1);
    }
    

    // Call other member functions to fill out the rest of the object
    facetData.resize(numFacets);

    setFacets(rho, spec);
    neighbors.resize(3*numFacets);

    setNeighbors();
    inView.resize(numFacets);

    setView();
    f2vlist.resize(numVerts);
    setF2V();

}






Body::Body(string bodyFile)
{
    ifstream inFile ( bodyFile );
    numFacets = 0;
    numVerts = 0;
    
    if ( !inFile.is_open() ) {
        // The file could not be opened
        cerr << "Requested body file doesn't exist\n"; // change this to throw an error
    }
    else {
        // Safely use the file stream
        string linestr, item;
        vector<string> elems;
        
        // Assume that the input file is a .obj file, so the format has lines that start with "v" or "f"
        // get a line
        while (getline(inFile, linestr)) {
            // split line based on white space
            stringstream ss(linestr);
            while (getline(ss, item, ' ')) {
                elems.push_back(item);
            }
            if (elems[0]=="v") {
                vertices.push_back(stod(elems[1]));
                vertices.push_back(stod(elems[2]));
                vertices.push_back(stod(elems[3]));
                numVerts++;
            } else if (elems[0]=="f") {
                facetList.push_back(stoi(elems[1]));
                facetList.push_back(stoi(elems[2]));
                facetList.push_back(stoi(elems[3]));
                std::cout << facetList[numFacets - 3] << " " << facetList[numFacets - 2] << " " << facetList[numFacets - 1] << std::endl; 
                numFacets++;
            }
            elems.clear();
        }
    }
    inFile.close();
    
    // Call other member functions to fill out the rest of the object
    facetData.resize(numFacets);
    setFacets(0.0, 0.0);
    neighbors.resize(3*numFacets);
    setNeighbors();
    inView.resize(numFacets);
    setView();
    setF2V();
}

Body::Body(string bodyFile, string optFile)
{
    ifstream inFile ( bodyFile );
    ifstream inFile2 ( optFile );
    numFacets = 0;
    numVerts = 0;
    
    if ( !inFile.is_open() ) {
        // The file could not be opened
        cerr << "Requested body file doesn't exist\n"; // change this to throw an error
    }
    else {
        // Safely use the file stream
        string linestr, item;
        vector<string> elems;
        
        // Assume that the input file is a .obj file, so the format has lines that start with "v" or "f"
        // get a line
        while (getline(inFile, linestr)) {
            // split line based on white space
            stringstream ss(linestr);
            while (getline(ss, item, ' ')) {
                elems.push_back(item);
            }
            if (elems[0]=="v") {
                vertices.push_back(stod(elems[1]));
                vertices.push_back(stod(elems[2]));
                vertices.push_back(stod(elems[3]));
                numVerts++;
            } else if (elems[0]=="f") {
                facetList.push_back(stoi(elems[1]));
                facetList.push_back(stoi(elems[2]));
                facetList.push_back(stoi(elems[3]));
                numFacets++;
            }
            elems.clear();
        }
    }
    inFile.close();
    
    vector<double> rhoVec(numFacets);
    vector<double> specVec(numFacets);
    
    if ( !inFile2.is_open() ) {
        // The file could not be opened
        cerr << "Requested optical parameter file doesn't exist\n"; // change this to throw an error
    }
    else {
        // Safely use the file stream
        string linestr, item;
        vector<string> elems;
        int counter = 0;
        
        // This file has one line for every facet, each line has two double values: reflectivity (rho) and the percent of that which is specular (s)
        while (getline(inFile2, linestr)) {
            // split line based on white space
            stringstream ss(linestr);
            while (getline(ss, item, ' ')) {
                elems.push_back(item);
            }
            rhoVec[counter] = stod(elems[0]);
            specVec[counter] = stod(elems[1]);
            counter++;
            elems.clear();
        }
    }
    inFile2.close();
    
    // Call other member functions to fill out the rest of the object
    facetData.resize(numFacets);
    setFacetsVec(&rhoVec, &specVec);
    neighbors.resize(3*numFacets);
    setNeighbors();
    inView.resize(numFacets);
    setView();
    f2vlist.resize(numVerts);
    setF2V();
}

void Body::setNeighbors()
{
    int i, j;
    int x, y, z, a, b, c;
    int xc, yc, zc;
    vector<int> numneigh(numFacets,0);
    
    for (i=0;i<numFacets-1;i++){
        x = facetList[3*i];
        y = facetList[3*i+1];
        z = facetList[3*i+2];
        for (j=i+1;j<numFacets;j++){
            if (numneigh[i]==3) break;
            a = facetList[3*j];
            b = facetList[3*j+1];
            c = facetList[3*j+2];
            xc = (x==a || x==b || x==c);
            yc = (y==a || y==b || y==c);
            zc = (z==a || z==b || z==c);
            if ((xc&&yc)||(xc&&zc)||(yc&&zc)) {
                neighbors[3*i+numneigh[i]] = j;
                numneigh[i]++;
                neighbors[3*j+numneigh[j]] = i;
                numneigh[j]++;
            }
        }
    }
    
    return;
}

void Body::setFacets(double rho, double s)
{
    int ii;
    int verts[3];

    
    
    for (ii=0; ii<numFacets; ii++) {
        verts[0] = facetList[3*ii];
        verts[1] = facetList[3*ii + 1];
        verts[2] = facetList[3*ii + 2];

        Facet tempFacet(&verts[0], &vertices, rho, s);

        tempFacet.computeOptical();

        facetData[ii] = tempFacet;
    }
}

void Body::setFacetsVec(vector<double>* rho, vector<double>* s)
{
    int ii;
    int verts[3];
    
    for (ii=0; ii<numFacets; ii++) {
        verts[0] = facetList[3*ii];
        verts[1] = facetList[3*ii + 1];
        verts[2] = facetList[3*ii + 2];
        Facet tempFacet(&verts[0], &vertices, (*rho)[ii], (*s)[ii]);
        tempFacet.computeOptical();
        facetData[ii] = tempFacet;
    }
}

void Body::setView()
{
    int ii, jj, kk, mm, nn;
    double rtemp[3], rmag;
    double* rc1;
    double* rc2;
    bool flag;
    
    for (ii=0; ii<numFacets; ii++){
        vector<int>* rowP = new vector<int>;
        //row.resize((numFacets-1)/2);
        //rowP->resize(2000); // 4530 per is max?
        inView[ii] = rowP;
    }
    
    // Check for facets that COULD see each other
    // Test relative to all vertices to account for all alignments; break as soon as a true statement is found
//    int vCount;
    for (ii=0; ii<numFacets; ii++){
//        rc1 = facetData[ii].getRc();

/*        vector<int>* rowP = new vector<int>;
        rowP->resize((numFacets-ii+5)/5);
        inView[ii] = rowP; */
//        vCount = 0;

//        for (jj=ii+1; jj<numFacets; jj++){
        for (jj=0; jj<numFacets; jj++){
//            rc2 = facetData[jj].getRc();
            if (ii==jj) {
                continue;
            }
            flag = false;
            for (mm=0; mm<3; mm++) {
                rc1 = &vertices[facetData[ii].getVert(mm)];
                for (nn=0; nn<3; nn++) {
                    rc2 = &vertices[facetData[jj].getVert(nn)];
                    for (kk=0; kk<3; kk++) {
                        rtemp[kk] = rc2[kk] - rc1[kk];
                    }
                    rmag = sqrt(rtemp[0]*rtemp[0] + rtemp[1]*rtemp[1] + rtemp[2]*rtemp[2]);
                    for (kk=0; kk<3; kk++) {
                        rtemp[kk] /= rmag;
                    }
                    if ((dot(&rtemp[0], facetData[ii].getNormal()) > 1.0e-10) && (dot(&rtemp[0], facetData[jj].getNormal()) < -1.0e-10)) { // JWM 2/10/16 took out >=/<= to get rid of planar facets seeing eachother
                        // JWM 2/24/16 normalized rtemp; set limit so vertex has to be > 1e-10 radians above horizon to count as in view
//                        try {
                            //if (vCount >= inView[ii]->size()) {
                            //    inView[ii]->resize(inView[ii]->size() + 2000);
                                //cout << "size inVew[" << ii << "] = " << inView[ii]->size() << "\n";
                            //}
                            //inView[ii]->at(vCount) = jj;
                        inView[ii]->push_back(jj);
//                            vCount++;
                            //inView[ii]->push_back(jj);
                            //inView[jj]->push_back(ii);
//                        } catch (...) {
//                            cout << "size inVew[" << ii << "] = " << inView[ii]->size() << "\n";
//                            cout << "size inVew[" << jj << "] = " << inView[jj]->size() << "\n";
//                            cout << "jay\n";
//                        }
                        flag = true;
                        break;
                    }
                }
                if (flag==true) {
                    break;
                }
            }
            
//            if (vCount > 0) {
//                break;
//            }
        }
        
        
        
        //inView[ii]->shrink_to_fit();
        
        /*if (ii==77){
            for (mm=0; mm<inView[ii]->size(); mm++) {
                cout << inView[ii]->at(mm) << "\n";
            }
            cout << "jay\n";
        }*/
        
    }
    
    
//    for (ii=0; ii<numFacets; ii++){
//        inView[ii]->shrink_to_fit();
//    }
    
    viewCheck = true;
    
}

void Body::setViewRc()
{
    // Check for intervening facets based only on center-to-center view
    int ii, jj, kk, test;
    double hit_pnt[3];
    double uhat[3], vmag, pt[3];
    double* rc1;
    double* rc2;
    vector<int> dummy;
    
    // If for some reason the view hasn't been set yet, do so now
    if (!viewCheck) {
        setView();
    }
    
    // Initialize inViewRc as a copy of inView
    inViewRc.resize(numFacets);
    for (ii=0; ii<numFacets; ii++){
        vector<int>* rowP = new vector<int>;
        inViewRc[ii] = rowP;
        for (jj=0; jj<(*inView[ii]).size(); jj++){
            inViewRc[ii]->push_back((*inView[ii])[jj]);
        }
    }
    
    // Check every center-to-center ray trace to see if it intersects the facet that is "inView"
    for (ii=0; ii<numFacets; ii++) {
        if (!inViewRc[ii]->empty()) {
            rc1 = facetData[ii].getRc();
            for (jj=0; jj<inViewRc[ii]->size(); jj++) {
                // skip if this connection already checked because these are symmetric
                if (ii > inViewRc[ii]->at(jj)) {
                    continue;
                }
                
                // set the test direction from facet ii to the facet in inViewRc[ii][jj]
                rc2 = facetData[inViewRc[ii]->at(jj)].getRc();
                for (kk=0; kk<3; kk++) {
                    uhat[kk] = rc2[kk] - rc1[kk];
                }
                vmag = sqrt(dot(&uhat[0], &uhat[0])); // JWM 2/10/16 moved normalization up here so the horizon checking would actually work!
                for (kk=0; kk<3; kk++) {
                    uhat[kk] = uhat[kk]/vmag;
                }
                
                // ray tracing doesn't detect if it goes below either facets' horizon - use previous test for this
                // JWM 2/24/16 modified to have to be > 1e-10 radians above horizon to count as in view
                if ((dot(&uhat[0], facetData[ii].getNormal()) < -1.0e-10) || (dot(&uhat[0], facetData[inViewRc[ii]->at(jj)].getNormal()) > 1.0e-10)) {
                    // This means the center connecting line goes through the facets, so erase
                    // Erase the mirror inViewRc entry also
                    inViewRc[inViewRc[ii]->at(jj)]->erase(find(inViewRc[inViewRc[ii]->at(jj)]->begin(), inViewRc[inViewRc[ii]->at(jj)]->end(), ii));
                    inViewRc[ii]->erase(inViewRc[ii]->begin()+jj);
                    jj--;
                } else {
                    for (kk=0; kk<3; kk++) {
                        pt[kk] = rc1[kk] + 0.001*vmag*uhat[kk]; // move starting point 0.1% of the distance along the connecting line so the ray tracing check doesn't just return the starting facet
                    }
                    test = bodyVox.ray_intersect_voxel(&pt[0], &uhat[0], &facetList, &vertices, &facetData, &hit_pnt[0], &dummy);
                    if (inViewRc[ii]->at(jj) !=  test) {
                        // This means we intersected a different facet before the targeted one, so need to erase
                        // Erase the mirror inViewRc entry also
                        // Note that we don't check for multiple facets because if it hits the facet it is supposed to, it won't hit multiple facets
                        inViewRc[inViewRc[ii]->at(jj)]->erase(find(inViewRc[inViewRc[ii]->at(jj)]->begin(), inViewRc[inViewRc[ii]->at(jj)]->end(), ii));
                        inViewRc[ii]->erase(inViewRc[ii]->begin()+jj);
                        jj--;
                    }
                }
            }
        }
    }
}

void Body::setViewGrid()
{
    // Using a grid requires extra outputs beyond just inViewGrid, which says what facets can see each other at all.
    // Also need to output how much they can see of one another, and the center of pressure offsets from the original center of facet.

    // Check for intervening facets based on vertex-to-vertex views, and compute visible portions
    int ii, jj, kk, test, mm, pp, vCount;
    double hit_pnt[3];
    double uhat[3], vmag, pt[3];
    double* rc1;
    double* rc2;
    vector<int> dummy;
    
    // If for some reason the view hasn't been set yet, do so now
    if (!viewCheck) {
        setView();
    }
    
    // If grid hasn't been initialized do so now
    
    inViewGrid = inView;
    
    // Check every vertex-to-vertex ray trace to see if it intersects the facet that is "inView"
    for (ii=0; ii<numFacets; ii++) {
        if (!inViewGrid[ii]->empty()) {

            for (jj=0; jj<inViewGrid[ii]->size(); jj++) {
                // skip if this connection already checked because these are symmetric
                if (ii > inViewGrid[ii]->at(jj)) {
                    continue;
                }
                vCount = 0;
                for (mm=0; mm<3; mm++) {
                    rc1 = &vertices[facetData[ii].getVert(mm)];
                    for (pp=0; pp<3; pp++) {
                        // set the test direction from facet ii to the facet in inViewGrid[ii][jj]
                        rc2 = &vertices[facetData[inViewGrid[ii]->at(jj)].getVert(pp)];
                        for (kk=0; kk<3; kk++) {
                            uhat[kk] = rc2[kk] - rc1[kk];
                        }
                        vmag = sqrt(dot(&uhat[0], &uhat[0]));
                        for (kk=0; kk<3; kk++) {
                            uhat[kk] = uhat[kk]/vmag;
                            pt[kk] = rc1[kk] + 0.001*vmag*uhat[kk]; // move starting point 0.1% of the distance along the connecting line so the ray tracing check doesn't just return the starting facet
                        }
                        test = bodyVox.ray_intersect_voxel(&pt[0], &uhat[0], &facetList, &vertices, &facetData, &hit_pnt[0], &dummy);
                        if (test == inViewGrid[ii]->at(jj)) {
                            ++vCount;
                        }
                    }
                }
                
                
                if (vCount == 0) {
                    // This means we intersected a different facet before the targeted one, so need to erase
                    // Erase the mirror inViewGrid entry also
                    inViewGrid[inViewGrid[ii]->at(jj)]->erase(find(inViewGrid[inViewGrid[ii]->at(jj)]->begin(), inViewGrid[inViewGrid[ii]->at(jj)]->end(), ii));
                    inViewGrid[ii]->erase(inViewGrid[ii]->begin()+jj);
                } else if (vCount < 9) {
                    // This means the view is partially blocked - compute grid outputs

                }
            }
        }
    }
}


void Body::setF2V()
{
    vector<int>::iterator iter_f;
    vector < vector<int> >::iterator iter_f2v;
    vector<int> row;
    int ii = 1; // starts at 1 because facetlist has vertex numbers starting at 1
    int jj;
    
    for (iter_f2v = f2vlist.begin(); iter_f2v != f2vlist.end(); ++iter_f2v) {
        *iter_f2v = row;
        iter_f = facetList.begin();
        while (iter_f != facetList.end()) {
            iter_f = find(iter_f, facetList.end(), ii);
            if (iter_f != facetList.end()) {
                jj = iter_f - facetList.begin();
                iter_f2v->push_back(jj/3); // will give the vertex number, starting at 0
                ++iter_f;
            }
        }
        ii++;
    }
    
}


void Body::setVoxelGrid(double xmax_in, double ymax_in, double zmax_in, int N)
{
    VoxelGrid bodyVoxTemp(xmax_in, ymax_in, zmax_in, N);
    bodyVox = bodyVoxTemp;
    bodyVox.fillGrid(&facetList, &vertices, numVerts, &f2vlist, &facetData);
}

int Body::ray_intersection(double* pt, double* uhat, double* hit_pnt, vector<int>* extraFacets)
{
    int test_hit;
    test_hit = bodyVox.ray_intersect_voxel(pt, uhat, &facetList, &vertices, &facetData, hit_pnt, extraFacets);
    return test_hit;
}

double Body::dot ( double* a, double* b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    
}

int Body::getNumFacets()
{
    return numFacets;
}

int Body::getNumVerts()
{
    return numVerts;
}

/*void Body::getFacetInfo(int fnum, Facet* fac)
{
    cout << &facetData.at(fnum);
    fac = &facetData.at(fnum);
    return;
}*/

Facet Body::getFacet( int fnum)
{
    return facetData[fnum];
}

double* Body::getFacetNormal(int fnum)
{
    return facetData[fnum].getNormal();
}

double* Body::getFacetRc(int fnum)
{
    return facetData[fnum].getRc();
}

double Body::getFacetArea(int fnum)
{
    return facetData[fnum].getArea();
}

int Body::getVinVox(int voxnum)
{
    return bodyVox.getNumVinVox(voxnum);
}

int Body::getFinVox(int voxnum)
{
    return bodyVox.getNumFinVox(voxnum);
}

int* Body::getNeighbors(int fnum)
{
    return &neighbors[3*fnum];
}

vector<int>* Body::getInView(int fnum)
{
    return inView[fnum];
}

vector<int>* Body::getInViewRc(int fnum)
{
    return inViewRc[fnum];
}

vector<int>* Body::getF2V(int fnum)
{
    return &f2vlist[fnum];
}

double Body::getMaxDimVoxGrid()
{
    return bodyVox.getMaxDim();
}

void Body::getOpticalFourier(int fnum, double& rho, double& s, double& a2, double nn[3][3], double Ar1[3], double Ar2[3][3], double Ar3[3][3], double Ar1_2[3], double Ar2_2[3][3], double Ar3_2[3][3])
{
//    double temp;
    a2 = facetData[fnum].geta2();
//    a2 = &temp;
//    temp = facetData[fnum].getRho();
//    rho = &temp;
    rho = facetData[fnum].getRho();
    s = facetData[fnum].getS();
//    s = &temp;
    
    facetData[fnum].getOpticalFourier(nn, Ar1, Ar2, Ar3, Ar1_2, Ar2_2, Ar3_2);
    
    return;
    
}

double Body::getMaxDim(int axis)
{
    // axis == 1 is x axis
    // axis == 2 is y axis
    // axis == 3 is z axis

    int offset;
    
    if (axis == 1) {
        offset = 0;
    } else if (axis == 2) {
        offset = 1;
    } else if (axis == 3) {
        offset = 2;
    }
    
    vector<double> data;
    for (int ii = offset; ii<3*numVerts; ii+=3) {
        data.push_back(abs(vertices[ii]));
    }
    
    return *max_element(data.begin(), data.end());
    
}