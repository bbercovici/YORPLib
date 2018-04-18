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
//  VoxelGrid.cpp
//  Spacecraft_SRP
//
//  Created by Jay McMahon on 5/19/14.
//  Copyright (c) 2014 Jay McMahon. All rights reserved.
//

#include "YORPLib/VoxelGrid.h"
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;
namespace YORPLib{

    VoxelGrid::VoxelGrid()
    {
        VoxelGrid(1.0,1.0,1.0,1);
    }

    VoxelGrid::VoxelGrid(double xmax_in, double ymax_in, double zmax_in, int N_in){
        N = N_in;
        dx = 2.0*xmax_in/(double)N;
        dy = 2.0*ymax_in/(double)N;
        dz = 2.0*zmax_in/(double)N;
        
    // Resize vectors in class
        int N3 = N*N*N;
        int N2 = N*N;
        xmax.resize(N3);
        xmin.resize(N3);
        ymax.resize(N3);
        ymin.resize(N3);
        zmax.resize(N3);
        zmin.resize(N3);
        posx.resize(N3);
        negx.resize(N3);
        posy.resize(N3);
        negy.resize(N3);
        posz.resize(N3);
        negz.resize(N3);
        
        flist.resize(N3);
        vlist.resize(N3);
        fnum.resize(N3);
        vnum.resize(N3);
        
        vector<int> row;
        
        
        int vox_ind = 0;
        
        int ii, jj, kk;
        
        for (ii = 0; ii<N; ii++)
        {
            for (jj = 0; jj<N; jj++)
            {
                for (kk = 0; kk<N; kk++)
                {
                    xmin[vox_ind] = -xmax_in + ((double)kk)*dx;
                    xmax[vox_ind] = -xmax_in + ((double)kk+1.0)*dx;
                    ymin[vox_ind] = -ymax_in + ((double)jj)*dy;
                    ymax[vox_ind] = -ymax_in + ((double)jj+1.0)*dy;
                    zmin[vox_ind] = -zmax_in + ((double)ii)*dz;
                    zmax[vox_ind] = -zmax_in + ((double)ii+1.0)*dz;
                    
                // Define voxel X neighbors
                    if (kk==0) {
                        posx[vox_ind] = vox_ind + 1;
                    negx[vox_ind] = N3; // N3 as an index will be out of bounds - indicates no neighbor in that direction
                } else if (kk==N-1) {
                    posx[vox_ind] = N3;
                    negx[vox_ind] = vox_ind - 1;
                } else {
                    posx[vox_ind] = vox_ind + 1;
                    negx[vox_ind] = vox_ind - 1;
                }
                
                // Define voxel Y neighbors
                if (jj==0) {
                    posy[vox_ind] = vox_ind + N;
                    negy[vox_ind] = N3;
                } else if (jj==N-1) {
                    posy[vox_ind] = N3;
                    negy[vox_ind] = vox_ind - N;
                } else {
                    posy[vox_ind] = vox_ind + N;
                    negy[vox_ind] = vox_ind - N;
                }
                
                // Define voxel Z neighbors
                if (ii==0) {
                    posz[vox_ind] = vox_ind + N2;
                    negz[vox_ind] = N3;
                } else if (ii==N-1) {
                    posz[vox_ind] = N3;
                    negz[vox_ind] = vox_ind - N2;
                } else {
                    posz[vox_ind] = vox_ind + N2;
                    negz[vox_ind] = vox_ind - N2;
                }
                
                // Initialize facet and vertex lists/numbers
                flist[vox_ind] = row;
                vlist[vox_ind] = row;
                fnum[vox_ind] = 0;
                vnum[vox_ind] = 0;
                
                // Increment index
                vox_ind++;
            }
        }
    }
}

void VoxelGrid::fillGrid(vector<int>* facets, vector<double>* vertices, int numVerts, vector< vector<int> >* f2vlist, vector<Facet>* fdata)
{
    int ii;
    vector<int> voxel_list;
    vector<int>::iterator iter_vl;
    vector< vector<int> > v2voxlist; // local list of voxels each vertex is in
    v2voxlist.resize(numVerts);
    
    //num_vertices = size(vertices,1);
    
    //f_2_v_list = facets_per_vertex(facets,num_vertices);
    ii = 0;
    for (vector<double>::iterator iter = vertices->begin(); iter != vertices->end(); iter+=3) {
        
        double v1[3];
        v1[0] = *iter;
        v1[1] = *(iter+1);
        v1[2] = *(iter+2);
        
        find_pt_in_voxel(&v1[0], &voxel_list);
        v2voxlist[ii] = voxel_list;
        
        for (iter_vl = voxel_list.begin(); iter_vl != voxel_list.end(); ++iter_vl) {
            
            // add vertex and associated facets to this voxel's list
            vlist[*iter_vl].push_back(ii);
            vnum[*iter_vl]++;
//            flist[*iter_vl].reserve(flist[*iter_vl].size() + f2vlist[ii].size());
            copy((*f2vlist)[ii].begin(), (*f2vlist)[ii].end(), back_inserter(flist[*iter_vl]));
            sort(flist[*iter_vl].begin(),flist[*iter_vl].end());
            flist[*iter_vl].erase(unique(flist[*iter_vl].begin(),flist[*iter_vl].end()),flist[*iter_vl].end());
            fnum[*iter_vl] = flist[*iter_vl].size();
        }
        voxel_list.erase(voxel_list.begin(),voxel_list.end());
        ii++;
    }
    
    // Check each facet's vertices to see if the facet spans more than one voxel
    vector <int> facet_voxs, checkList;
    int temp = 0, temp2 = 0, temp3 = 0, jj, kk, mm;
    int N2 = N*N;
    int N3 = N2*N;
    int fxmin,fxmax,fymin,fymax,fzmin,fzmax;
    ii=0; // index of facets (for flist/fnum)
    for (vector<int>::iterator iterf = facets->begin(); iterf != facets->end(); iterf+=3) {
        // Recall that facetList contains the numbers from the input file, so the value referenced by iterf needs to have 1 subtracted to start at index 0
        facet_voxs = v2voxlist[*(iterf)-1];
        copy(v2voxlist[*(iterf+1)-1].begin(), v2voxlist[*(iterf+1)-1].end(), back_inserter(facet_voxs));
        copy(v2voxlist[*(iterf+2)-1].begin(), v2voxlist[*(iterf+2)-1].end(), back_inserter(facet_voxs));
        sort(facet_voxs.begin(), facet_voxs.end());
        facet_voxs.erase(unique(facet_voxs.begin(), facet_voxs.end()), facet_voxs.end());
        if (facet_voxs.size()>1) {
            if (facet_voxs.size()==2) {
                // if there are only two voxels and they share a face that is all possible
                temp = abs(facet_voxs[1] - facet_voxs[0]);
                if ((temp==1) || (temp==N) || (temp == N2)) {
                    ii++;
                    continue;
                }
            }
            // otherwise, check all voxels in between the max/min on each axis to see if the triangle intersects
            fxmin = N3;
            fxmax = -1;
            fymin = N3;
            fymax = -1;
            fzmin = N3;
            fzmax = -1;
            for (jj=0; jj<facet_voxs.size(); jj++) {
                getVoxInd(facet_voxs[jj],&temp,&temp2,&temp3);
                fxmin = min(temp,fxmin);
                fxmax = max(temp,fxmax);
                fymin = min(temp2,fymin);
                fymax = max(temp2,fymax);
                fzmin = min(temp3,fzmin);
                fzmax = max(temp3,fzmax);
            }
            // create list of voxels to check
            checkList.erase(checkList.begin(),checkList.end());
            for (jj=fzmin; jj<=fzmax; jj++) {
                for (kk=fymin; kk<=fymax; kk++) {
                    for (mm=fxmin; mm<=fxmax; mm++) {
                        temp = mm + kk*N + jj*N2;
                        if (find(facet_voxs.begin(),facet_voxs.end(),temp)==facet_voxs.end()) {
                            checkList.push_back(temp);
                        }
                    }
                }
            }
            // Check to see if the current facet intersects each voxel in checkList
            for (vector<int>::iterator iterc = checkList.begin(); iterc != checkList.end(); ++iterc) {
                // function returns true if an intersection
                if (intersectTriVox(*iterc, &(*vertices)[(*fdata)[ii].getVert(0)], &(*vertices)[(*fdata)[ii].getVert(1)], &(*vertices)[(*fdata)[ii].getVert(2)], (*fdata)[ii].getNormal())) {
                    // there is an intersection; add facet to appropriate voxel
                    flist[*iterc].push_back(ii);
                    fnum[*iterc]++;
                }
            }
        }
        //facet_voxs.erase(facet_voxs.begin(),facet_voxs.end());
        ii++;
    }
    
    // Remove duplicate facet entries
    int vox_ind = 0;
    for (ii = 0; ii<N; ii++) {
        for (jj = 0; jj<N; jj++) {
            for (kk = 0; kk<N; kk++) {
                sort(flist[vox_ind].begin(),flist[vox_ind].end());
                flist[vox_ind].erase(unique(flist[vox_ind].begin(), flist[vox_ind].end()), flist[vox_ind].end());
                fnum[vox_ind] = flist[vox_ind].size();
                vox_ind++;
            }
        }
    }
    
    return;
}

void VoxelGrid::getVoxInd(int voxnum, int* x, int* y, int* z)
{
    // Get the x, y, and z indices of a given voxel
    *x = voxnum % N;
    *y = ((voxnum%(N*N))-*x)/N;
    *z = (voxnum - *x - *y*N)/(N*N);
    return;
}

void VoxelGrid::find_pt_in_voxel(double* pt, vector<int>* voxel_list)
{
    // Test to ensure the pt is inside the total voxel_grid box
    if ((pt[0] < xmin[0]) || (pt[0] > abs(xmin[0])) || (pt[1] < ymin[0]) || (pt[1] > abs(ymin[0])) || (pt[2] < zmin[0]) || (pt[2] > abs(zmin[0]))) {
        return;
    }
    
    // Find voxels this pt is inside of (or on the border)
    // I know that the first voxel is at the negative end of all three axes
    // Account for numerical rounding by comparing to floor
    double num_accuracy = 1e-14;
    double xgrid, ygrid, zgrid;
    int xind, yind, zind;
    
    xgrid = (pt[0] - xmin[0])/dx;
    if (abs(xgrid - round(xgrid)) < num_accuracy) {
        xgrid = round(xgrid);
    }
    xind = (int)floor(xgrid);
    if (xind == N) {
        xind--;
    }
    ygrid = (pt[1] - ymin[0])/dy;
    if (abs(ygrid - round(ygrid)) < num_accuracy) {
        ygrid = round(ygrid);
    }
    yind = (int)floor(ygrid);
    if (yind == N) {
        yind--;
    }
    zgrid = (pt[2] - zmin[0])/dz;
    if (abs(zgrid - round(zgrid)) < num_accuracy) {
        zgrid = round(zgrid);
    }
    zind = (int)floor(zgrid);
    if (zind == N) {
        zind--;
    }
    
    voxel_list->push_back(xind + yind*N + zind*N*N);
    
    // check neighboring voxels
    if ((xgrid == xind) && (xind != 0)) {
        voxel_list->push_back((*voxel_list)[0] - 1);
        if ((ygrid == yind) && (yind != 0)) {
            voxel_list->push_back((*voxel_list)[0] - N - 1);
            if ((zgrid == zind) && (zind != 0)) {
                voxel_list->push_back((*voxel_list)[0] - N*N - N - 1);
            }
        }
        if ((zgrid == zind) && (zind != 0)) {
            voxel_list->push_back((*voxel_list)[0] - N*N - 1);
        }
    }
    if ((ygrid == yind) && (yind != 0)) {
        voxel_list->push_back((*voxel_list)[0] - N);
        if ((zgrid == zind) && (zind != 0)) {
            voxel_list->push_back((*voxel_list)[0] - N*N - N);
        }
    }
    if ((zgrid == zind) && (zind != 0)) {
        voxel_list->push_back((*voxel_list)[0] - N*N);
    }
    
    return;
}

int VoxelGrid::ray_intersect_voxel(double* pt, double* uhat, vector<int>* facets, vector<double>* vertices, vector <Facet>* facetData, double* hit_pnt, vector<int>* extraFacets)
{
    int hit_facet = -1;
    vector<int> voxel_list;
    double tX, tY, tZ, t;
    int current_voxel;
    
    // Find which voxel ray starts in
    find_pt_in_voxel(pt, &voxel_list);
    
    // If starts outside voxel_grid, find where the ray intersects
    while (voxel_list.size() == 0) {
        // Figure out where first intersection occurs
        t = 0.0;
        if (abs(uhat[0]) < 1e-16) {
            tX = -1.0;
        } else {
            if (pt[0] > 0) {
                tX = (abs(xmin[0]) - pt[0])/uhat[0];
            } else {
                tX = (xmin[0] - pt[0])/uhat[0];
            }
        }
        if (tX > 0.0) {
            t = tX;
        }
        if (abs(uhat[1]) < 1e-16) {
            tY = -1.0;
        } else {
            if (pt[1] > 0.0) {
                tY = (abs(ymin[0]) - pt[1])/uhat[1];
            } else {
                tY = (ymin[0] - pt[1])/uhat[1];
            }
        }
        if ((tY > 0.0) && ((t==0.0) || (tY < t))) {
            t = tY;
        }
        if (abs(uhat[2]) < 1e-16) {
            tZ = -1.0;
        } else {
            if (pt[2] > 0.0) {
                tZ = (abs(zmin[0]) - pt[2])/uhat[2];
            } else {
                tZ = (zmin[0] - pt[2])/uhat[2];
            }
        }
        if ((tZ > 0.0) && ((t==0.0) || (tZ < t))) {
            t = tZ;
        }
        
        // If ray is pointed away from the voxel_grid, return no intersection
        if (t == 0.0) {
            return hit_facet; // returns -1 if no intersection
        }
        
        // Update pt
        pt[0] = pt[0] + t*uhat[0];
        pt[1] = pt[1] + t*uhat[1];
        pt[2] = pt[2] + t*uhat[2];
        
        // Check to see if this intersect a voxel now
        find_pt_in_voxel(pt, &voxel_list);
        
    }
    
    if (voxel_list.size() > 1) {
        // the ray starts on a border - figure out which of the list is the
        // voxel the ray is going to travel through
        // It is possible that the ray will stay on the border of two voxels -
        // in that case, just choose one as the equalities for testing facets
        // will work it out appropriately.
        double stepsize = min(min(dx, dy), dz)/10.0;
        double pt2[3];
        vector<int> voxel_list2;
        pt2[0] = pt[0] + stepsize*uhat[0];
        pt2[1] = pt[1] + stepsize*uhat[1];
        pt2[2] = pt[2] + stepsize*uhat[2];
        find_pt_in_voxel(pt2, &voxel_list2);
        
        if ((abs(uhat[0])< 1e-13) || (abs(uhat[1])< 1e-13) || (abs(uhat[2])< 1e-13)) {
            // ray will travel on some border - decide appropriate facet to test
            while (find(voxel_list.begin(),voxel_list.end(),voxel_list2[0])==voxel_list.end()) {
                stepsize = stepsize/10.0;
                pt2[0] = pt[0] + stepsize*uhat[0];
                pt2[1] = pt[1] + stepsize*uhat[1];
                pt2[2] = pt[2] + stepsize*uhat[2];
                voxel_list2.erase(voxel_list2.begin(),voxel_list2.end());
                find_pt_in_voxel(pt2, &voxel_list2);
            }
            current_voxel = voxel_list2[0];

        } else {
            // ray will enter one of the current voxels uniquely - find it
            // ensure that we find one unique voxel, and that it is one of the originally bordered voxels
            while ((voxel_list2.size() > 1) || (find(voxel_list.begin(),voxel_list.end(),voxel_list2[0])==voxel_list.end())) {
                stepsize = stepsize/10.0;
                pt2[0] = pt[0] + stepsize*uhat[0];
                pt2[1] = pt[1] + stepsize*uhat[1];
                pt2[2] = pt[2] + stepsize*uhat[2];
                voxel_list2.erase(voxel_list2.begin(),voxel_list2.end());
                find_pt_in_voxel(pt2, &voxel_list2);
            }
            current_voxel = voxel_list2[0];
        }
        
    } else {
        current_voxel = voxel_list[0];
    }
    
    // Get step sizes to cross current voxel and full voxels in each direction
    double tDeltaX, tDeltaY, tDeltaZ;
    double tMaxX, tMaxY, tMaxZ;
    
    if (abs(uhat[0]) < 1e-16) {
        tDeltaX = 2.0*dx*N;
        tMaxX = 2.0*dx*N;
    } else {
        tDeltaX = dx/abs(uhat[0]);
        if (uhat[0] > 0.0) {
            tMaxX = (xmax[current_voxel] - pt[0])/uhat[0];
        } else {
            tMaxX = (xmin[current_voxel] - pt[0])/uhat[0];
        }
    }
    if (abs(uhat[1]) < 1e-16) {
        tDeltaY = 2.0*dy*N;
        tMaxY = 2.0*dy*N;
    } else {
        tDeltaY = dy/abs(uhat[1]);
        if (uhat[1] > 0.0) {
            tMaxY = (ymax[current_voxel] - pt[1])/uhat[1];
        } else {
            tMaxY = (ymin[current_voxel] - pt[1])/uhat[1];
        }
    }
    if (abs(uhat[2]) < 1e-16) {
        tDeltaZ = 2.0*dz*N;
        tMaxZ = 2.0*dz*N;
    } else {
        tDeltaZ = dz/abs(uhat[2]);
        if (uhat[2] > 0.0) {
            tMaxZ = (zmax[current_voxel] - pt[2])/uhat[2];
        } else {
            tMaxZ = (zmin[current_voxel] - pt[2])/uhat[2];
        }
    }
    
    // Traverse through the voxel grid
    int N3 = N*N*N;
    while (current_voxel != N3) {
        
        // Conduct test - if the test is satisfied then return
        // Note that we have to make sure the test function returns an empty
        // hit_facet if there is no hit
        // If current voxel doesn't contain any facets, just traverse to next
        // voxel. No need to take the time to call the function
        if (fnum[current_voxel] != 0) {
            hit_facet = test_intersection(pt, uhat, facets, vertices, facetData, hit_pnt, current_voxel, extraFacets);
            
            if (hit_facet != -1) {
                return hit_facet;
            }
        }
        
        // Traverse to next voxel
        if (tMaxX < tMaxY) {
            if (tMaxX < tMaxZ) {
                if (uhat[0] > 0.0) {
                    current_voxel = posx[current_voxel];
                } else {
                    current_voxel = negx[current_voxel];
                }
                tMaxX= tMaxX + tDeltaX;
            } else {
                if (uhat[2] > 0.0) {
                    current_voxel = posz[current_voxel];
                } else {
                    current_voxel = negz[current_voxel];
                }
                tMaxZ= tMaxZ + tDeltaZ;
            }
        } else {
            if (tMaxY < tMaxZ) {
                if (uhat[1] > 0.0) {
                    current_voxel = posy[current_voxel];
                } else {
                    current_voxel = negy[current_voxel];
                }
                tMaxY= tMaxY + tDeltaY;
            } else {
                if (uhat[2] > 0.0) {
                    current_voxel = posz[current_voxel];
                } else {
                    current_voxel = negz[current_voxel];
                }
                tMaxZ= tMaxZ + tDeltaZ;
            }
        }
        
//        if (hit_facet != -1) {
//            current_voxel = N3;
//        }
    }
    
    // If we made it here, it is because the ray left the grid without intersecting any facets
    // hit_facet will still be = -1, and hit_pnt will be empty
    return hit_facet;
}

int VoxelGrid::test_intersection(double* pt, double* uhat, vector<int>* facets, vector<double>* vertices, vector <Facet>* facetData, double* hit_pnt, int current_voxel, vector<int>* extraFacets)
{
    int hit_facet = -1;
    
    // Just need to find the first facet hit, so have to look at all of the
    // facets in the group, and select the shortest intersection if there are
    // more than one.
    
    // Pay no attention to which "side" of the triangle is intersected; just
    // pick the first one along the path
    
    // Note that I could change this to output t, which is the range to the
    // intersection along the ray
    
    int num_facets = fnum[current_voxel];
    int ii, current_facet;
    
    double tmin = 10.0*N*max(max(dx,dy),dz); // larger than the total voxel grid to start
    double pt_int[3], t, temp1, temp2;
    double* nhat;
    double* v1;
    double* v2;
    double* v3;
    bool hit;
    
    for (ii = 0; ii < num_facets; ii++) {
        
        // Solve for t in ray equation that gives a point in the triangle plane
        // dot((pt + t.*uhat) - v, normal) == 0 where v is one of the vertices in the triangle
        // dot(pt - v,normal) + t*dot(uhat,normal) == 0
        // t = dot(v - pt,normal)/dot(uhat,normal)
        current_facet = flist[current_voxel][ii];
        nhat = (*facetData)[current_facet].getNormal();
        v1 = &(*vertices)[(*facetData)[current_facet].getVert(0)];
        v2 = &(*vertices)[(*facetData)[current_facet].getVert(1)];
        v3 = &(*vertices)[(*facetData)[current_facet].getVert(2)];
        pt_int[0] = v1[0] - pt[0];
        pt_int[1] = v1[1] - pt[1];
        pt_int[2] = v1[2] - pt[2];
        temp1 = dot(&pt_int[0],nhat);
        temp2 = dot(uhat,nhat);
        t = temp1/temp2;
        
        // Ray is directed; don't travel backwards to an intersection
        if (t >= 0.0) {
            pt_int[0] = pt[0] + t*uhat[0];
            pt_int[1] = pt[1] + t*uhat[1];
            pt_int[2] = pt[2] + t*uhat[2];
            
            // test this point to see if its inside the current facet
            hit = in_triangle_3D(v1, v2, v3, &pt_int[0], nhat);
            
            if (hit) {
                if (t < tmin) {
                    hit_facet = current_facet;
                    hit_pnt[0] = pt_int[0];
                    hit_pnt[1] = pt_int[1];
                    hit_pnt[2] = pt_int[2];
                    tmin = t;
                }
                else if (abs(t - tmin) < 1e-14) {
                    extraFacets->push_back(current_facet);
                }
            }
        }
    }
    
    return hit_facet;
}

void VoxelGrid::inverse ( double* answer, double* x, double* y, double* n)
{
    double denom;
    
    denom = x[0]*(y[1]*n[2]-n[1]*y[2]) - y[0]*(x[1]*n[2]-n[1]*x[2]) + n[0]*(x[1]*y[2]-y[1]*x[2]);
    
    answer[0] = (y[1]*n[2]-n[1]*y[2])/denom;
    answer[1] = (y[2]*n[0]-n[2]*y[0])/denom;
    answer[2] = (y[0]*n[1]-n[0]*y[1])/denom;
    answer[3] = (x[2]*n[1]-n[2]*x[1])/denom;
    answer[4] = (x[0]*n[2]-n[0]*x[2])/denom;
    answer[5] = (x[1]*n[0]-n[1]*x[0])/denom;
    //    answer[6] = (y[2]*x[1]-x[2]*y[1])/denom;
    //    answer[7] = (y[0]*x[2]-x[0]*y[2])/denom;
    //    answer[8] = (y[1]*x[0]-x[1]*y[0])/denom;
    
    return;
}

bool VoxelGrid::in_triangle_3D( double *A, double *B, double *C, double *P, double *nhat)
{
    // In triangle stuff
    double v0[3], v1[3], v2[3], answer[6], uv[2];
    int i;
    
    for (i=0;  i<3; i++){
        v0[i] = C[i] - A[i];
        v1[i] = B[i] - A[i];
        v2[i] = P[i] - A[i];
    }
    
    inverse( &answer[0], &v1[0], &v0[0], nhat );
    
    uv[0] = answer[0]*v2[0] + answer[1]*v2[1] + answer[2]*v2[2];
    uv[1] = answer[3]*v2[0] + answer[4]*v2[1] + answer[5]*v2[2];
    
    return ((uv[0] >= -1e-14) && (uv[1] >= -1e-14) && (uv[0] + uv[1] <= (1.0+1e-14)));
}

double VoxelGrid::dot ( double* a, double* b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    
}

int VoxelGrid::getNumFinVox(int VoxNum)
{
    return fnum[VoxNum];
}

int VoxelGrid::getNumVinVox(int VoxNum)
{
    return vnum[VoxNum];
}

double VoxelGrid::getMaxDim()
{
    double x = 0.5*dx*N;
    double y = 0.5*dy*N;
    double z = 0.5*dz*N;
    
    return sqrt(x*x + y*y + z*z);
}

bool VoxelGrid::intersectTriVox(int voxNum, double* v1, double* v2, double* v3, double* nhat)
{
    // test for an intersection of the current voxel and facet using the separating axis theorem
    // returns true if there IS an intersection, false if no intersection
    double axis[3], rtemp;
    double hx = dx/2.0;
    double hy = dy/2.0;
    double hz = dz/2.0;
    double center[3], vc1[3], vc2[3], vc3[3];
    center[0] = xmin[voxNum] + hx;
    center[1] = ymin[voxNum] + hy;
    center[2] = zmin[voxNum] + hz;
    int ii;
    // translate facet to voxel centered frame
    for (ii = 0; ii<3; ii++) {
        vc1[ii] = v1[ii] - center[ii];
        vc2[ii] = v2[ii] - center[ii];
        vc3[ii] = v3[ii] - center[ii];
    }
    
    // First test the three voxel face normal directions, which are assumed to be axis aligned
    axis[0] = 1.0; axis[1] = 0.0; axis[2] = 0.0;
    if (satAABBTri(hx, &axis[0], &vc1[0], &vc2[0], &vc3[0])) {
        return false;
    }
    axis[0] = 0.0; axis[1] = 1.0; axis[2] = 0.0;
    if (satAABBTri(hy, &axis[0], &vc1[0], &vc2[0], &vc3[0])) {
        return false;
    }
    axis[0] = 0.0; axis[1] = 0.0; axis[2] = 1.0;
    if (satAABBTri(hz, &axis[0], &vc1[0], &vc2[0], &vc3[0])) {
        return false;
    }
    
    // Test facet normal
    rtemp = hx*abs(nhat[0]) + hy*abs(nhat[1]) + hz*abs(nhat[2]);
    if (satAABBTri(rtemp, nhat, &vc1[0], &vc2[0], &vc3[0])) {
        return false;
    }
    
    // Test facet edges crossed with voxel edges
    double f0[3], f1[3], f2[3], e[3];
    for (ii = 0; ii<3; ii++) {
        f0[ii] = v2[ii]-v1[ii];
        f1[ii] = v3[ii]-v2[ii];
        f2[ii] = v1[ii]-v3[ii];
    }
    
    e[0]= 0.0; e[1] = 0.0; e[2] = 0.0;
    for (int jj=0; jj<3; jj++) {
        e[jj] = 1.0;
        cross(&e[0], &f0[0], &axis[0]);
        rtemp = hx*abs(axis[0]) + hy*abs(axis[1]) + hz*abs(axis[2]);
        if (satAABBTri(rtemp, &axis[0], &vc1[0], &vc2[0], &vc3[0])) {
            return false;
        }
        e[jj] = 0.0;
    }
    for (int jj=0; jj<3; jj++) {
        e[jj] = 1.0;
        cross(&e[0], &f1[0], &axis[0]);
        rtemp = hx*abs(axis[0]) + hy*abs(axis[1]) + hz*abs(axis[2]);
        if (satAABBTri(rtemp, &axis[0], &vc1[0], &vc2[0], &vc3[0])) {
            return false;
        }
        e[jj] = 0.0;
    }
    for (int jj=0; jj<3; jj++) {
        e[jj] = 1.0;
        cross(&e[0], &f2[0], &axis[0]);
        rtemp = hx*abs(axis[0]) + hy*abs(axis[1]) + hz*abs(axis[2]);
        if (satAABBTri(rtemp, &axis[0], &vc1[0], &vc2[0], &vc3[0])) {
            return false;
        }
        e[jj] = 0.0;
    }
    
    // No separating axis found, so they intersect
    return true;
    
}

bool VoxelGrid::satAABBTri(double rBox, double* axis, double* v1, double* v2, double* v3)
{
    // returns TRUE if there IS a separating axis, indicating NO INTERSECTION
    
    double p1, p2, p3;
    p1 = dot(axis, v1);
    p2 = dot(axis, v2);
    p3 = dot(axis, v3);
    
    if ((min(p1,min(p2,p3)) > rBox) || (max(p1,max(p2,p3)) < -rBox)) {
        return true;
    } else {
        return false;
    }
}

void VoxelGrid::cross(double* a, double* b, double* c)
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
    return;
}
}