/*=========================================================================
  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkCurvatures1.cxx,v $
  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.
     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.
=========================================================================*/
#include "vtkCurvatures1.h"

#include "vtkCellArray.h"
#include "vtkDoubleArray.h"
#include "vtkFieldData.h"
#include "vtkFloatArray.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataNormals.h"
#include "vtkPolygon.h"
#include "vtkTriangle.h"
#include "vtkIdList.h"

#include <algorithm>
#include <cmath>
#include <valarray>
#include <math.h>
#include <vector>
#include <tuple>
#include <iostream>
#include <fstream>

using std::valarray;
using std::vector;
using std::tuple;
using std::make_tuple;
using std::make_pair;
using std::max;
using std::min;

vtkStandardNewMacro(vtkCurvatures1);

extern "C" {
  extern void ssytrf_(char *uplo, int *n, float *a, int *lda, int *ipiv,
		      float *work, int *lwork, int *info);
  extern void ssytrs_(char *uplo, int *n, int *nrhs, float *a, int *lda,
		      int *ipiv, float *b, int *ldb, int *info);
  extern void ssyev_(char *jobz, char *uplo, int *n, float *a, int *lda,
		     float *w, float *work, int *lwork, int *info);
};

//------------------------------------------------------------------------------
#if VTK3
vtkCurvatures1* vtkCurvatures1::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkCurvatures1");
  if(ret)
    {
    return (vtkCurvatures1*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkCurvatures1;
}
#endif
//-------------------------------------------------------//
vtkCurvatures1::vtkCurvatures1()
{
  this->CurvatureType = 0;
  this->InvertMeanCurvature = 0;
  this->ndepth = 2;
  this->dx = 1;
  this->dy = 1;
  this->dz = 1;
}
//-------------------------------------------------------//

vtkCurvatures1::~vtkCurvatures1() {
}

void vtkCurvatures1::getPlane(double &a, double &b, double& c, double&d, const valarray<double> p, const valarray<double> n) {
    a = n[0];
    b = n[1];
    c = n[2];
    d = -(a*p[0] + b*p[1] + c*p[2]);
}

void vtkCurvatures1::getBasisVectors(valarray<double>& b1, valarray<double>& b2, valarray<double>& b3, const valarray<double>& n) {
    double b1x, b1y, b1z;
    double b2x, b2y, b2z;
    double b3x, b3y, b3z;

    double nx = n[0];
    double ny = n[1];
    double nz = n[2];

    if(fabs(nx) > TOLERANCE) {
      b1x = -1.0 / nx * (ny + nz);
      b1y = b1z = 1.0;
    }
    else if(fabs(ny) > TOLERANCE) {
      b1x = b1z = 1.0;
      b1y = -1.0 / ny * (nx + nz);
    }
    else if(fabs(nz) > TOLERANCE) {
      b1x = b1y = 1.0;
      b1z = -1.0 / ny * (nx + ny);
    }

    b1 = {b1x, b1y, b1z};
    double mag = vtkCurvatures1::myNorm(b1);
    if (mag == 0.0) mag = 1.0;
    b1 /= mag;

    b3 = {nx, ny, nz};

    valarray<double> temp_b2(3);
    vtkCurvatures1::myCross(b3, b1, temp_b2);
    mag = vtkCurvatures1::myNorm(temp_b2);
    if(mag == 0.0) mag = 1.0;
    b2 = {temp_b2[0], temp_b2[1], temp_b2[2]};
    b2 /= mag;
}

void vtkCurvatures1::genNeighborhoods(vtkDataSet *inData, int ndepth) {

    int i;
    vtkIdType a, b, c;
    if (ndepth <= 0) return;

    // clear neighborhoods
    this->neighbors.clear();

    // initialize neighborhoods with new vtkIdLists
    for (int i = 0; i < this->numPoints; i++) {
        this->neighbors.push_back(vtkSmartPointer<vtkIdList>::New());
    }

    vtkPolyData *pd = vtkPolyData::SafeDownCast( inData );
    vtkIdType ncells = pd->GetNumberOfCells();

    // iterate through cells on mesh
    for (i = 0; i < ncells; i++) {
        vtkIdType npts;
        const vtkIdType* pts;

        // extract vertices that comprise current cell
        pd->GetCellPoints(i, npts, pts);
        if (npts != 3) {
            cerr << "number of points that comprise cell is not 3" << endl;
        }
        a = pts[0];
        b = pts[1];
        c = pts[2];

        if (addNeighbor(this->neighbors[a], b, true))  // check for duplicates
            addNeighbor(this->neighbors[b], a, false); // don't need to check
        if (addNeighbor(this->neighbors[a], c, true))  // ...
            addNeighbor(this->neighbors[c], a, false); // ...
        if (addNeighbor(this->neighbors[b], c, true))  // ...
            addNeighbor(this->neighbors[c], b, false); // ...
    }

    // for ndepth > 1, add neighbors of neighbors
    for (int j = 1; j < ndepth; j++) {
        // create and initialize new bigger neighborhood
        vector<vtkSmartPointer<vtkIdList>> bigNghbd;
        copyNeighbors(this->neighbors, bigNghbd);

        // iterate through vertices
        for (i = 0; i < this->numPoints; i++) {
          int numn = this->neighbors[i]->GetNumberOfIds();

          // iterate through each neighbor of ith vertex
          for (int k = 0; k < numn; k++) {
              vtkIdType nei = this->neighbors[i]->GetId(k);
              int newnumn = this->neighbors[nei]->GetNumberOfIds();

              // iterate through neighbors of current neighbor (k)
              // and add to neighbor list of vertex i
              for (int l = 0; l < newnumn; l++) {
                  vtkIdType newnei = this->neighbors[nei]->GetId(l);
                  addNeighbor(bigNghbd[i], newnei, true);
              }
          }
        }

        this->neighbors.clear(); // do we need this?
        this->neighbors = bigNghbd;
    }
}

/* deep copy contents of orig into copy */
void vtkCurvatures1::copyNeighbors(vector<vtkSmartPointer<vtkIdList>> orig,
    vector<vtkSmartPointer<vtkIdList>>& copy) {
        for (int i = 0; i < this->numPoints; i++) {
            copy.push_back(vtkSmartPointer<vtkIdList>::New());
            for (int j = 0; j < orig[i]->GetNumberOfIds(); j++) {
                copy[i]->InsertNextId(orig[i]->GetId(j));
            }
        }
    }

bool vtkCurvatures1::addNeighbor(vtkSmartPointer<vtkIdList>& list, vtkIdType nbr, bool check) {
    // if we need to check for duplicates
    if (check) {
        if (list->IsId(nbr) == -1) {    // if vertex ID is not already in list
            list->InsertNextId(nbr);
        } else {                        // if vertex ID is  already in list
            return false;
        }
    } else {
        list->InsertNextId(nbr);
    }
    return true;
}

void vtkCurvatures1::GetPrincipalCurvature(vtkPolyData *mesh, int ndepth, int dx, int dy, int dz) {

    this->prinCurvature.clear();

    // ensure that ndepth isn't too large or too small
    if (ndepth < 1 || ndepth > 1000) {
        ndepth = 2;
    }

    mesh->BuildLinks();

    // generate neighborhood of depth ndepth for each vertex
    genNeighborhoods(mesh, ndepth);

    // if currently no unit normals, then generate them
    if (!this->hasUnitNormals()) {
        this->genUnitNormals(mesh);
        if (!this->hasUnitNormals()) {
            cerr << "ERROR: cannot generate neighborhoods" << endl;
            return;
        }
    }

    // get max number of neighbors to initialize the size of the following arrays
    int maxn = getMaxNeighbors();

    maxn++;

    // used for shape operator calculation
    valarray<float> h(maxn);
    valarray<float> u(maxn);
    valarray<float> v(maxn);
    valarray<float> uu(maxn);
    valarray<float> vv(maxn);
    valarray<float> two_uv(maxn);

    // used for shape operator calculation
    // and eigenvalue calculation via lapack libraries
    float tmp[3],U[3][3],BU[3],work[3],workc[5];
    float rval;
    int   IPIV[4];
    float C[2][2],W[2];
    char  JOBZ, UPLO;
    int   INFO,  LDA, LDB, LWORK, N, NRHS;
    int   LDC, LWORKC, NC;
    int   i,j;

    // iterate through each vertex on surface and compute principal curvatures
    for (int idx = 0; idx < this->numPoints; idx++) {

        // initialize variables to store tangent plane parameters
        double a, b, c, d;

        // initialize arrays to store plane basis vectors
        valarray<double> b1(3);
        valarray<double> b2(3);
        valarray<double> b3(3);

        // extract x, y, z coordinates for current idx
        double ptx, pty, ptz;
        double pt[3];
        mesh->GetPoint(idx, pt);

        ptx = pt[0] * dx;
        pty = pt[1] * dy;
        ptz = pt[2] * dz;

        valarray<double> pp = {ptx, pty, ptz};
        valarray<double> nn = this->unitNormals[idx];

        vtkCurvatures1::getPlane(a, b, c, d, pp, nn);
        vtkCurvatures1::getBasisVectors(b1, b2, b3, nn);

        double nx = nn[0];
        double ny = nn[1];
        double nz = nn[2];

        int num_nei = this->neighbors[idx]->GetNumberOfIds();

        // iterate through neighbors
        for (i = 0; i < num_nei; i++) {
            vtkIdType nei = this->neighbors[idx]->GetId(i);
            if (nei < 0) continue;

            // extract x, y, z coordinates of ith neighbor
            double px, py, pz;
            double pt_nei[3];
            mesh->GetPoint(nei, pt_nei);

            px = pt_nei[0] * dx;
            py = pt_nei[1] * dy;
            pz = pt_nei[2] * dz;

            // h[i] is the distance from the neighbor point to the tangent plane
            // defined by original point normal
            h[i] = a*px + b*py + c*pz + d;

            // tmp holds x, y, z coords of neighbor point projected onto tangent
            // plane - original point
            tmp[0] = px - h[i]*nx - ptx;
            tmp[1] = py - h[i]*ny - pty;
            tmp[2] = pz - h[i]*nz - ptz;

            u[i] = tmp[0]*b1[0] + tmp[1]*b1[1] + tmp[2]*b1[2];
            v[i] = tmp[0]*b2[0] + tmp[1]*b2[1] + tmp[2]*b2[2];

            two_uv[i] = 2.0*u[i]*v[i];
            uu[i] = u[i]*u[i];
            vv[i] = v[i]*v[i];
        }

        // initialize BU and U to 0
        for (i = 0; i < 3; i++) {
            BU[i] = 0;
            for (j = 0; j < 3; j++) {
                U[i][j] = 0;
            }
        }

        // create matrix to plug into equation solver
        for(i=0; i < num_nei; i++) {
            U[0][0] += (uu[i]*uu[i]);
            U[0][1] += (uu[i]*two_uv[i]);
            U[0][2] += (uu[i]*vv[i]);
            U[1][1] += (two_uv[i]*two_uv[i]);
            U[1][2] += (two_uv[i]*vv[i]);
            U[2][2] += (vv[i]*vv[i]);
            BU[0]   += (uu[i]*2.0*h[i]);
            BU[1]   += (two_uv[i]*2.0*h[i]);
            BU[2]   += (vv[i]*2.0*h[i]);
        }

        // matrix is symmetric
        U[1][0] = U[0][1];
        U[2][0] = U[0][2];
        U[2][1] = U[1][2];

        UPLO  = 'L';
        INFO  = 1;
        LDA   = 3;
        LDB   = 3;
        LWORK = 3;
        N     = 3;
        NRHS  = 1;

        // compute factorization of real symmetric matrix U
        ssytrf_(&UPLO, &N, (float*)U, &LDA, IPIV, work, &LWORK, &INFO);

        if(INFO == 0) {

            // solve system of linear equations to obtain shape operator
            ssytrs_(&UPLO, &N, &NRHS, (float*)U, &LDA, IPIV, BU,
              &LDB , &INFO);

            // matrix C is the shape operator for this vertex
            C[0][0] = BU[0];
            C[0][1] = BU[1];
            C[1][0] = BU[1];
            C[1][1] = BU[2];

            /* get the eigenvalue of the matrix c */
            JOBZ   = 'N';
            LDC    = 2;
            NC     = 2;
            LWORKC = 3*NC-1;
            ssyev_(&JOBZ, &UPLO, &NC, (float*)C, &LDC, W, workc,
             &LWORKC,&INFO);

            this->prinCurvature.push_back(make_pair(W[0], W[1]));

        } else {
            this->prinCurvature.push_back(make_pair(0.0, 0.0));
            cerr << "ERROR: systrs failed in GetPrincipalCurvature()" << endl;
        }
    }
}

int vtkCurvatures1::getMaxNeighbors() {
    int maxn = 0;
    for (int m = 0; m < this->numPoints; m++) {
        if (this->neighbors[m]->GetNumberOfIds() > maxn) {
            maxn = this->neighbors[m]->GetNumberOfIds();
        }
    }
    return maxn;
}

bool vtkCurvatures1::hasNormals() {
    return this->normals.size() > 0;
}

bool vtkCurvatures1::hasUnitNormals() {
    return this->unitNormals.size() > 0;
}

void vtkCurvatures1::genUnitNormals(vtkPolyData* mesh) {
    // check if normals exist for this surface
    if (!this->hasNormals()) {
        // if not, then generate normals first
        this->genNormals(mesh);
    }
    // iterate through vertices
    for (int i = 0; i < this->numPoints; i++) {
        valarray<double> normal = this->normals[i];
        double norm = myNorm(normal);
        if (norm != 0) {    // if norm is not 0, then normalize the normal vector
            double x = normal[0]/norm;
            double y = normal[1]/norm;
            double z = normal[2]/norm;
            valarray<double> unitNorm = {x, y, z};
            this->unitNormals.push_back(unitNorm);
        } else {
            this->unitNormals.push_back(normal);
        }
    }
}

void vtkCurvatures1::genNormals(vtkPolyData* mesh) {

    // initializing normals array
    for (int i = 0; i < this->numPoints; i++) {
        valarray<double> init = {0, 0, 0};
        this->normals.push_back(init);
    }
    mesh->BuildLinks();
    const int F = mesh->GetNumberOfCells();

    // iterate through triangles
    // and accumulate the contribution from that triangle to each of its vertices
    for (int f = 0; f < F ; f++) {
        vtkNew<vtkIdList> vertices;
        double p0[3];
        double p1[3];
        double p2[3];

        // gets three points that comprise triangle and stores in vertices
        mesh->GetCellPoints(f, vertices);

        // get vertex IDs that comprise the current cell
        int v0 = vertices->GetId(0);
        int v1 = vertices->GetId(1);
        int v2 = vertices->GetId(2);

        // get the x,y,z coordinates that correspond to each vertex
        mesh->GetPoint(v0, p0);
        mesh->GetPoint(v1, p1);
        mesh->GetPoint(v2, p2);

        // create valarray vectors for each point so that we can use
        // element-wise operations
        valarray<double> p0_1 = {p0[0], p0[1], p0[2]};
        valarray<double> p1_1 = {p1[0], p1[1], p1[2]};
        valarray<double> p2_1 = {p2[0], p2[1], p2[2]};

        valarray<double> cross0(3);
        myCross(p1_1 - p0_1, p2_1 - p0_1, cross0);
        this->normals[v0] += cross0;

        valarray<double> cross1(3);
        myCross(p2_1 - p1_1, p0_1 - p1_1, cross1);
        this->normals[v1] += cross1;

        valarray<double> cross2(3);
        myCross(p0_1 - p2_1, p1_1 - p2_1, cross2);
        this->normals[v2] += cross2;
    }

    for (int i = 0; i < this->numPoints; i++) {
        this->normals[i] = this->normals[i] / 6.0;
    }
}

double vtkCurvatures1::myNorm(valarray<double> temp) {
    return sqrt(temp[0]*temp[0] + temp[1]*temp[1] + temp[2]*temp[2]);
}

 void vtkCurvatures1::myCross(valarray<double> a, valarray<double> b, valarray<double>& c) {
    double Cx = a[1]*b[2] - a[2]*b[1];
    double Cy = a[2]*b[0] - a[0]*b[2];
    double Cz = a[0]*b[1] - a[1]*b[0];
    c[0] = Cx;
    c[1] = Cy;
    c[2] = Cz;
 }

double vtkCurvatures1::checkCurv(double curv) {
    if (fabs(curv) > MAX_CURV) {
        if (curv < 0) curv = -MAX_CURV;
        else curv = MAX_CURV;
    }
    return curv;
}

void vtkCurvatures1::GetMeanCurvature(vtkPolyData *mesh) {
    this->GetPrincipalCurvature(mesh, this->ndepth, this->dx, this->dy, this->dz);

    const vtkNew<vtkDoubleArray> gaussCurvature;
    gaussCurvature->SetName("Mean_Curvature");
    gaussCurvature->SetNumberOfComponents(1);
    gaussCurvature->SetNumberOfTuples(this->numPoints);
    // Get the array so we can write to it directly
    double *gaussCurvatureData = gaussCurvature->GetPointer(0);

    for (int i = 0; i < this->prinCurvature.size(); i++) {
        double temp = (this->prinCurvature[i].first + this->prinCurvature[i].second)/2;
        gaussCurvatureData[i] = -1 * checkCurv(temp);
    }

    mesh->GetPointData()->AddArray(gaussCurvature);
    mesh->GetPointData()->SetActiveScalars("Mean_Curvature");

    vtkDebugMacro("Set Values of Gauss Curvature: Done");
}
//--------------------------------------------
void vtkCurvatures1::GetGaussCurvature(vtkPolyData *mesh) {
    this->GetPrincipalCurvature(mesh, this->ndepth, this->dx, this->dy, this->dz);
    const vtkNew<vtkDoubleArray> gaussCurvature;
    gaussCurvature->SetName("Gauss_Curvature");
    gaussCurvature->SetNumberOfComponents(1);
    gaussCurvature->SetNumberOfTuples(this->numPoints);
    // Get the array so we can write to it directly
    double *gaussCurvatureData = gaussCurvature->GetPointer(0);

    for (int i = 0; i < this->prinCurvature.size(); i++) {
        double temp = (this->prinCurvature[i].first * this->prinCurvature[i].second);
        gaussCurvatureData[i] = -1 * checkCurv(temp);
    }

    mesh->GetPointData()->AddArray(gaussCurvature);
    mesh->GetPointData()->SetActiveScalars("Gauss_Curvature");

    vtkDebugMacro("Set Values of Gauss Curvature: Done");
}

void vtkCurvatures1::GetMaximumCurvature(vtkPolyData *mesh) {
    this->GetPrincipalCurvature(mesh, this->ndepth, this->dx, this->dy, this->dz);
    const vtkNew<vtkDoubleArray> maxCurvature;
    maxCurvature->SetName("Maximum_Curvature");
    maxCurvature->SetNumberOfComponents(1);
    maxCurvature->SetNumberOfTuples(this->numPoints);
    // Get the array so we can write to it directly
    double *maxCurvatureData = maxCurvature->GetPointer(0);

    for (int i = 0; i < this->prinCurvature.size(); i++) {
        double temp;
        if (fabs(this->prinCurvature[i].first) > fabs(this->prinCurvature[i].second)) {
            temp = this->prinCurvature[i].first;
        } else {
            temp = this->prinCurvature[i].second;
        }
        maxCurvatureData[i] = -1 * checkCurv(temp);
    }

    mesh->GetPointData()->AddArray(maxCurvature);
    mesh->GetPointData()->SetActiveScalars("Maximum_Curvature");

    vtkDebugMacro("Set Values of Maximum Curvature: Done");
}

void vtkCurvatures1::GetMinimumCurvature(vtkPolyData *mesh) {
    this->GetPrincipalCurvature(mesh, this->ndepth, this->dx, this->dy, this->dz);
    const vtkNew<vtkDoubleArray> minCurvature;
    minCurvature->SetName("Minimum_Curvature");
    minCurvature->SetNumberOfComponents(1);
    minCurvature->SetNumberOfTuples(this->numPoints);
    // Get the array so we can write to it directly
    double *minCurvatureData = minCurvature->GetPointer(0);

    for (int i = 0; i < this->prinCurvature.size(); i++) {
        double temp;
        if (fabs(this->prinCurvature[i].first) < fabs(this->prinCurvature[i].second)) {
            temp = this->prinCurvature[i].first;
        } else {
            temp = this->prinCurvature[i].second;
        }
        minCurvatureData[i] = -1 * checkCurv(temp);

    }

    mesh->GetPointData()->AddArray(minCurvature);
    mesh->GetPointData()->SetActiveScalars("Minimum_Curvature");

    vtkDebugMacro("Set Values of Maximum Curvature: Done");
}

//-------------------------------------------------------
int vtkCurvatures1::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
    // get the info objects
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    // get the input and ouptut
    vtkPolyData *input = vtkPolyData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

    // Null input check
    if (!input) return 0;

    output->CopyStructure(input);
    output->GetPointData()->PassData(input->GetPointData());
    output->GetFieldData()->PassData(input->GetFieldData());

    this->numPoints = output->GetNumberOfPoints();

    //-------------------------------------------------------//
    //    Set Curvatures as PointData  Scalars               //
    //-------------------------------------------------------//

    if ( this->CurvatureType == VTK_CURVATURE_GAUSS ) {
        this->GetGaussCurvature(output);
    }
    else if ( this->CurvatureType == VTK_CURVATURE_MEAN ) {
        this->GetMeanCurvature(output);
    }
    else if ( this->CurvatureType ==  VTK_CURVATURE_MAXIMUM ) {
        this->GetMaximumCurvature(output);
    }
    else if ( this->CurvatureType ==  VTK_CURVATURE_MINIMUM ) {
        this->GetMinimumCurvature(output);
    }
    else {
        vtkErrorMacro("Only Gauss, Mean, Max, and Min Curvature type available");
        return 1;
    }

    return 1;
}
/*-------------------------------------------------------*/
void vtkCurvatures1::PrintSelf(ostream& os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os,indent);
    os << indent << "CurvatureType: " << this->CurvatureType << "\n";
    os << indent << "InvertMeanCurvature: " << this->InvertMeanCurvature << "\n";
}
