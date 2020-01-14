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

#include <algorithm>
#include <cmath>
#include <valarray>
#include <math.h>
#include <vector>
#include <tuple>

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
}
//-------------------------------------------------------//

vtkCurvatures1::~vtkCurvatures1() {
    // if (this->prinCurvature) {
    //     this->prinCurvature->Delete();
    // }
}

static void getPlane(double &a, double &b, double& c, double&d, const valarray<double> p, const valarray<double> n) {
    a = n[0];
    b = n[1];
    c = n[2];
    d = -(a*p[0] + b*p[1] + c*p[2]);
}

static void getBasisVectors(valarray<double>& b1, valarray<double>& b2, valarray<double>& b3, const valarray<double>& n) {
    double b1x, b1y, b1z;
    double b2x, b2y, b2z;
    double b3x, b3y, b3z;

    double nx = n[0];
    double ny = n[1];
    double nz = n[2];

    if(fabs(nx > TOLERANCE)) {
      b1x = -1.0 / nx * (ny + nz);
      b1y = b1z = 1.0;
    }
    else if(fabs(ny > TOLERANCE)) {
      b1x = b1z = 1.0;
      b1y = -1.0 / ny * (nx + nz);
    }
    else if(fabs(nz > TOLERANCE)) {
      b1x = b1y = 1.0;
      b1z = -1.0 / ny * (nx + ny);
    }

    b1 = {b1x, b1y, b1z};
    double mag = vtkMath::Norm(vtkCurvatures1::getArr(b1));
    if (mag == 0.0) mag = 1.0;
    b1 /= mag;

    b3 = {nx, ny, nz};

    double temp_b2[3];
    vtkMath::Cross(vtkCurvatures1::getArr(b3), vtkCurvatures1::getArr(b1), temp_b2);
    mag = vtkMath::Norm(temp_b2);
    if(mag == 0.0) mag = 1.0;
    b2 = vtkCurvatures1::getValArr(temp_b2);
    b2 /= mag;
}

void vtkCurvatures1::genNeighborhoods(vtkDataSet *inData, int ndepth) {

    // clear neighborhoods
    this->neighbors.clear();

    // initialize neighborhoods
    for (int i = 0; i < this->numPoints; i++) {
        this->neighbors.push_back(vtkSmartPointer<vtkIdList>::New());
    }

    int a, b, c;

    vtkPolyData *pd = vtkPolyData::SafeDownCast( inData );
    vtkIdType ncells = pd->GetNumberOfCells();
    for ( vtkIdType i = 0; i < ncells; i++) {

        // vtkIdType ctype = pd->GetCellType(i);
        //
        // if (ctype == VTK_POLYGON || ctype == VTK_TRIANGLE || ctype == VTK_LINE) {
        //     vtkIdType *pts;
        //     vtkIdType npts;
        //     pd->GetCellPoints(i, npts, pts);
        //     double cost;
        //
        //     for (int j = 0; j < npts; ++j) {
        //         vtkIdType u = pts[j];
        //         vtkIdType v = pts[(( j + 1 ) % npts)];
        //
        //         vtkSmartPointer<vtkIdList> mu = this->neighbors[u];
        //         if ( mu->IsId(v) == -1 ) {
        //             mu->InsertNextId(v);
        //         }
        //
        //         vtkSmartPointer<vtkIdList> mv = this->neighbors[v];
        //         if ( mv->IsId(u) == -1 ) {
        //             mv->InsertNextId(u);
        //         }
        //     }
        // }
        int npts;
        const vtkNew<vtkIdList> pts;
        pd->GetCellPoints(i, npts, pts);
        if (npts != 3) {
            cerr << "number of points that comprise cell is not 3" << endl;
        }
        a = pts[0];
        b = pts[1];
        c = pts[2];

        if (addNeighbor(this->neighbors[a], b, true))
            addNeighbor(this->neighbors[b], a, false));
        if (addNeighbor(this->neighbors[a], c, true))
            addNeighbor(this->neighbors[c], a, false));
        if (addNeighbor(this->neighbors[b], c, true))
            addNeighbor(this->neighbors[c], b, false));
    }

    // for ndepth > 1, add neighbors of neighbors
    for (int j = 1; j < ndepth; j++) {
      vector<vtkSmartPointer<vtkIdList>> bigNghbd = this->neighbors;
      for (i = 0; i < this->numPoints; i++) {
          int numn = this->neighbors[i]->GetNumberOfIds();
      }

      for (int k = 0; k < numn; k++) {
          vtkIdType nei = this->neighbors[i]->GetId(k);
          int newnumn = this->neighbors[nei]->GetNumberOfIds();
          for (int l = 0; l < newnumn; l++) {
              vtkIdType newnei = this->neighbors[nei]->GetId(l);
              addNeighbor(bigNghbd[i], newnei, true);
          }
      }

      this->neighbors.clear(); // do we need this?
      this->neighbors = bigNghbd;
    }
  //this->AdjacencyBuildTime.Modified();???
}
static bool addNeighbor(vtkSmartPointer<vtkIdList> list, vtkIdType nbr, bool check) {
    if (check) {
        if (list->IsId(nbr) != -1) list->InsertNextId(nbr);
        else return false;
    } else {
        list->InsertNextId(nbr);
    }
    return true;
}

void vtkCurvatures1::GetPrincipalCurvature(vtkPolyData *mesh, int ndepth, int dx, int dy, int dz) {
    // ensure that ndepth isn't too large or too small
    if (ndepth < 1 or ndepth > 1000) {
        ndepth = 2;
    }

    genNeighborhoods(mesh, ndepth); // need to write this function later
    // have check to see if neighborhoods generated properly
    if (!this->hasUnitNormals()) {
        this->genUnitNormals(mesh);
    }

    int maxn = getMaxNeighbors();
    maxn++;

    valarray<float> h(maxn);
    valarray<float> u(maxn);
    valarray<float> v(maxn);
    valarray<float> uu(maxn);
    valarray<float> vv(maxn);
    valarray<float> two_uv(maxn);

    float tmp[3],U[3][3],BU[3],work[3],workc[5];
    float rval;
    int   IPIV[4];
    float C[2][2],W[2];
    char  JOBZ, UPLO;
    int   INFO,  LDA, LDB, LWORK, N, NRHS;
    int   LDC, LWORKC, NC;
    int   i,j;

    mesh->BuildLinks();

    for (int idx = 0; idx < this->numPoints; i++) {
        double a, b, c, d;
        valarray<double> b1(3), b2(3), b3(3);

        int num_nei = this->neighbors[idx]->GetNumberOfIds();

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

        for (i = 0; i < num_nei; i++) {
            vtkIdType nei = this->neighbors[idx]->GetId(i);
            if (nei < 0) continue;

            double px, py, pz;
            double pt_nei[3];
            mesh->GetPoint(nei, pt_nei);

            px = pt_nei[0] * dx;
            py = pt_nei[1] * dy;
            pz = pt_nei[2] * dz;

            h[i] = a*px + b*py + c*pz + d;
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

        ssytrf_(&UPLO, &N, (float*)U, &LDA, IPIV, work, &LWORK, &INFO);

        if(INFO == 0) {
            ssytrs_(&UPLO, &N, &NRHS, (float*)U, &LDA, IPIV, BU,
              &LDB , &INFO);
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
    if (!this->hasNormals()) {
        this->genNormals(mesh);
    }
    for (int i = 0; i < this->numPoints; i++) {
        valarray<double> normal = this->normals[i];
        double norm = vtkMath::Norm(getArr(normal));
        if (norm != 0) {
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
    for (int i = 0; i < this->numPoints; i++) {
        // this->normals.push_back(make_tuple(0, 0, 0));
        valarray<double> init = {0, 0, 0};
        this->normals[i] = init;
    }
    mesh->BuildLinks();
    const int F = mesh->GetNumberOfCells();
    double p0[3];
    double p1[3];
    double p2[3];
    const vtkNew<vtkIdList> vertices;
    for (int f = 0; f < F; f++) {
        // need to reset this each time?
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

        double cross[3];
        // kinda messy - is there a better way to do this?
        vtkMath::Cross(getArr(getValArr(p1) - getValArr(p0)),
            getArr(getValArr(p2) - getValArr(p0)), cross);
        this->normals[v0] += getValArr(cross);

        vtkMath::Cross(getArr(getValArr(p2)-getValArr(p1)),
            getArr(getValArr(p0) -getValArr(p1)), cross);
        this->normals[v1] += getValArr(cross);


        vtkMath::Cross(getArr(getValArr(p0) - getValArr(p2)),
            getArr(getValArr(p1) - getValArr(p2)), cross);
        this->normals[v2] += getValArr(cross);
    }

    for (int i = 0; i < this->numPoints; i++) {
        this->normals[i] = this->normals[i]/6.0;
    }
}

static double* getArr(const valarray<double> temp) {
    static double ans[3];
    for (int i = 0; i < 3; i++) {
        ans[i] = temp[i];
    }
    return ans;
}

static valarray<double> getValArr(const double temp[3]) {
    valarray<double> ans;
    for (int i = 0; i < 3; i++) {
        ans[i] = temp[i];
    }
    return ans;
}

static double checkCurv(double curv) {
    if (fabs(curv) > MAX_CURV) {
        if (curv < 0) curv = -MAX_CURV;
        else curv = MAX_CURV;
    }
    return curv;
}
void vtkCurvatures1::GetMeanCurvature(vtkPolyData *mesh) {
    if (prinCurvature.size() == 0) {
        this->GetPrincipalCurvature(mesh, 2, 1, 1, 1);
    }
    const vtkNew<vtkDoubleArray> meanCurvature;
    meanCurvature->SetName("Mean_Curvature");
    meanCurvature->SetNumberOfComponents(1);
    meanCurvature->SetNumberOfTuples(this->numPoints);
    // Get the array so we can write to it directly
    double *meanCurvatureData = meanCurvature->GetPointer(0);

    for (int i = 0; i < this->prinCurvature.size(); i++) {
        double temp = (this->prinCurvature[i].first + this->prinCurvature[i].second)/2;
        meanCurvatureData[i] = -1 * checkCurv(temp);
    }

    mesh->GetPointData()->AddArray(meanCurvature);
    mesh->GetPointData()->SetActiveScalars("Mean_Curvature");

    vtkDebugMacro("Set Values of Mean Curvature: Done");

    if (meanCurvature) meanCurvature->Delete();
}
//--------------------------------------------
void vtkCurvatures1::GetGaussCurvature(vtkPolyData *mesh) {
    if (prinCurvature.size() == 0) {
        this->GetPrincipalCurvature(mesh, 2, 1, 1, 1);
    }
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

    if (gaussCurvature) gaussCurvature->Delete();
}

void vtkCurvatures1::GetMaximumCurvature(vtkPolyData *mesh) {
    if (prinCurvature.size() == 0) {
        this->GetPrincipalCurvature(mesh, 2, 1, 1, 1);
    }
    const vtkNew<vtkDoubleArray> maxCurvature;
    maxCurvature->SetName("Maximum_Curvature");
    maxCurvature->SetNumberOfComponents(1);
    maxCurvature->SetNumberOfTuples(this->numPoints);
    // Get the array so we can write to it directly
    double *maxCurvatureData = maxCurvature->GetPointer(0);

    for (int i = 0; i < this->prinCurvature.size(); i++) {
        double temp = max(fabs(this->prinCurvature[i].first), fabs(this->prinCurvature[i].second));
        maxCurvatureData[i] = -1 * checkCurv(temp);

    }

    mesh->GetPointData()->AddArray(maxCurvature);
    mesh->GetPointData()->SetActiveScalars("Maximum_Curvature");

    vtkDebugMacro("Set Values of Maximum Curvature: Done");

    if (maxCurvature) maxCurvature->Delete();
}

void vtkCurvatures1::GetMinimumCurvature(vtkPolyData *mesh) {
    if (prinCurvature.size() == 0) {
        this->GetPrincipalCurvature(mesh, 2, 1, 1, 1);
    }
    const vtkNew<vtkDoubleArray> minCurvature;
    minCurvature->SetName("Minimum_Curvature");
    minCurvature->SetNumberOfComponents(1);
    minCurvature->SetNumberOfTuples(this->numPoints);
    // Get the array so we can write to it directly
    double *minCurvatureData = minCurvature->GetPointer(0);

    for (int i = 0; i < this->prinCurvature.size(); i++) {
        double temp = min(fabs(this->prinCurvature[i].first), fabs(this->prinCurvature[i].second));
        minCurvatureData[i] = -1 * checkCurv(temp);

    }

    mesh->GetPointData()->AddArray(minCurvature);
    mesh->GetPointData()->SetActiveScalars("Minimum_Curvature");

    vtkDebugMacro("Set Values of Maximum Curvature: Done");

    if (minCurvature) minCurvature->Delete();
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

    //-------------------------------------------------------//
    //    Set Curvatures as PointData  Scalars               //
    //-------------------------------------------------------//

    // eventually these will be global variables that people can configure
    // via property panel?

    this->numPoints = output->GetNumberOfPoints();
    this->numPolys = output->GetNumberOfPolys();

    int ndepth = 2;
    int dx = 1;
    int dy = 1;
    int dz = 1;

    this->GetPrincipalCurvature(output, ndepth, dx, dy, dz);

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
