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

#include <cmath>
#include <valarray>
#include <math.h>
#include <vector>
#include <tuple>

using std::valarray;
using std::vector;
using std::tuple;
using std::make_tuple;

vtkStandardNewMacro(vtkCurvatures1);

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
  this->prinCurvature = vtkDoubleArray::New();
  // this->normals = vtkDoubleArray::New();
  // this->unitNormals = vtkDoubleArray::New();
  // initialize valarrays?
}
//-------------------------------------------------------//

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

void vtkCurvatures1::GetPrincipalCurvature(vtkPolyData *mesh, int ndepth, int dx, int dy, int dz) {
    // ensure that ndepth isn't too large or too small
    if (ndepth < 1 or ndepth > 1000) {
        ndepth = 2;
    }
    // genNeighborhoods(ndepth); // need to write this function later
    // have check to see if neighborhoods generated properly
    if (!this->hasUnitNormals()) {
        this->genUnitNormals(mesh);
    }
    this->prinCurvature->SetNumberOfTuples(this->numPoints);

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

        for (int i = 0; i < num_nei; i++) {
            int nei = this->neighbors[idx]->GetId(i);
            if (nei < 0) continue;

            double px, py, pz;
            double pt_nei[3];
            mesh->GetPoint(nei, pt_nei);
            
            px = pt_nei[0] * dx;
            py = pt_nei[1] * dy;
            pz = pt_nei[2] * dz;

            h[i] = a*px + b*py + c*pz + d;
            tmp[0] = px - h[i] * nx - ptx;
            tmp[1] = py - h[i] * ny - pty;
            tmp[2] = pz - h[i] * nz - ptz;

            u[i] = tmp[0]*b1.x()+tmp[1]*b1.y()+tmp[2]*b1.z();
            v[i] = tmp[0]*b2.x()+tmp[1]*b2.y()+tmp[2]*b2.z();

            two_uv[i] = 2.0*u[i]*v[i];
            uu[i] = u[i]*u[i];
            vv[i] = v[i]*v[i];

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
        // vtkMath::cross(diff(p1, p0), diff(p2, p0), cross);
        vtkMath::Cross(getArr(getValArr(p1) - getValArr(p0)),
            getArr(getValArr(p2) - getValArr(p0)), cross);

        // kinda messy - is there a better way to do this?
        // this->normals[v0] = getTuple(sum(getArray(this->normals[v0]), cross));
        this->normals[v0] += getValArr(cross);

        // vtkMath::cross(diff(p2, p1), diff(p0, p1), cross);
        vtkMath::Cross(getArr(getValArr(p2)-getValArr(p1)),
            getArr(getValArr(p0) -getValArr(p1)), cross);
        // this->normals[v1] = getTuple(sum(getArray(this->normals[v1]), cross));
        this->normals[v1] += getValArr(cross);

        // vtkMath::cross(diff(p0, p2), diff(p1, p2), cross);
        vtkMath::Cross(getArr(getValArr(p0) - getValArr(p2)),
            getArr(getValArr(p1) - getValArr(p2)), cross);
        // this->normals[v2] = getTuple(sum(getArray(this->normals[v1]), cross));
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
    // return valarray<double> ans(temp, 3);
}

// static tuple getTuple(double x[3]) {
//     return make_tuple(double[0], double[1], double[2]);
// }
//
// static double[3] getArray(tuple t) {
//     double ans[3];
//     for (int i = 0; i < 3; i++) {
//         ans[i] = get<i>(tuple);
//     }
//     return ans;
// }

// static double[3] diff(double x[3], double y[3]) {
//     double ans[3];
//     for (int i = 0; i < 3; i++) {
//         ans[i] = x[i] - y[i];
//     }
//     return ans;
// }
//
// static double[3] sum(double x[3], double y[3]) {
//     double ans[3];
//     for (int i = 0; i < 3; i++) {
//         ans[i] = x[i] + y[i];
//     }
//     return ans;
// }

// void vtkCurvatures1::getNorm(tuple temp) {
//     double x = get<0>(temp); double y = get<1>(temp); double z = get<2>(temp);
//     return sqrt(x*x + y*y + z*z);
// }


void vtkCurvatures1::GetMeanCurvature(vtkPolyData *mesh)
{
    cout << "in mean curv" << endl;
    vtkDebugMacro("Start vtkCurvatures1::GetMeanCurvature");

    // Empty array check
    if (mesh->GetNumberOfPolys()==0 || mesh->GetNumberOfPoints()==0)
      {
      vtkErrorMacro("No points/cells to operate on");
      return;
      }

    int numPts = mesh->GetNumberOfPoints();

    // vtkData
    vtkIdList* vertices, *vertices_n, *neighbours;

    vtkTriangle* facet;
    vtkTriangle* neighbour;
    //     create-allocate
    vertices   = vtkIdList::New();
    vertices_n = vtkIdList::New();
    neighbours = vtkIdList::New();
    facet      = vtkTriangle::New();
    neighbour  = vtkTriangle::New();
    vtkDoubleArray* meanCurvature = vtkDoubleArray::New();
    meanCurvature->SetName("Mean_Curvature");
    meanCurvature->SetNumberOfComponents(1);
    meanCurvature->SetNumberOfTuples(numPts);
    // Get the array so we can write to it directly
    double *meanCurvatureData = meanCurvature->GetPointer(0);
    //     data
    int v, v_l, v_r, v_o,  f, F, n, nv;// n short for neighbor

    //     create-allocate
    double n_f[3]; // normal of facet (could be stored for later?)
    double n_n[3]; // normal of edge
    double t[3];   // to store the cross product of n_f n_n
    double ore[3]; // origin of e
    double end[3]; // end of e
    double oth[3]; //     third vertex necessary for comp of n
    double vn0[3];
    double vn1[3]; // vertices for computation of neighbour's n
    double vn2[3];
    double e[3];   // edge (oriented)

    double cs, sn;    // cs: cos; sn sin
    double angle, length, Af, Hf;  // temporary store

    mesh->BuildLinks();
    //data init
    f = 0;
    F = mesh->GetNumberOfCells();
    // init, preallocate the mean curvature
    int* num_neighb = new int[numPts];
    for (v = 0; v < numPts; v++)
      {
      meanCurvatureData[v] = 0.0;
      num_neighb[v] = 0;
      }

    //     main loop
    vtkDebugMacro(<<"Main loop: loop over facets such that id > id of neighb");
    vtkDebugMacro(<<"so that every edge comes only once");
    //
    for (f = 0; f < F; f++)
      {
      mesh->GetCellPoints(f,vertices);
      nv = vertices->GetNumberOfIds();

      for (v = 0; v < nv; v++)
        {
        // get neighbour
        v_l = vertices->GetId(v);
        v_r = vertices->GetId((v+1) % nv);
        v_o = vertices->GetId((v+2) % nv);
        mesh->GetCellEdgeNeighbors(f,v_l,v_r,neighbours);

        // compute only if there is really ONE neighbour
        // AND meanCurvature has not been computed yet!
        // (ensured by n > f)
        if (neighbours->GetNumberOfIds() == 1 && (n = neighbours->GetId(0)) > f)
          {
          // find 3 corners of f: in order!
          mesh->GetPoint(v_l,ore);
          mesh->GetPoint(v_r,end);
          mesh->GetPoint(v_o,oth);
          // compute normal of f
          facet->ComputeNormal(ore,end,oth,n_f);
          // compute common edge
          e[0] = end[0]; e[1] = end[1]; e[2] = end[2];
          e[0] -= ore[0]; e[1] -= ore[1]; e[2] -= ore[2];
          length = double(vtkMath::Normalize(e));
          Af = double(facet->TriangleArea(ore,end,oth));
          // find 3 corners of n: in order!
          mesh->GetCellPoints(n,vertices_n);
          mesh->GetPoint(vertices_n->GetId(0),vn0);
          mesh->GetPoint(vertices_n->GetId(1),vn1);
          mesh->GetPoint(vertices_n->GetId(2),vn2);
          Af += double(facet->TriangleArea(vn0,vn1,vn2));
          // compute normal of n
          neighbour->ComputeNormal(vn0,vn1,vn2,n_n);
          // the cosine is n_f * n_n
          cs = double(vtkMath::Dot(n_f,n_n));
          // the sin is (n_f x n_n) * e
          vtkMath::Cross(n_f,n_n,t);
          sn = double(vtkMath::Dot(t,e));
          // signed angle in [-pi,pi]
          if (sn!=0.0 || cs!=0.0)
            {
            angle = atan2(sn,cs);
            Hf    = length*angle;
            }
          else
            {
            Hf = 0.0;
            }
          // add weighted Hf to scalar at v_l and v_r
          if (Af!=0.0)
            {
            (Hf /= Af) *=3.0;
            }
          meanCurvatureData[v_l] += Hf;
          meanCurvatureData[v_r] += Hf;
          num_neighb[v_l] += 1;
          num_neighb[v_r] += 1;
          }
        }
      }

    // put curvature in vtkArray
    for (v = 0; v < numPts; v++)
      {
        if (num_neighb[v]>0)
          {
          Hf = 0.5*meanCurvatureData[v]/(double)num_neighb[v];
          if (this->InvertMeanCurvature)
            {
            meanCurvatureData[v] = -Hf;
            }
          else
            {
            meanCurvatureData[v] = Hf;
            }
          }
        else
          {
          meanCurvatureData[v] = 0.0;
          }
      }

    mesh->GetPointData()->AddArray(meanCurvature);
    mesh->GetPointData()->SetActiveScalars("Mean_Curvature");

    vtkDebugMacro("Set Values of Mean Curvature: Done");
    // clean
    vertices  ->Delete();
    vertices_n->Delete();
    neighbours->Delete();
    facet     ->Delete();
    neighbour ->Delete();

    if (meanCurvature) meanCurvature->Delete();
    if (num_neighb) delete [] num_neighb;

    cout << "out mean curv" << endl;
};
//--------------------------------------------
#define CLAMP_MACRO(v)    ((v)<(-1) ? (-1) : (v) > (1) ? (1) : v)
void vtkCurvatures1::GetGaussCurvature(vtkPolyData *output)
{
    cout << "woah" << endl;
    vtkDebugMacro("Start vtkCurvatures1::GetGaussCurvature()");
    // vtk data
    vtkCellArray* facets = output->GetPolys();

    // Empty array check
    if (output->GetNumberOfPolys()==0 || output->GetNumberOfPoints()==0)
      {
      vtkErrorMacro("No points/cells to operate on");
      return;
      }

    vtkTriangle* facet = vtkTriangle::New();

    // other data
    vtkIdType Nv   = output->GetNumberOfPoints();

    double* K = new double[Nv];
    double* dA = new double[Nv];
    double pi2 = 2.0*vtkMath::Pi();
    for (int k = 0; k < Nv; k++)
      {
      K[k]  = pi2;
      dA[k] = 0.0;
      }

    double v0[3], v1[3], v2[3], e0[3], e1[3], e2[3];

    double A, alpha0, alpha1, alpha2;

    vtkIdType f, *vert=0;
    facets->InitTraversal();
    while (facets->GetNextCell(f,vert))
      {
      output->GetPoint(vert[0],v0);
      output->GetPoint(vert[1],v1);
      output->GetPoint(vert[2],v2);
      // edges
      e0[0] = v1[0] ; e0[1] = v1[1] ; e0[2] = v1[2] ;
      e0[0] -= v0[0]; e0[1] -= v0[1]; e0[2] -= v0[2];

      e1[0] = v2[0] ; e1[1] = v2[1] ; e1[2] = v2[2] ;
      e1[0] -= v1[0]; e1[1] -= v1[1]; e1[2] -= v1[2];

      e2[0] = v0[0] ; e2[1] = v0[1] ; e2[2] = v0[2] ;
      e2[0] -= v2[0]; e2[1] -= v2[1]; e2[2] -= v2[2];

      // normalise
      vtkMath::Normalize(e0); vtkMath::Normalize(e1); vtkMath::Normalize(e2);
      // angles
      // I get lots of acos domain errors so clamp the value to +/-1 as the
      // normalize function can return 1.000000001 etc (I think)
      double ac1 = vtkMath::Dot(e1,e2);
      double ac2 = vtkMath::Dot(e2,e0);
      double ac3 = vtkMath::Dot(e0,e1);
      alpha0 = acos(-CLAMP_MACRO(ac1));
      alpha1 = acos(-CLAMP_MACRO(ac2));
      alpha2 = acos(-CLAMP_MACRO(ac3));

      // surf. area
      A = double(facet->TriangleArea(v0,v1,v2));
      // UPDATE
      dA[vert[0]] += A;
      dA[vert[1]] += A;
      dA[vert[2]] += A;
      K[vert[0]] -= alpha1;
      K[vert[1]] -= alpha2;
      K[vert[2]] -= alpha0;
      }

    int numPts = output->GetNumberOfPoints();
    // put curvature in vtkArray
    vtkDoubleArray* gaussCurvature = vtkDoubleArray::New();
    gaussCurvature->SetName("Gauss_Curvature");
    gaussCurvature->SetNumberOfComponents(1);
    gaussCurvature->SetNumberOfTuples(numPts);
    double *gaussCurvatureData = gaussCurvature->GetPointer(0);

    for (int v = 0; v < Nv; v++)
      {
      if (dA[v]>0.0)
        {
        gaussCurvatureData[v] = 3.0*K[v]/dA[v];
        }
      else
        {
        gaussCurvatureData[v] = 0.0;
        }
      }

    output->GetPointData()->AddArray(gaussCurvature);
    output->GetPointData()->SetActiveScalars("Gauss_Curvature");

    vtkDebugMacro("Set Values of Gauss Curvature: Done");
    /*******************************************************/
    if (facet) facet->Delete();
    if (K)              delete [] K;
    if (dA)             delete [] dA;
    if (gaussCurvature) gaussCurvature->Delete();
    /*******************************************************/
};

void vtkCurvatures1::GetMaximumCurvature(vtkPolyData *input,vtkPolyData *output)
{
  cout << "in max curv" << endl;
  this->GetGaussCurvature(output);
  this->GetMeanCurvature(output);

  vtkIdType numPts = input->GetNumberOfPoints();

  vtkDoubleArray *maximumCurvature = vtkDoubleArray::New();
  maximumCurvature->SetNumberOfComponents(1);
  maximumCurvature->SetNumberOfTuples(numPts);
  maximumCurvature->SetName("Maximum_Curvature");
  output->GetPointData()->AddArray(maximumCurvature);
  output->GetPointData()->SetActiveScalars("Maximum_Curvature");
  maximumCurvature->Delete();

  vtkDoubleArray *gauss = (vtkDoubleArray *)output->GetPointData()->GetArray("Gauss_Curvature");
  vtkDoubleArray *mean = (vtkDoubleArray *)output->GetPointData()->GetArray("Mean_Curvature");
  double k, h, k_max,tmp;

  for (vtkIdType i = 0; i<numPts; i++)
    {
    k = gauss->GetComponent(i,0);
    h = mean->GetComponent(i,0);
    tmp = h*h - k;
    if (tmp >= 0)
      {
      k_max = h + sqrt(tmp);
      }
    else
      {
      vtkDebugMacro(<< "Maximum Curvature undefined at point: " << i);
      // k_max can be any real number. Undefined points will be indistinguishable
      // from points that actually have a k_max == 0
      k_max = 0;
      }
    maximumCurvature->SetComponent(i, 0, k_max);
    }

    cout << "out max curv" << endl;
}

void vtkCurvatures1::GetMinimumCurvature(vtkPolyData *input,vtkPolyData *output)
{
  cout << "segfault1" << endl;
  this->GetGaussCurvature(output);
  this->GetMeanCurvature(output);

  vtkIdType numPts = input->GetNumberOfPoints();

  vtkDoubleArray *minimumCurvature = vtkDoubleArray::New();
  minimumCurvature->SetNumberOfComponents(1);
  minimumCurvature->SetNumberOfTuples(numPts);
  minimumCurvature->SetName("Minimum_Curvature");
  output->GetPointData()->AddArray(minimumCurvature);
  output->GetPointData()->SetActiveScalars("Minimum_Curvature");
  minimumCurvature->Delete();

  cout << "segfault2" << endl;

  vtkDoubleArray *gauss = (vtkDoubleArray *)output->GetPointData()->GetArray("Gauss_Curvature");
  vtkDoubleArray *mean = (vtkDoubleArray *)output->GetPointData()->GetArray("Mean_Curvature");
  double k, h, k_max,tmp;

  for (vtkIdType i = 0; i<numPts; i++)
    {
    k = gauss->GetComponent(i,0);
    h = mean->GetComponent(i,0);
    tmp = h*h - k;
    if (tmp >= 0)
      {
      k_max = h + sqrt(tmp);
      }
    else
      {
      vtkDebugMacro(<< "Minimum Curvature undefined at point: " << i);
      // k_max can be any real number. Undefined points will be indistinguishable
      // from points that actually have a k_max == 0
      k_max = 0;
      }
    minimumCurvature->SetComponent(i, 0, k_max);
    }

    cout << "segfault3" << endl;
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
  if (!input)
    {
    return 0;
    }

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

  if ( this->CurvatureType == VTK_CURVATURE_GAUSS )
    {
      this->GetGaussCurvature(output);
    }
  else if ( this->CurvatureType == VTK_CURVATURE_MEAN )
    {
      this->GetMeanCurvature(output);
    }
  else if ( this->CurvatureType ==  VTK_CURVATURE_MAXIMUM )
    {
      this->GetMaximumCurvature(input, output);
    }
  else if ( this->CurvatureType ==  VTK_CURVATURE_MINIMUM )
    {
      this->GetMinimumCurvature(input, output);
    }
  else
    {
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
