/*=========================================================================
  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkCurvatures1.h,v $
  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.
     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#ifndef __vtkCurvatures1_h
#define __vtkCurvatures1_h

#include "vtkPolyDataAlgorithm.h"
#include "vtkIdList.h"
#include <vector>
#include <tuple>
#include <utility>
#include <valarray>
// #include <ADL/Lapack.h>

using std::vector;
using std::tuple;
using std::valarray;
using std::pair;

#define TOLERANCE 0.0001                   // for curvature calculations
#define MAX_CURV 20
#define VTK_CURVATURE_GAUSS 0
#define VTK_CURVATURE_MEAN  1
#define VTK_CURVATURE_MAXIMUM 2
#define VTK_CURVATURE_MINIMUM 3

class VTK_EXPORT vtkCurvatures1 : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(vtkCurvatures1,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Construct with curvature type set to Gauss
  static vtkCurvatures1 *New();

  static void getPlane(double&, double&, double& , double&, const valarray<double>, const valarray<double>);
  static void getBasisVectors(valarray<double>&, valarray<double>&, valarray<double>&, const valarray<double>&);
  static double checkCurv(double);
  static bool addNeighbor(vtkSmartPointer<vtkIdList>&, vtkIdType, bool);

  // static void myCross(const double a[3], const double b[3], double c[3]);
  static void myCross(valarray<double> a, valarray<double> b, valarray<double>& c);
  static double myNorm(valarray<double> temp);

  // Description:
  // Set/Get Curvature type
  // VTK_CURVATURE_GAUSS: Gaussian curvature, stored as
  // DataArray "Gauss_Curvature"
  // VTK_CURVATURE_MEAN : Mean curvature, stored as
  // DataArray "Mean_Curvature"
  vtkSetMacro(CurvatureType,int);
  vtkGetMacro(CurvatureType,int);
  void SetCurvatureTypeToGaussian()
  { this->SetCurvatureType(VTK_CURVATURE_GAUSS); }
  void SetCurvatureTypeToMean()
  { this->SetCurvatureType(VTK_CURVATURE_MEAN); }
  void SetCurvatureTypeToMaximum()
  { this->SetCurvatureType(VTK_CURVATURE_MAXIMUM); }
  void SetCurvatureTypeToMinimum()
  { this->SetCurvatureType(VTK_CURVATURE_MINIMUM); }

  // Description:
  // Set/Get the flag which inverts the mean curvature calculation for
  // meshes with inward pointing normals (default false)
  vtkSetMacro(InvertMeanCurvature,int);
  vtkGetMacro(InvertMeanCurvature,int);
  vtkBooleanMacro(InvertMeanCurvature,int);

  vtkSetMacro(ndepth, int);
  vtkGetMacro(ndepth, int);

  vtkSetMacro(dx, double);
  vtkGetMacro(dx, double);

  vtkSetMacro(dy, double);
  vtkGetMacro(dy, double);

  vtkSetMacro(dz, double);
  vtkGetMacro(dz, double);

protected:
  vtkCurvatures1();
  ~vtkCurvatures1() override;

  // Usual data generation method
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);


  // Description: Principal Curvatures
  // Principal curvatures obtained by extracting eigenvalues from shape operator
  // lambda_1 and lambda_2
  void GetPrincipalCurvature(vtkPolyData* mesh, int ndepth, int dx, int dy, int dz);

  // Description: Gauss Curvature
  // Determinant of the shape operator or lambda_1 * lambda_2
  void GetGaussCurvature(vtkPolyData *);

  // Description: Mean Curvature
  // Average of eigenvalues of shape operator (lambda_1 and lambda_2)
  void GetMeanCurvature(vtkPolyData *);

  // Description: Max Curvature
  // Max of principal curvatures
  void GetMaximumCurvature(vtkPolyData *);

  // Description: Min Curvature
  // min of principal curvatures
  void GetMinimumCurvature(vtkPolyData *);

  void genNeighborhoods(vtkDataSet*, int);

  int getMaxNeighbors();

  bool hasNormals();

  bool hasUnitNormals();

  void genUnitNormals(vtkPolyData* mesh);

  void genNormals(vtkPolyData* mesh);

  void copyNeighbors(vector<vtkSmartPointer<vtkIdList>> orig,
      vector<vtkSmartPointer<vtkIdList>>& copy);
  // Vars
  int CurvatureType;
  int InvertMeanCurvature;
  int ndepth;
  double dx, dy, dz;
  int numPoints;
  int numPolys;

  // vtkDoubleArray* prinCurvature;
  vector<pair<double, double>> prinCurvature;

  vector<vtkSmartPointer<vtkIdList>> neighbors;

  // vector<tuple<double, double, double>> normals;
  vector<valarray<double>> normals;
  // vector<tuple<double, double, double>> unitNormals;
  vector<valarray<double>> unitNormals;

  // vtkPolyData* output_mesh;

private:
  vtkCurvatures1(const vtkCurvatures1&);  // Not implemented.
  void operator=(const vtkCurvatures1&);  // Not implemented.

};

#endif
