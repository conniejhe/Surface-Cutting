/*=========================================================================
  Program:   Visualization Toolkit
  Module:    vtkDijkstraGraphGeodesicPath1.h
  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.
     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.
=========================================================================*/
/**
 * @class   vtkDijkstraGraphGeodesicPath1
 * @brief   Dijkstra algorithm to compute the graph geodesic.
 *
 * Takes as input a polygonal mesh and performs a single source shortest
 * path calculation. Dijkstra's algorithm is used. The implementation is
 * similar to the one described in Introduction to Algorithms (Second Edition)
 * by Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest, and
 * Cliff Stein, published by MIT Press and McGraw-Hill. Some minor
 * enhancement are added though. All vertices are not pushed on the heap
 * at start, instead a front set is maintained. The heap is implemented as
 * a binary heap. The output of the filter is a set of lines describing
 * the shortest path from StartVertex to EndVertex.
 *
 * @warning
 * The input polydata must have only triangle cells.
 *
 * @par Thanks:
 * The class was contributed by Rasmus Paulsen.
 * www.imm.dtu.dk/~rrp/VTK . Also thanks to Alexandre Gouaillard and Shoaib
 * Ghias for bug fixes and enhancements.
*/

#ifndef vtkDijkstraGraphGeodesicPath1_h
#define vtkDijkstraGraphGeodesicPath1_h

#include "vtkFiltersModelingModule.h" // for export macro
#include "vtkGraphGeodesicPath.h"

class vtkDijkstraGraphInternals;
class vtkIdList;

#define VTK_LINE_GEODESIC 0
#define VTK_LINE_SULCUS 1
#define VTK_LINE_GYRUS 2

#define VTK_CURVATURE_GAUSS 0
#define VTK_CURVATURE_MEAN  1
#define VTK_CURVATURE_MAXIMUM 2
#define VTK_CURVATURE_MINIMUM 3

class VTKFILTERSMODELING_EXPORT vtkDijkstraGraphGeodesicPath1 :
                           public vtkGraphGeodesicPath
{
public:

  /**
   * Instantiate the class
   */
  static vtkDijkstraGraphGeodesicPath1 *New();

  //@{
  /**
   * Standard methods for printing and determining type information.
   */
  vtkTypeMacro(vtkDijkstraGraphGeodesicPath1,vtkGraphGeodesicPath);
  void PrintSelf(ostream& os, vtkIndent indent) override;
  //@}

  vtkSetMacro(LineType, int);
  vtkGetMacro(LineType, int);
  void SetLineTypeToGeodesic() {
    this->SetLineType(VTK_LINE_GEODESIC);
  }
  void SetLineTypeToSulcus() {
    this->SetLineType(VTK_LINE_SULCUS);
  }
  void SetLineTypeToGyrus() {
    this->SetLineType(VTK_LINE_GYRUS);
  }

  //@{
  /**
   * Set/Get Curvature type
   * VTK_CURVATURE_GAUSS: Gaussian curvature, stored as
   * DataArray "Gauss_Curvature"
   * VTK_CURVATURE_MEAN : Mean curvature, stored as
   * DataArray "Mean_Curvature"
   */
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

  vtkSetMacro(ndepth, int);
  vtkGetMacro(ndepth, int);

  vtkSetMacro(dx, double);
  vtkGetMacro(dx, double);

  vtkSetMacro(dy, double);
  vtkGetMacro(dy, double);

  vtkSetMacro(dz, double);
  vtkGetMacro(dz, double);

  //@{
  /**
   * The vertex ids (of the input polydata) on the shortest path
   */
  vtkGetObjectMacro(IdList, vtkIdList);
  //@}

  //@{
  /**
   * Stop when the end vertex is reached
   * or calculate shortest path to all vertices
   */
  vtkSetMacro(StopWhenEndReached, vtkTypeBool);
  vtkGetMacro(StopWhenEndReached, vtkTypeBool);
  vtkBooleanMacro(StopWhenEndReached, vtkTypeBool);
  //@}

  //@{
  /**
   * Use scalar values in the edge weight (experimental)
   */
  vtkSetMacro(UseScalarWeights, vtkTypeBool);
  vtkGetMacro(UseScalarWeights, vtkTypeBool);
  vtkBooleanMacro(UseScalarWeights, vtkTypeBool);
  //@}

  //@{
  /**
   * Use the input point to repel the path by assigning high costs.
   */
  vtkSetMacro(RepelPathFromVertices, vtkTypeBool);
  vtkGetMacro(RepelPathFromVertices, vtkTypeBool);
  vtkBooleanMacro(RepelPathFromVertices, vtkTypeBool);
  //@}

  //@{
  /**
   * Specify vtkPoints to use to repel the path from.
   */
  virtual void SetRepelVertices(vtkPoints*);
  vtkGetObjectMacro(RepelVertices, vtkPoints);
  //@}

  /**
   * Fill the array with the cumulative weights.
   */
  virtual void GetCumulativeWeights(vtkDoubleArray *weights);

protected:
  vtkDijkstraGraphGeodesicPath1();
  ~vtkDijkstraGraphGeodesicPath1() override;

  int RequestData(vtkInformation *, vtkInformationVector **,
                          vtkInformationVector *) override;

  // Build a graph description of the input.
  virtual void BuildAdjacency( vtkDataSet *inData );

  vtkTimeStamp AdjacencyBuildTime;

  // The fixed cost going from vertex u to v.
  virtual double CalculateStaticEdgeCost( vtkDataSet *inData, vtkIdType u, vtkIdType v);

  //updated cost calculation to account for curvature
  virtual double CalculateCostWithCurv(vtkIdType i, vtkIdType j, double dist);

  // The cost going from vertex u to v that may depend on one or more vertices
  //that precede u.
  virtual double CalculateDynamicEdgeCost( vtkDataSet *, vtkIdType , vtkIdType )
  { return 0.0; }

  void GenCurvature(vtkPolyData *);

  vtkDoubleArray* GetCurvature();

  void CalcMinMaxCurv();

  void Initialize( vtkDataSet *inData );

  void Reset();

  // Calculate shortest path from vertex startv to vertex endv.
  virtual void ShortestPath( vtkDataSet *inData, int startv, int endv );

  // Relax edge u,v with weight w.
  void Relax(const int& u, const int& v, const double& w);

  // Backtrace the shortest path
  void TraceShortestPath( vtkDataSet* inData, vtkPolyData* outPoly,
               vtkIdType startv, vtkIdType endv);

  // The number of vertices.
  int NumberOfVertices;

  // The vertex ids on the shortest path.
  vtkIdList *IdList;

  //Internalized STL containers.
  vtkDijkstraGraphInternals *Internals;

  vtkTypeBool StopWhenEndReached;
  vtkTypeBool UseScalarWeights;
  vtkTypeBool RepelPathFromVertices;
  //vtkSmartPointer<vtkDoubleArray> Curvature;
  vtkDoubleArray* Curvature;
  //double* Curvature;
  vtkPoints* RepelVertices;

  double maxCurv;
  double minCurv;

  int LineType;
  int CurvatureType;
  int ndepth;
  double dx, dy, dz;

private:
  vtkDijkstraGraphGeodesicPath1(const vtkDijkstraGraphGeodesicPath1&) = delete;
  void operator=(const vtkDijkstraGraphGeodesicPath1&) = delete;


};

#endif
