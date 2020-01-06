/*=========================================================================
  Program:   Visualization Toolkit
  Module:  vtkDijkstraGraphGeodesicPath1.cxx
  Language:  C++
  Made by Rasmus Paulsen
  email:  rrp(at)imm.dtu.dk
  web:    www.imm.dtu.dk/~rrp/VTK
  This class is not mature enough to enter the official VTK release.
=========================================================================*/
#include "vtkDijkstraGraphGeodesicPath1.h"

#include "vtkCellArray.h"
#include "vtkDijkstraGraphInternals.h"
#include "vtkExecutive.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkIdList.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkCurvatures.h"

#define ARRAY_SIZE(array) (sizeof((array))/sizeof((array[0])))
vtkStandardNewMacro(vtkDijkstraGraphGeodesicPath1);
vtkCxxSetObjectMacro(vtkDijkstraGraphGeodesicPath1,RepelVertices,vtkPoints);

//----------------------------------------------------------------------------
vtkDijkstraGraphGeodesicPath1::vtkDijkstraGraphGeodesicPath1()
{
  this->IdList = vtkIdList::New();
  this->Curvature = vtkDoubleArray::New();
  //this->Curvature = new double[1000];
  this->Internals = new vtkDijkstraGraphInternals;
  this->StopWhenEndReached = 0;
  this->UseScalarWeights = 0;
  this->NumberOfVertices = 0;
  this->RepelPathFromVertices = 0;
  this->RepelVertices = nullptr;
  this->maxCurv = -1E20;
  this->minCurv = 1E20;
  this->LineType = 0;
}

//----------------------------------------------------------------------------
vtkDijkstraGraphGeodesicPath1::~vtkDijkstraGraphGeodesicPath1()
{
  if (this->IdList)
  {
    this->IdList->Delete();
  }
  if (this->Curvature)
  {
    //delete[] this->Curvature;
    this->Curvature->Delete();
  }
  delete this->Internals;
  this->SetRepelVertices(nullptr);
}

//----------------------------------------------------------------------------
void vtkDijkstraGraphGeodesicPath1::GetCumulativeWeights(vtkDoubleArray *weights)
{
  if (!weights)
  {
    return;
  }

  weights->Initialize();
  double *weightsArray = new double[this->Internals->CumulativeWeights.size()];
  std::copy(this->Internals->CumulativeWeights.begin(),
    this->Internals->CumulativeWeights.end(), weightsArray);
  weights->SetArray(weightsArray, static_cast<vtkIdType>(this->Internals->CumulativeWeights.size()), 0);
}

//----------------------------------------------------------------------------
int vtkDijkstraGraphGeodesicPath1::RequestData(
  vtkInformation *           vtkNotUsed( request ),
  vtkInformationVector **    inputVector,
  vtkInformationVector *     outputVector)
{
  vtkInformation * inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo =   outputVector->GetInformationObject(0);

  vtkPolyData *input = vtkPolyData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
    cout << "getting input" << endl;
  if (!input)
  {
    return 0;
  }

  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));
    cout << "getting output" << endl;
  if (!output)
  {
    return 0;
  }
  this->NumberOfVertices = input->GetNumberOfPoints();
  this->Curvature->SetNumberOfValues(this->NumberOfVertices);
  cout << "generating curvature" << endl;
  this->GenCurvature(input);
  cout << "calculating min and max curvature" << endl;
  this->CalcMinMaxCurv();

//what is this supposed to do lol
  if ( this->AdjacencyBuildTime.GetMTime() < input->GetMTime() )
  {
    cout << "about to initialize..." << endl;
    this->Initialize( input );
    cout << "initialized" << endl;
  }
  else
  {
    //this->Reset();
    cout << "about to initialize..." << endl;
    this->Initialize( input );
    cout << "initialized" << endl;

  }

  if (this->NumberOfVertices == 0)
  {
    return 0;
  }

  this->ShortestPath( input, this->StartVertex, this->EndVertex );
  this->TraceShortestPath( input, output, this->StartVertex, this->EndVertex );
  return 1;
}

void vtkDijkstraGraphGeodesicPath1::GenCurvature(vtkPolyData *in) {
  vtkSmartPointer<vtkCurvatures> curv =
    vtkSmartPointer<vtkCurvatures>::New();
  curv->SetInputData(in);
  curv->SetCurvatureTypeToMaximum();
  curv->Update();

  cout << "calculating curvature" << endl;

  vtkSmartPointer<vtkPolyData> curvOutput =
    vtkSmartPointer<vtkPolyData>::New();
  curvOutput->ShallowCopy(curv->GetOutput());

  //temp->SetNumberOfValues(this->NumberOfVertices);
  vtkDoubleArray* temp = vtkDoubleArray::SafeDownCast(curvOutput->GetPointData()->GetArray("Maximum_Curvature"));
  //this->Curvature = curvOutput->GetPointData()->GetArray("Maximum_Curvature");
  //this->Curvature = static_cast<vtkDoubleArray *>(
    //curvOutput->GetPointData()->GetArray("Mean_Curvature"));
  //for debugging
  for (int i = 0; i < this->NumberOfVertices; i++) {
    this->Curvature->SetValue(i, temp->GetValue(i));
  }
  //temp->Delete();

}

void vtkDijkstraGraphGeodesicPath1::CalcMinMaxCurv() {
  double temp = 0.0;
  for (int i = 0; i < this->NumberOfVertices; i++) {
    temp = this->Curvature->GetValue(i);
    if (temp > this->maxCurv) {
      this->maxCurv = temp;
    }
    if (temp < this->minCurv) {
      this->minCurv = temp;
    }
  }
  cout << "max curv:" << this->maxCurv << endl;
  cout << "min curv:" << this->minCurv << endl;

  cout << "calculated max and min curvature values" << endl;
}

//----------------------------------------------------------------------------
void vtkDijkstraGraphGeodesicPath1::Initialize( vtkDataSet *inData )
{

  this->Internals->CumulativeWeights.resize( this->NumberOfVertices );
  this->Internals->Predecessors.resize( this->NumberOfVertices );
  this->Internals->OpenVertices.resize( this->NumberOfVertices );
  this->Internals->ClosedVertices.resize( this->NumberOfVertices );
  this->Internals->Adjacency.clear( );
  this->Internals->Adjacency.resize( this->NumberOfVertices );
  this->Internals->BlockedVertices.resize( this->NumberOfVertices );
  // The heap has elements from 1 to n
  this->Internals->InitializeHeap( this->NumberOfVertices );

  this->Reset();
  this->BuildAdjacency( inData );
}

//----------------------------------------------------------------------------
void vtkDijkstraGraphGeodesicPath1::Reset()
{
  std::fill( this->Internals->CumulativeWeights.begin(),
    this->Internals->CumulativeWeights.end(), -1.0 );
  std::fill( this->Internals->Predecessors.begin(),
    this->Internals->Predecessors.end(), -1 );
  std::fill( this->Internals->OpenVertices.begin(),
    this->Internals->OpenVertices.end(), false );
  std::fill( this->Internals->ClosedVertices.begin(),
    this->Internals->ClosedVertices.end(), false );
  if( this->RepelPathFromVertices )
  {
    std::fill( this->Internals->BlockedVertices.begin(),
      this->Internals->BlockedVertices.end(), false );
  }

  this->IdList->Reset();
  this->Internals->ResetHeap();
}

//----------------------------------------------------------------------------
double vtkDijkstraGraphGeodesicPath1::CalculateStaticEdgeCost(
     vtkDataSet *inData, vtkIdType u, vtkIdType v)
{
  double p1[3];
  inData->GetPoint(u,p1);
  double p2[3];
  inData->GetPoint(v,p2);

  //assuming this just calculates euclidean distance
  double dist = sqrt(vtkMath::Distance2BetweenPoints(p1, p2));

  if (this->UseScalarWeights)
  {
    // Note this edge cost is not symmetric!
    vtkFloatArray *scalars =
      static_cast<vtkFloatArray*>(inData->GetPointData()->GetScalars());
    //    float s1 = scalars->GetValue(u);
    double s2 = static_cast<double>(scalars->GetValue(v));

    double wt = s2*s2;
    if (wt != 0.0)
    {
      dist  /= wt;
    }
  }
  return dist;
}

//----------------------------------------------------------------------------
// This is probably a horribly inefficient way to do it.
void vtkDijkstraGraphGeodesicPath1::BuildAdjacency(vtkDataSet *inData)
{

  vtkPolyData *pd = vtkPolyData::SafeDownCast( inData );
  vtkIdType ncells = pd->GetNumberOfCells();
  cout << "starting to build adjacency matrix..." << endl;
  for ( vtkIdType i = 0; i < ncells; i++)
  {
    // Possible types
    //    VTK_VERTEX, VTK_POLY_VERTEX, VTK_LINE,
    //    VTK_POLY_LINE,VTK_TRIANGLE, VTK_QUAD,
    //    VTK_POLYGON, or VTK_TRIANGLE_STRIP.

    vtkIdType ctype = pd->GetCellType(i);

    // Until now only handle polys and triangles
    // TODO: All types
    if (ctype == VTK_POLYGON || ctype == VTK_TRIANGLE || ctype == VTK_LINE)
    {
      vtkIdType *pts;
      vtkIdType npts;
      pd->GetCellPoints(i, npts, pts);
      double cost;

      for (int j = 0; j < npts; ++j)
      {
        vtkIdType u = pts[j];
        vtkIdType v = pts[(( j + 1 ) % npts)];
        double curvu, curvv, ff, dist;

        std::map<int,double>& mu = this->Internals->Adjacency[u];
        if ( mu.find(v) == mu.end() )
        {
          //might not be symmetric
          dist = this->CalculateStaticEdgeCost( inData, u, v );
          //is symmetric
          cost = CalculateCostWithCurv(u, v, dist);
          mu.insert( std::pair<int,double>( v, cost ) );
        }

        std::map<int,double>& mv = this->Internals->Adjacency[v];
        if ( mv.find(u) == mv.end() )
        {
          dist = this->CalculateStaticEdgeCost( inData, v, u );
          cost = CalculateCostWithCurv(v, u, dist);
          mv.insert( std::pair<int,double>( u, cost ) );
        }
      }
    }
  }

  this->AdjacencyBuildTime.Modified();
}

double vtkDijkstraGraphGeodesicPath1::CalculateCostWithCurv(vtkIdType i, vtkIdType j, double dist) {
  double curvi = this->Curvature->GetValue(i);
  double curvj = this->Curvature->GetValue(j);
  double ff = 0.0;
  double cost = 0.0;
  switch(LineType) {
    //sulci
    case 1:
      ff = (this->minCurv - 0.5 * (curvi + curvj))/0.1;
      cost = dist * ff * ff;
      //cout << "sulci" << endl;
      break;
    //gyri
    case 2:
      ff = (this->maxCurv - 0.5 * (curvi + curvj))/0.1;
      cost = dist * ff * ff;
      //cout << "gyri" << endl;
      break;
    //geodesic
    case 0:
      cost = dist;
      //cout << "geodesic" << endl;
      break;
  }
  return cost;
}

//----------------------------------------------------------------------------
void vtkDijkstraGraphGeodesicPath1::TraceShortestPath(
               vtkDataSet *inData, vtkPolyData *outPoly,
               vtkIdType startv, vtkIdType endv)
{
  vtkPoints   *points = vtkPoints::New();
  vtkCellArray *lines = vtkCellArray::New();

  // n is far to many. Adjusted later
  lines->InsertNextCell(this->NumberOfVertices);

  // trace backward
  vtkIdType v = endv;
  double pt[3];
  vtkIdType id;
  while (v != startv)
  {
    IdList->InsertNextId(v);

    inData->GetPoint(v,pt);
    id = points->InsertNextPoint(pt);
    lines->InsertCellPoint(id);

    v = this->Internals->Predecessors[v];
  }

  this->IdList->InsertNextId(v);

  inData->GetPoint(v,pt);
  id = points->InsertNextPoint(pt);
  lines->InsertCellPoint(id);

  lines->UpdateCellCount( points->GetNumberOfPoints() );
  outPoly->SetPoints(points);
  points->Delete();
  outPoly->SetLines(lines);
  lines->Delete();
}

//----------------------------------------------------------------------------
void vtkDijkstraGraphGeodesicPath1::Relax(const int& u, const int& v, const double& w)
{
  double du = this->Internals->CumulativeWeights[u] + w;
  if (this->Internals->CumulativeWeights[v] > du)
  {
    this->Internals->CumulativeWeights[v] = du;
    this->Internals->Predecessors[v] = u;

    this->Internals->HeapDecreaseKey(v);
  }
}

//----------------------------------------------------------------------------
void vtkDijkstraGraphGeodesicPath1::ShortestPath( vtkDataSet *inData,
                                                int startv, int endv )
{
  int u, v;

  if( this->RepelPathFromVertices && this->RepelVertices )
  {
    // loop over the pts and if they are in the image
    // get the associated index for that point and mark it as blocked
    for( int i = 0; i < this->RepelVertices->GetNumberOfPoints(); ++i )
    {
        double* pt = this->RepelVertices->GetPoint( i );
        u = inData->FindPoint( pt );
        if ( u < 0 || u == startv || u == endv )
        {
          continue;
        }
        this->Internals->BlockedVertices[u] = true;
    }
  }

  this->Internals->CumulativeWeights[startv] = 0;

  this->Internals->HeapInsert(startv);
  this->Internals->OpenVertices[startv] = true;

  bool stop = false;
  while ((u = this->Internals->HeapExtractMin()) >= 0 && !stop)
  {
    // u is now in ClosedVertices since the shortest path to u is determined
    this->Internals->ClosedVertices[u] = true;
    // remove u from OpenVertices
    this->Internals->OpenVertices[u] = false;

    if (u == endv && this->StopWhenEndReached)
    {
      stop = true;
    }

    std::map<int,double>::iterator it = this->Internals->Adjacency[u].begin();

    // Update all vertices v adjacent to u
    for ( ; it != this->Internals->Adjacency[u].end(); ++it )
    {
      v = (*it).first;

      // ClosedVertices is the set of vertices with determined shortest path...
      // do not use them again
      if ( !this->Internals->ClosedVertices[v] )
      {
        // Only relax edges where the end is not in ClosedVertices
        // and edge is in OpenVertices
        double w;
        if ( this->Internals->BlockedVertices[v] )
        {
          w = VTK_FLOAT_MAX;
        }
        else
        {
          w = (*it).second + this->CalculateDynamicEdgeCost( inData, u, v );
        }

        if ( this->Internals->OpenVertices[v] )
        {
          this->Relax(u, v, w);
        }
        // add edge v to OpenVertices
        else
        {
          this->Internals->OpenVertices[v] = true;
          //set new cumulative weight of v to be the weight of start vertex + cost between them
          this->Internals->CumulativeWeights[v] = this->Internals->CumulativeWeights[u] + w;

          // Set Predecessor of v to be u
          this->Internals->Predecessors[v] = u;
          this->Internals->HeapInsert(v);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------
void vtkDijkstraGraphGeodesicPath1::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "StopWhenEndReached: ";
  if ( this->StopWhenEndReached )
  {
    os << "On\n";
  }
  else
  {
    os << "Off\n";
  }
  os << indent << "UseScalarWeights: ";
  if ( this->UseScalarWeights )
  {
    os << "On\n";
  }
  else
  {
    os << "Off\n";
  }
  os << indent << "RepelPathFromVertices: ";
  if ( this->RepelPathFromVertices )
  {
    os << "On\n";
  }
  else
  {
    os << "Off\n";
  }
  os << indent << "RepelVertices: " << this->RepelVertices << endl;
  os << indent << "IdList: " << this->IdList << endl;
  os << indent << "Number of vertices in input data: " << this->NumberOfVertices << endl;
}
