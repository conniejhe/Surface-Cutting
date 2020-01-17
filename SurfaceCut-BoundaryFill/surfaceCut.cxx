#include "surfaceCut.h"
#include "vtkIdList.h"
#include "vtkIdTypeArray.h"
#include "vtkIntArray.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkDijkstraGraphGeodesicPath.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkImplicitSelectionLoop.h"
#include "vtkClipPolyData.h"
#include "vtkSelectPolyData.h"
#include "vtkCleanPolyData.h"
#include "vtkAppendPolyData.h"
#include <queue>
#include <unordered_map>
#include <vector>
#include <iostream>

using std::vector;
using std::cout;
using std::queue;
using std::unordered_map;

vtkStandardNewMacro(surfaceCut);

//constructor
surfaceCut::surfaceCut() {
    this->SetNumberOfInputPorts(3);
    this->UserPoints = vtkIdList::New();
    this->colorArray = vtkIntArray::New();
}

surfaceCut::~surfaceCut() {
    if (this->UserPoints) {
        this->UserPoints->Delete();
    }
    if (this->colorArray) {
        this->colorArray->Delete();
    }
}

int surfaceCut::FillInputPortInformation(int port, vtkInformation* info)
{
    if (port == 0) {
        info->Remove(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE());
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
        return 1;
    } else if (port == 1) {
        info->Remove(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE());
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
        return 1;
    } else if (port == 2) {
        info->Remove(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE());
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
        return 1;
    }
    return 0;
}


int surfaceCut::RequestData(vtkInformation* vtkNotUsed(request),
    vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *lineInfo = inputVector[1]->GetInformationObject(0);
    vtkInformation *selectionInfo = inputVector[2]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    vtkPolyData *input = vtkPolyData::SafeDownCast(
        inInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkPolyData *line = vtkPolyData::SafeDownCast(
        lineInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkUnstructuredGrid *sel = vtkUnstructuredGrid::SafeDownCast(
        selectionInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkPolyData *output = vtkPolyData::SafeDownCast(
        outInfo->Get(vtkDataObject::DATA_OBJECT()));

    cout << "Obtained Data Objects." << endl;

    if (!output || !input || !line) {
    return 0;
    }

    this->NumberOfVertices = input->GetNumberOfPoints();

    /* Do we need to do this? */
    for (int i = 0; i < this->NumberOfVertices; i++) {
        this->adjacencyMatrix.push_back(vtkSmartPointer<vtkIdList>::New());
    }

    vtkIdTypeArray* selIds = vtkIdTypeArray::SafeDownCast(
        sel->GetPointData()->GetArray("vtkOriginalPointIds"));

    this->insidePoint = selIds->GetValue(0);

    vtkIdTypeArray* origIds = vtkIdTypeArray::SafeDownCast(
        line->GetPointData()->GetArray("vtkOriginalPointIds"));

    cout << "Got Point Ids in loop." << endl;

    for (vtkIdType i = 0; i < origIds->GetNumberOfTuples(); i++) {
      this->UserPoints->InsertNextId(origIds->GetValue(i));
    }

    cout << "Extracted Point IDs from loop." << endl;

    vtkPolyData* temp = vtkPolyData::New();

    cout << "created a new polydata" << endl;

    temp->ShallowCopy(input);

    cout << "shallow copied input polydata object" << endl;



    cout << "built Adjacency" << endl;

    // splitComponents(temp);

    ColorBoundary();

    BuildAdjacency(temp);

    FillBoundary(temp, this->insidePoint, 1, 2);

    CutSurface(temp, output);

    // vtkIdList* reachableNodes = findReachableNodes(temp);
    //
    // cut(temp, output, reachableNodes);

    return 1;
}

void surfaceCut::ColorBoundary() {
    // initialize color boundary array
    this->colorArray->SetNumberOfComponents(1);
    this->colorArray->SetNumberOfTuples(this->NumberOfVertices);

    for (int i = 0; i < this->UserPoints->GetNumberOfIds(); i++) {
        vtkIdType curr_id = this->UserPoints->GetId(i);
        this->colorArray->SetValue(curr_id, 1);
    }
}

void surfaceCut::FillBoundary(vtkDataSet *inData, vtkIdType i,
    int bound_color, int fill_color) {
    vtkPolyData *pd = vtkPolyData::SafeDownCast( inData );
    if (this->colorArray->GetValue(i) != bound_color &&
        this->colorArray->GetValue(i) != fill_color)
    {
        this->colorArray->SetValue(i, fill_color);
        int length = this->adjacencyMatrix[i]->GetNumberOfIds();

        for (vtkIdType n = 0; n < length; n++) {
            vtkIdType neighbor = this->adjacencyMatrix[i]->GetId(n);
            FillBoundary(inData, neighbor, bound_color, fill_color);
        }
    }
}

void surfaceCut::CutSurface(vtkPolyData* in, vtkPolyData* out) {
    cout << "starting to cut" << endl;

    // out->CopyStructure(in);
    // out->GetPointData()->PassData(in->GetPointData());
    // out->GetFieldData()->PassData(in->GetFieldData());

    vtkIdType F = in->GetNumberOfCells();

    in->BuildLinks();

    vtkIdType a, b, c;

    for (vtkIdType f = 0; f < F; f++) {
        cout << "on cell: " << f << endl;
        vtkIdType npts;
        vtkIdType* pts;

        in->GetCellPoints(f, npts, pts);

        a = pts[0];
        b = pts[1];
        c = pts[2];

        cout << "wtf" << endl;

        if (!(this->colorArray->GetValue(a)) || !(this->colorArray->GetValue(b))
            || !(this->colorArray->GetValue(c))) {
            cout << "deleting cell: " << f << endl;
            in->DeleteCell(f);
        }

        // set colorArray[a], [b], [c] to 1??
        // or have separate array that basically marks these cells for deletion
        // but doesn't actually delete until after in a separate loop?
    }

    in->RemoveDeletedCells();

    cout << "done with cut" << endl;

    vtkCleanPolyData *Clean = vtkCleanPolyData::New();
    Clean->SetInputData(in);
    Clean->Update();
    out->ShallowCopy(Clean->GetOutput());
}

// modeled off of buildAdjacency method in vtkDijkstraGraphGeodesicPath
void surfaceCut::BuildAdjacency(vtkDataSet *inData)
{

  vtkPolyData *pd = vtkPolyData::SafeDownCast( inData );
  vtkIdType ncells = pd->GetNumberOfCells();
  cout << "starting to build adjacency matrix..." << endl;
  for ( vtkIdType i = 0; i < ncells; i++)
  {

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

        vtkSmartPointer<vtkIdList> mu = this->adjacencyMatrix[u];
        if ( mu->IsId(v) == -1 )
        {
          mu->InsertNextId(v);
        }

        vtkSmartPointer<vtkIdList> mv = this->adjacencyMatrix[v];
        if ( mv->IsId(u) == -1 )
        {
          mv->InsertNextId(u);
        }
      }
    }
  }

  //this->AdjacencyBuildTime.Modified();???
}

// vtkIdList* surfaceCut::BFS(int compNum, vtkIdType src, vtkIntArray* visited) {
//     queue<vtkIdType> queue;
//
//     queue.push(src);
//
//     visited->SetValue(src, 1);
//
//     vtkIdList* reachableNodes = vtkIdList::New();
//
//     while(!queue.empty())
//     {
//         vtkIdType u = queue.front();
//         queue.pop();
//
//         reachableNodes->InsertNextId(u);
//
//         // get number of neighboring vertices
//         int length = this->adjacencyMatrix[u]->GetNumberOfIds();
//
//         // Get all adjacent vertices of the dequeued
//         // vertex u. If a adjacent has not been visited,
//         // then mark it visited and enqueue it
//         for (vtkIdType i = 0; i < length; i++) {
//
//             vtkIdType temp = this->adjacencyMatrix[u]->GetId(i);
//
//             if (!(visited->GetValue(temp))) {
//                 visited->SetValue(temp, 1);
//                 queue.push(temp);
//             }
//
//         }
//     }
//     return reachableNodes;
//
// }
//
// vtkIdList* surfaceCut::findReachableNodes(vtkPolyData* in) {
//
//     vtkSmartPointer<vtkIntArray> visited = vtkIntArray::New();
//     visited->SetNumberOfValues(this->NumberOfVertices + 1);
//
//     // Initialize all vertices with 0 in visited array (i.e. not visited yet)
//     for (int i = 0; i <= this->NumberOfVertices; i++) {
//         visited->SetValue(i, 0);
//     }
//
//     // vtkIdList to store reachable nodes
//     vtkIdList* reachedNodes;
//
//     int compNum = 0;
//
//     vtkIdType u = this->outsidePoint;
//
//     if(!(visited->GetValue(u))) {
//         compNum++;
//         reachedNodes = this->BFS(compNum, u, visited);
//         cout << "called BFS" << endl;
//     }
//     for (vtkIdType i = 0; i < reachedNodes->GetNumberOfIds(); i++) {
//         cout << "reachable node: " << reachedNodes->GetId(i) << endl;
//     }
//
//     return reachedNodes;
//
// }
//
// void surfaceCut::cut(vtkPolyData* in, vtkPolyData* out, vtkIdList* list) {
//     cout << "starting to cut" << endl;
//     // for (vtkIdType i = 0; i < list->GetNumberOfIds(); i++) {
//     //     in->BuildLinks();
//     //     cout << "deleted point: " << list->GetId(i) << endl;
//     //     in->DeletePoint(list->GetId(i));
//     //
//     // }
//
//     for (vtkIdType i = 0; i < list->GetNumberOfIds(); i++) {
//         in->BuildLinks();
//         // initialize list of cells
//         vtkSmartPointer<vtkIdList> cellIdList = vtkSmartPointer<vtkIdList>::New();
//         // get all cells that contain the current vertex
//         in->GetPointCells(list->GetId(i), cellIdList);
//         for (vtkIdType j = 0; j < cellIdList->GetNumberOfIds(); j++) {
//             // delete each of these cells
//             in->DeleteCell(cellIdList->GetId(j));
//         }
//         in->RemoveDeletedCells();
//     }
//     vtkCleanPolyData *Clean = vtkCleanPolyData::New();
//     Clean->SetInputData(in);
//     Clean->Update();
//     out->ShallowCopy(Clean->GetOutput());
// }
//
// /* splits poly data into two components by deleting all cells along the
//     user-specified path */
// void surfaceCut::splitComponents(vtkPolyData* in) {
//
//     // can delete this first for loop???
//     for (vtkIdType i = 0; i < this->UserPoints->GetNumberOfIds(); i++) {
//         vtkSmartPointer<vtkIdList> cellIdList = vtkSmartPointer<vtkIdList>::New();
//         // get all cells that contain the current vertex
//         in->GetPointCells(UserPoints->GetId(i), cellIdList);
//     }
//     // iterate through
//     for (vtkIdType i = 0; i < this->UserPoints->GetNumberOfIds(); i++) {
//         in->BuildLinks();
//         // initialize list of cells
//         vtkSmartPointer<vtkIdList> cellIdList = vtkSmartPointer<vtkIdList>::New();
//         // get all cells that contain the current vertex
//         in->GetPointCells(UserPoints->GetId(i), cellIdList);
//         for (vtkIdType j = 0; j < cellIdList->GetNumberOfIds(); j++) {
//             in->DeleteCell(cellIdList->GetId(j));
//         }
//         in->RemoveDeletedCells();
//     }
//
// }

void surfaceCut::PrintSelf(ostream& os, vtkIndent indent) {
    this->Superclass::PrintSelf(os, indent);

    os << indent << "User Selected Points: " << this->UserPoints << endl;
}
