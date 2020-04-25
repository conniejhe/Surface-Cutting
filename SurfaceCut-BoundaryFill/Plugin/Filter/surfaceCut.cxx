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
#include "vtkCellData.h"
#include "vtkCleanPolyData.h"
#include "vtkUnstructuredGrid.h"

#include <vector>
#include <iostream>

using std::vector;
using std::cout;

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

    vtkPolyData* output = vtkPolyData::New();
    vtkPolyData *data_output = vtkPolyData::SafeDownCast(
        outInfo->Get(vtkDataObject::DATA_OBJECT()));

    if (!output || !input || !line) {
    return 0;
    }

    output->CopyStructure(input);
    output->GetPointData()->PassData(input->GetPointData());
    output->GetCellData()->PassData(input->GetCellData());
    output->GetFieldData()->PassData(input->GetFieldData());

    this->NumberOfVertices = input->GetNumberOfPoints();

    for (int i = 0; i < this->NumberOfVertices; i++) {
        this->adjacencyMatrix.push_back(vtkSmartPointer<vtkIdList>::New());
    }

    vtkIdTypeArray* selIds = vtkIdTypeArray::SafeDownCast(
        sel->GetPointData()->GetArray("vtkOriginalPointIds"));

    if (selIds == NULL) {
        cerr << "ERROR: Need to select at least one inside point." << endl;
        vtkErrorMacro(<< "ERROR: Need to select at least one inside point.");
        return 0;
    }

    this->insidePoint = selIds->GetValue(0);

    vtkIdTypeArray* origIds = vtkIdTypeArray::SafeDownCast(
        line->GetPointData()->GetArray("vtkOriginalPointIds"));

    if (origIds == NULL) {
        cerr << "ERROR: Loop does not contain any points." << endl;
        vtkErrorMacro(<< "ERROR: Loop does not contain any points.");
        return 0;
    }

    for (vtkIdType i = 0; i < origIds->GetNumberOfTuples(); i++) {
      this->UserPoints->InsertNextId(origIds->GetValue(i));
    }

    int numIds = this->UserPoints->GetNumberOfIds();

    if (this->UserPoints->GetId(0) != this->UserPoints->GetId(numIds - 1)) {
        vtkErrorMacro("ERROR: loop is not closed");
        return 0;
    }

    ColorBoundary();

    BuildAdjacency(output);

    FillBoundary(output, this->insidePoint, 1, 2);

    CutSurface(output,data_output);

    return 1;
}

void surfaceCut::ColorBoundary() {
    // initialize color boundary array
    this->colorArray->SetNumberOfComponents(1);
    this->colorArray->SetNumberOfTuples(this->NumberOfVertices);

    // initialize
    for (int j = 0; j < this->NumberOfVertices; j++) {
        this->colorArray->SetValue(j, 0);
    }

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

    vtkIdType F = in->GetNumberOfCells();

    in->BuildLinks();

    vtkIdType a, b, c;

    for (vtkIdType f = 0; f < F; f++) {
        vtkIdType npts;
        const vtkIdType* pts;

        in->GetCellPoints(f, npts, pts);
        a = pts[0];
        b = pts[1];
        c = pts[2];

        // keeps cell only if all three vertices are colored
        if (!(this->colorArray->GetValue(a)) || !(this->colorArray->GetValue(b))
            || !(this->colorArray->GetValue(c))) {

            vtkIdType not_colored = 0;
            if (!(this->colorArray->GetValue(a))) {
                not_colored = a;
            } else if (!(this->colorArray->GetValue(b))) {
                not_colored = b;
            } else {
                not_colored = c;
            }

            int num_nei = this->adjacencyMatrix[not_colored]->GetNumberOfIds();
            int colored_nei = 0;
            for (vtkIdType n = 0; n < num_nei; n++) {
                if(this->colorArray->GetValue(n)) {
                    colored_nei++;
                }
            }

            if (colored_nei < num_nei) {
                in->DeleteCell(f);
            }
        }
    }

    in->RemoveDeletedCells();

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
  for ( vtkIdType i = 0; i < ncells; i++)
  {

    vtkIdType ctype = pd->GetCellType(i);

    // Until now only handle polys and triangles
    // TODO: All types
    if (ctype == VTK_POLYGON || ctype == VTK_TRIANGLE || ctype == VTK_LINE)
    {
      const vtkIdType *pts;
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

}

void surfaceCut::PrintSelf(ostream& os, vtkIndent indent) {
    this->Superclass::PrintSelf(os, indent);

    os << indent << "User Selected Points: " << this->UserPoints << endl;
}
