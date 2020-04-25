/*=========================================================================
  Program:   Paraview
  Module:  vtkMyDGGP.cxx
  Language:  C++
  Author: Connie He
=========================================================================*/

#include "vtkMyDGGP.h"
#include "vtkIdList.h"
#include "vtkIdTypeArray.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "../../../SurfaceTracker-TextEntry/Plugin/FiltervtkDijkstraGraphGeodesicPath1.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkImplicitSelectionLoop.h"
#include "vtkClipPolyData.h"
#include "vtkSelectPolyData.h"
#include "vtkCleanPolyData.h"
#include "vtkAppendPolyData.h"
#include <vector>
#include <iostream>

using std::vector;
using std::cout;

vtkStandardNewMacro(vtkMyDGGP);

//constructor
vtkMyDGGP::vtkMyDGGP() {
  this->SetNumberOfInputPorts(2);
  this->UserPoints = vtkIdList::New();
  //vtkSmartPointer<vtkUserPoints> UserPoints = vtkSmartPointer<vtkUserPoints>::New();
  this->IdList = vtkIdList::New();
  this->LineType = VTK_LINE_GEODESIC;
}

vtkMyDGGP::~vtkMyDGGP() {
  if (this->UserPoints) {
      this->UserPoints->Delete();
  }
  if (this->IdList) {
      this->IdList->Delete();
  }
}

int vtkMyDGGP::FillInputPortInformation(int port, vtkInformation* info)
{
  if (port == 0) {
    info->Remove(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE());
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    return 1;
  }
  else if (port == 1) {
    info->Remove(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE());
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}


int vtkMyDGGP::RequestData(vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{

  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *selInfo = inputVector[1]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  vtkPolyData *input = vtkPolyData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkUnstructuredGrid *sel = vtkUnstructuredGrid::SafeDownCast(
    selInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  if (!output || !input || !sel) {
    return 0;
  }

  vtkSmartPointer<vtkIdTypeArray> origIds = vtkIdTypeArray::SafeDownCast(
    sel->GetPointData()->GetArray("vtkOriginalPointIds"));

    // put this in so that it resets every time you change line type
    if (this->UserPoints->GetNumberOfIds() > 0) {
        this->UserPoints->Reset();
        this->UserPoints->Squeeze();
    }

    if (origIds == NULL) {
        cerr << "ERROR: Not enough points selected." << endl;
        vtkErrorMacro(<< "ERROR: Not enough points selected.");
        return 0;
    }

  for (vtkIdType i = 0; i < origIds->GetNumberOfTuples(); i++) {
      this->UserPoints->InsertNextId(origIds->GetValue(i));
      cout << origIds->GetValue(i) << endl;
  }
  //appender to collect all geodesic paths
  vtkSmartPointer<vtkAppendPolyData> appender =
    vtkSmartPointer<vtkAppendPolyData>::New();
  int length = this->UserPoints->GetNumberOfIds();

  if (this->IdList->GetNumberOfIds() > 0) {
      this->IdList->Reset();
      this->IdList->Squeeze();
  }
  for (vtkIdType i = 0; i < (length - 1); i++) {
    //instantiate this outside of loop?
    vtkSmartPointer<vtkDijkstraGraphGeodesicPath1> dijkstra =
      vtkSmartPointer<vtkDijkstraGraphGeodesicPath1>::New();

    if (this->LineType == VTK_LINE_GEODESIC ) {
      dijkstra->SetLineTypeToGeodesic();
    } else if (this->LineType == VTK_LINE_SULCUS ) {
      dijkstra->SetLineTypeToSulcus();
    } else if (this->LineType == VTK_LINE_GYRUS) {
      dijkstra->SetLineTypeToGyrus();
    }
    dijkstra->SetInputData(input);
    dijkstra->SetStartVertex(this->UserPoints->GetId(i));
    dijkstra->SetEndVertex(this->UserPoints->GetId(i + 1));
    dijkstra->Update();
    vtkSmartPointer<vtkIdList> tempIdList = dijkstra->GetIdList();

    for (vtkIdType j = tempIdList->GetNumberOfIds() - 1; j >= 0; j--) {
      this->IdList->InsertUniqueId(tempIdList->GetId(j));
    }

    vtkSmartPointer<vtkPolyData> pathOutput =
      vtkSmartPointer<vtkPolyData>::New();
    pathOutput->ShallowCopy(dijkstra->GetOutput());
    appender->AddInputData(pathOutput);
  }

  vtkSmartPointer<vtkCleanPolyData> cleaner =
    vtkSmartPointer<vtkCleanPolyData>::New();
  cleaner->PointMergingOn();
  cleaner->SetInputConnection(appender->GetOutputPort());
  cleaner->Update();

  output->ShallowCopy(cleaner->GetOutput());

  vtkSmartPointer<vtkPointData> outputPD = output->GetPointData();

  vtkSmartPointer<vtkIdTypeArray> originalPointIds = vtkSmartPointer<vtkIdTypeArray>::New();

  originalPointIds->SetNumberOfComponents(1);
  originalPointIds->SetName("vtkOriginalPointIds");

  for (vtkIdType i = 0; i < this->IdList->GetNumberOfIds(); i++) {
    originalPointIds->InsertNextValue(this->IdList->GetId(i));
  }
  outputPD->AddArray(originalPointIds);

  return 1;
}

void vtkMyDGGP::PrintSelf(ostream& os, vtkIndent indent) {
  this->Superclass::PrintSelf(os, indent);

  os << indent << "User Selected Points: " << this->UserPoints << endl;
  os << indent << "IdList: " << this->IdList << endl;
}
