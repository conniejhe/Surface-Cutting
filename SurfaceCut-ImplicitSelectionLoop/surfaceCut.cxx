#include "surfaceCut.h"
#include "vtkIdList.h"
#include "vtkIdTypeArray.h"
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
#include <vector>
#include <iostream>

using std::vector;
using std::cout;

vtkStandardNewMacro(surfaceCut);

//constructor
surfaceCut::surfaceCut() {
  this->SetNumberOfInputPorts(2);
  this->UserPoints = vtkIdList::New();
}

surfaceCut::~surfaceCut() {
  if (this->UserPoints) {
      this->UserPoints->Delete();
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
  }
  return 0;
}


int surfaceCut::RequestData(vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{

  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *lineInfo = inputVector[1]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  vtkPolyData *input = vtkPolyData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkPolyData *line = vtkPolyData::SafeDownCast(
    lineInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  cout << "Obtained Data Objects." << endl;

  if (!output || !input || !line) {
    return 0;
  }

  vtkIdTypeArray* origIds = vtkIdTypeArray::SafeDownCast(
    line->GetPointData()->GetArray("vtkOriginalPointIds"));

  cout << "Got Point Ids in loop." << endl;

  for (vtkIdType i = 0; i < origIds->GetNumberOfTuples(); i++) {
      this->UserPoints->InsertNextId(origIds->GetValue(i));
  }

  cout << "Extracted Point IDs from loop." << endl;

  vtkSmartPointer<vtkPoints> selectionPoints =
    vtkSmartPointer<vtkPoints>::New();
  selectionPoints->SetNumberOfPoints(this->UserPoints->GetNumberOfIds());

  vtkSmartPointer<vtkPoints> temp = input->GetPoints();
  temp->GetPoints(this->UserPoints, selectionPoints);

  cout << "Extracted point coordinates from point IDs" << endl;

  vtkSmartPointer<vtkImplicitSelectionLoop> loop =
    vtkSmartPointer<vtkImplicitSelectionLoop>::New();
  loop->SetLoop(selectionPoints);

  vtkSmartPointer<vtkClipPolyData> clip =
    vtkSmartPointer<vtkClipPolyData>::New();
  clip->GenerateClippedOutputOn();
  clip->SetInputData(input);
  clip->SetClipFunction(loop);
  clip->SetValue(0.0);
  clip->Update();

  cout << "Done clipping." << endl;

  vtkSmartPointer<vtkPolyData> clipOutput =
    vtkSmartPointer<vtkPolyData>::New();
  clipOutput->ShallowCopy(clip->GetClippedOutput());

  output->ShallowCopy(clipOutput);

  return 1;
}

void surfaceCut::PrintSelf(ostream& os, vtkIndent indent) {
  this->Superclass::PrintSelf(os, indent);

  os << indent << "User Selected Points: " << this->UserPoints << endl;

}
