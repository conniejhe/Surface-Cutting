#ifndef __vtkMyDGGP_h
#define __vtkMyDGGP_h

#include "vtkPolyDataAlgorithm.h"

class vtkIdList;

#define VTK_LINE_GEODESIC 0
#define VTK_LINE_SULCUS 1
#define VTK_LINE_GYRUS 2

class VTK_EXPORT vtkMyDGGP : public vtkPolyDataAlgorithm
{

public:
  static vtkMyDGGP* New();

  // Description:
  // Standard methods for printing and determining type information.
  vtkTypeMacro(vtkMyDGGP, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) VTK_OVERRIDE;

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

  //Get and set macros
  //vtkGetObjectMacro(IdList, vtkIdList);

protected:
  vtkMyDGGP();
  ~vtkMyDGGP() override;

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) VTK_OVERRIDE;

  int FillInputPortInformation(int port, vtkInformation* info) VTK_OVERRIDE;

  //member variables

  //user-selected points
  vtkIdList *UserPoints;

  //vertex ids on combined shortest path
  vtkIdList *IdList;

  int LineType;

private:
  vtkMyDGGP(const vtkMyDGGP&) = delete;
  void operator=(const vtkMyDGGP&) = delete;

};

#endif
