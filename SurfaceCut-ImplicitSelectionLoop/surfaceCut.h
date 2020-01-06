#ifndef __surfaceCut_h
#define __surfaceCut_h

#include "vtkPolyDataAlgorithm.h"

class vtkIdList;

class VTK_EXPORT surfaceCut : public vtkPolyDataAlgorithm
{

public:
  static surfaceCut* New();

  // Description:
  // Standard methods for printing and determining type information.
  vtkTypeMacro(surfaceCut, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) VTK_OVERRIDE;

  //Get and set macros
  //vtkGetObjectMacro(IdList, vtkIdList);

protected:
  surfaceCut();
  ~surfaceCut() override;

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) VTK_OVERRIDE;

  int FillInputPortInformation(int port, vtkInformation* info) VTK_OVERRIDE;

  //member variables

  //user-selected points
  vtkIdList *UserPoints;

private:
  surfaceCut(const surfaceCut&) = delete;
  void operator=(const surfaceCut&) = delete;

};

#endif
