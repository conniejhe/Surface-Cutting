#ifndef __surfaceCut_h
#define __surfaceCut_h

#include "vtkPolyDataAlgorithm.h"
#include <vector>

using std::vector;

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

    // void splitComponents(vtkPolyData *in);
    //
    // vtkIdList* findReachableNodes(vtkPolyData* in);
    //
    // vtkIdList* BFS(int compNum, vtkIdType src, vtkIntArray* visited);
    //
    // void cut(vtkPolyData* in, vtkPolyData* out, vtkIdList* list);

    void BuildAdjacency(vtkDataSet *inData);

    void ColorBoundary();

    void FillBoundary(vtkDataSet *inData, vtkIdType i,
      int bound_color, int fill_color);

    void CutSurface(vtkPolyData* in, vtkPolyData* out);

    //member variables

    //user-selected points
    vtkIdList *UserPoints;

    vtkIdType insidePoint;

    int NumberOfVertices;

    vector<vtkSmartPointer<vtkIdList>> adjacencyMatrix;

    vtkIntArray* colorArray;

private:
    surfaceCut(const surfaceCut&) = delete;
    void operator=(const surfaceCut&) = delete;

};

#endif
