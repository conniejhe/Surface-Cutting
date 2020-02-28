# Surface-Cutting

## Table of Contents ##
- [Surface-Cutting](#surface-cutting)
  * [Introduction](#introduction)
  * [Plugins Overview](#plugins-overview)
    + [Curvature](#Curvature)
    + [SurfaceTracker-TextEntry](#surfacetracker-textentry)
    + [SurfaceTracker-ManualSelection](#surfacetracker-manualselection)
    + [SurfaceCut-ImplicitSelectionLoop](#surfacecut-implicitselectionloop)
    + [SurfaceCut-ConnectedComponents](#surfacecut-connectedcomponents)
    + [SurfaceCut-BoundaryFill](#surfacecut-boundaryfill)
  * [Compiling and Loading Custom Plugins](#compiling-and-loading-custom-plugins)

## Introduction ##
This repository contains the source code for the surface cutting plugin project for Paraview. The following plugins have been developed specifically for use by the Center of Imaging Science. All plugins were developed using Paraview's Plugin Development framework, which is described in greater detail here: https://www.paraview.org/Wiki/ParaView/Plugin_HowTo#Adding_plugins_to_ParaView_source

Developing custom plugins for Paraview requires a Paraview version that is compiled and built from source. Details on how to do that are here: https://www.paraview.org/Wiki/ParaView:Build_And_Install

All development thus far has been done on Bofur (CIS machine), using the version of Paraview located in this path:   
```
/export/bofur/akolasny/paraview/paraview_build/bin/paraview
```
CIS users can access this version of Paraview (with the plugins autoloaded) using this command:
```
paraview-connie
```

The surface cutting project consists of two different parts (and hence two main plugins):
  1. **Surface Tracking**: This filter allows the user to select points on the brain surface and connects them with a continuous curve to form a closed loop.
  2. **Surface Cutting**: Given a closed loop of points (obtained using the surface tracking plugin), this plugin extracts the user-selected region formed by the closed loop as a separate PolyData object.
  
All questions should be directed to Connie He (<conbonhe98@gmail.com>).

## Plugins Overview ##
Each plugin consists of (at minimum) four different files:
  1. Source code (.cxx file) for the filter 
  2. Source code header (.h)
  3. Server Manager XML configuration (.xml) file which determines what properties are displayed in the property panel and the required input connections (along with what datatype they should be), essentially configures the plugin front-end
  4. CMakeLists.txt: basically a make file which tells cmake which classes to build
  
### Curvature ###
This plugin calculates the curvature of the brain surface at each vertex. It implements Hamann's algorithm, which derives a tangent plane at each point and uses it to compute the local shape operator.

  **Input**: 
  - *pipeline browser*: Brain surface that is being processed (vtkPolyData); usually rendered from a BYU file  
  - *property panel*: 
    - set curvature calculation type (mean, gauss, max, min)
    - set neighorhood depth (int)
    - set voxel dimensions [dx, dy, dz] (double)
    
  **Output**: A vtkPolyData object that is the same as the input but with an extra Curvature array in Point Data.
 
### SurfaceTracker-TextEntry ###
This plugin connects two user-inputted points and has three different modes which determines how the points are connected:
  1. **Geodesic**: uses Dijkstra's algorithm - simple shortest path dynamic programming approach
  2. **Gyrus**: uses Dijkstra's algorithm with a modified cost function to ensure that the generated curve is along the Gyri of the brain
  3. **Sulcus**: similar idea as Gyrus but for sulci

Paraview has an existing filter that calculates geodesic shortest path using Dijkstra's algorithm on a graph which much of the code was based on. The original filter source code can be found here: https://github.com/Kitware/VTK/blob/master/Filters/Modeling/vtkDijkstraGraphGeodesicPath.cxx, along with the documentation here: https://vtk.org/doc/nightly/html/classvtkDijkstraGraphGeodesicPath.html
  
  **Input**: 
  - *pipeline browser*: Brain surface that is being processed (vtkPolyData); usually rendered from a BYU file  
  - *property panel*: 
    - set start vertex (double)
    - set end vertex (double)
    - set line type (geodesic, gyrus, and sulcus) 
    - set curvature type used for cost function (mean, gauss, max, min)
      - may change this later because should really only be using max curvature
    - set neighborhood depth (int)
    - set voxel dimensions [dx, dy, dz] (double)
    
  **Output**: A set of lines corresponding to the curve generated to connect the two points.
  
**Notes**: 
 - May consider adding curvature filter output as an input into this plugin to reduce redundancy.

### SurfaceTracker-ManualSelection ###
This plugin has similar functionality to the SurfaceTracker-TextEntry plugin. The difference is that instead of inputting the start and end vertices, the user can freely select points along the surface and this filter will connect all of these vertices with curves (the three modes of geodesic, gyrus, and sulcus exist as well).

  **Input**: 
  - *pipeline browser*: 
    - Brain surface that is being processed (vtkPolyData); usually rendered from a BYU file
    - User-selected vertices (vtkUnstructuredGrid): use the "interactive Select Points On" toolbar and select points on the surface, extract these points using the "Extract Selection" filter. This will return a vtkUnstructuredGrid that is visible in the pipeline which can be used as an input into this filter.
  - *property panel*: 
    - set line type (geodesic, gyrus, and sulcus)  
    
  **Output**: A set of lines corresponding to the curve generated to connect the set of points.  
  
**Notes**: 
 - This is supposed to work on a string of consecutive points (not just two points) so it is better than the previous method in that way. 
 - If you want to connect different segments using two different modes (geodesic, gyrus, and sulcus), you will need to use the filter twice on two different extracted point selections and then combine them later one with the "Append Geometry" filter (this is a bit of a hassle, so I will work on finding a better solution for this)
 - May set curvature filter output as an input to this plugin as well to reduce redundancy (of property panel attributes)

### SurfaceCut-ImplicitSelectionLoop ###
This plugin cuts the brain surface along the curves generated by the surface tracking filter. It uses vtk filter Implicit Selection Loop along with the clip filter to extract the inner region of the loop. Documentation for implicit selection loop can be found here: https://vtk.org/doc/nightly/html/classvtkImplicitSelectionLoop.html, and this is an example using it: https://lorensen.github.io/VTKExamples/site/Cxx/PolyData/ImplicitSelectionLoop/  

  **Input**: 
  - *pipeline browser*: 
    - Brain surface that is being processed (vtkPolyData); usually rendered from a BYU file
    - Output of the Surface Tracker filter (vtkPolyData) that consists of a closed loop of points  
    
  **Output**: A vtkPolyData that is the extracted region.
  
**Notes**: 
 - This method currently doesn't work too well. It relies on the calculation of the implicit function value of each point. There is a clip function which when set to 0, will clip out the positive region (cells outside the loop will have positive implicit function values). However, I suspect this doesn't work too well because the surface itself is irregular and this may only work well for conical or spherical-shaped objects.
 - The results of this filter are variable, but it often doesn't capture the entire region you specified or captures too much.

### SurfaceCut-ConnectedComponents ###
This plugin is similar in functionality to the SurfaceCut-ImplicitSelectionLoop plugin, but it uses a different algorithm to extract the region inside the loop. We first build an adjacency list to keep track of each vertex's neighbors. Next, we split the brain surface into two components by deleting the cells that are in contact with the vertices along the user-specified path. The plugin also takes as input a vertex that isn't inside the loop, so that we can find all reachable vertices from that outside vertex and remove those from the graph.

  **Input**: 
  - *pipeline browser*: 
    - Brain surface that is being processed (vtkPolyData); usually rendered from a BYU file
    - Output of the Surface Tracker filter (vtkPolyData) that consists of a closed loop of points
    - An outside point (vtkUnstructuredGrid): obtained by selecting a point oustide of the desired region (using the interactive select points on tool) and extracting this selection  
    
  **Output**: A vtkPolyData that is the extracted region.
  
**Notes**:
 - This method is more consistent than the implicit selection loop + clip method, but it takes longer to run (around 30 seconds). 
 - Another downside is that it truncates the region of interest slightly, but still retains most of its general outline.
 - You can use the loop subdivison filter to increase the granularity of the mesh (split each existing triangle into more triangles), but then the plugin takes significantly longer to run.
 
 ### SurfaceCut-BoundaryFill ###
This plugin is similar in the functionality to the two previous SurfaceCut plugins, but it uses a different algorithm to extract the region inside the loop. It uses the Boundary Fill algorithm adopted from this page: https://www.geeksforgeeks.org/boundary-fill-algorithm/. It is similar to how the 'fill' command in MS Paint works. We first build an adjacency list to keep track of each vertex's neighbors. Next, we start from a user-selected vertex that is inside the desired region and recursively visit neighbors until the boundary (defined by surface tracker loop) is reached. Lastly, we remove the cells that have at least one vertex that wasn't visited, and only keep the cells where all three vertices were visited.

  **Input**: 
  - *pipeline browser*: 
    - Brain surface that is being processed (vtkPolyData); usually rendered from a BYU file
    - Output of the Surface Tracker filter (vtkPolyData) that consists of a closed loop of points
    - An inside point (vtkUnstructuredGrid): obtained by selecting a point inside of the desired region (using the interactive select points on tool) and extracting this selection  
    
  **Output**: A vtkPolyData that is the extracted region.
  
**Notes**:
 - This method is the most effective surface cut algorithm that I've experimented with so far.
 - It is effective in clipping the desired region most of time and completes within a second.

## Compiling and Loading Custom Plugins ##
Before using one of these plugins, you must compile them into a build directory.
All of the above plugins can be built using these commands:

 1. Make a new build directory in the plugin folder and navigate into it (should be empty)
    1. example: cis/home/che/Documents/CPPFilterPlugin/build
 2. Type following series of commands: ccmake .. -> c -> g -> make
    1. will show EMPTY CACHE at first
 3. need to set ParaView_DIR: /export/bofur/akolasny/paraview/paraview_build
 4. need to set Qt_DIR: /usr/local/qt/Qt-5.11.2/5.11.2/gcc_64/lib/cmake/Qt5
 
These commands only need to be executed the first time you compile a plugin. Otherwise, after making modifications to the code you can compile again by just running "make" in the build folder. 

The plugin will be compiled as a .so file that can be loaded into Paraview. From the Paraview menu, select Tools -> Manage Plugins... -> Load New...  
Navigate into build folder for plugin and select the .so file, and it should be good to go.

------------------------------------------------------------------------------------------------------------------------------
Last Edited by: Connie He (2020 Feb 27)
