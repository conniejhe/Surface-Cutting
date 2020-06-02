# Surface-Cutting

## Table of Contents ##
- [Surface-Cutting](#surface-cutting)
  * [Introduction](#introduction)
  * [Version Compatibility and Dependencies](#version-compatibility-and-dependencies)
  * [Plugins Overview](#plugins-overview)
  * [Compiling and Loading Custom Plugins](#compiling-and-loading-custom-plugins)


## Introduction ##
This repository contains the source code for the surface cutting plugin project for Paraview. The following plugins have been developed specifically for use by the Center of Imaging Science. All plugins were developed using Paraview's Plugin Development framework, which is described in greater detail here: [Plugins](https://www.paraview.org/Wiki/ParaView/Plugin_HowTo#Adding_plugins_to_ParaView_source)

Developing custom plugins for Paraview requires a Paraview version that is compiled and built from source. Details on how to do that are here: [Build and Install](https://www.paraview.org/Wiki/ParaView:Build_And_Install). Make sure to turn the 'PARAVIEW_USE_QT' and 'PARAVIEW_USE_PYTHON' variables on in ccmake.

All development thus far has been done on Bofur (CIS machine), using the version of Paraview located in this path:   
```
$ /export/bofur/akolasny/paraview/paraview_build/bin/paraview
```
CIS users can access this version of Paraview (with the plugins autoloaded) using this command:
```
$ paraview-connie
```
The surface cutting project consists of two different parts (and hence two main plugins):
  1. **Surface Tracking**: This filter allows the user to select points on the brain surface and connects them with a continuous curve to form a closed loop.
  2. **Surface Cutting**: Given a closed loop of points (obtained using the surface tracking plugin), this plugin extracts the user-selected region formed by the closed loop as a separate PolyData object.
  
All questions should be directed to Connie He (<conbonhe98@gmail.com>).
  
## Version Compatibility and Dependencies ##
Currently, the plugins have been developed to be compatible with versions 5.6.0 and 5.8.0 of Paraview and for MacOS and Linux (e.g. Ubuntu) distributions. Windows is not supported as of now. The plugins on branch 'master' are compatible with Paraview version 5.8.0 (MacOS and Linux), and the plugins on branch 'v5.6.0-modified' are compatible with Paraview version 5.6.0 (MacOS). Plugins on branch 'v5.6.0' were developed for Paraview version 5.6.0 (Linux).

The Curvature and Surface Tracking (Text and Manual) plugins require Lapack and Blas as dependencies, as they utilize Lapack subroutines to perform linear algebra computations. Thus, the CMakeLists.txt files should be modified to reflect the appropriate LAPACK and BLAS installation directories for your machine. 

For version 5.8.0, the CMakeLists.txt file for the Curvature plugin that must be modified is under the following directory: Curvature/Plugin/Filter/CMakeLists.txt. It is similar for the Surface Tracker plugins. For version 5.6.0, there is only one CMakeLists.txt file per plugin and that is the file in which the file should be modified.

## Plugins Overview ##
The files required for the plugin and file hierarchy varies depending on which version of Paraview you use.
Each plugin consists of (at minimum) four different files:
  1. Source code (.cxx file) for the filter 
  2. Source code header (.h)
  3. Server Manager XML configuration (.xml) file which determines what properties are displayed in the property panel and the required input connections (along with what datatype they should be), essentially configures the plugin front-end
  4. CMakeLists.txt: basically a make file which tells cmake which classes to build
  

<details> 
 <summary> <strong> Curvature </strong> </summary>
This plugin calculates the curvature of the brain surface at each vertex. It implements Hamann's algorithm, which derives a tangent plane at each point and uses it to compute the local shape operator.

  **Input**: 
  - *pipeline browser*: Brain surface that is being processed (vtkPolyData); usually rendered from a BYU file  
  - *property panel*: 
    - set curvature calculation type (mean, gauss, max, min)
    - set neighorhood depth (int)
    - set voxel dimensions [dx, dy, dz] (double)
    
  **Output**: A vtkPolyData object that is the same as the input but with an extra Curvature array in Point Data. 

**Notes**: 
 - Requires Lapack as a dependency.
</details>

<details> 
 <summary> <strong> SurfaceTracker-TextEntry </strong> </summary>
This plugin connects two user-inputted points and has three different modes which determines how the points are connected:
  <ol>
   <li> <strong> Geodesic </strong>: uses Dijkstra's algorithm - simple shortest path dynamic programming approach </li>
   <li> <strong> Gyrus </strong>: uses Dijkstra's algorithm with a modified cost function to ensure that the generated curve is along the Gyri of the brain </li>
   <li> <strong> Sulcus </strong>: similar idea as Gyrus but for sulci </li>
 </ol>
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
</details>

<details> 
 <summary> <strong> SurfaceTracker-ManualSelection </strong> </summary>
This plugin has similar functionality to the SurfaceTracker-TextEntry plugin. The difference is that instead of inputting the start and end vertices, the user can freely select points along the surface and this filter will connect all of these vertices with curves (the three modes of geodesic, gyrus, and sulcus exist as well).

  **Input**: 
  - *pipeline browser*: 
    - Brain surface that is being processed (vtkPolyData); usually rendered from a BYU file
    - User-selected vertices (vtkUnstructuredGrid): use the "interactive Select Points On" toolbar and select points on the surface, extract these points using the "Extract Selection" filter. This will return a vtkUnstructuredGrid that is visible in the pipeline which can be used as an input into this filter.
  - *property panel*: 
    - set line type (geodesic, gyrus, and sulcus)  
    
  **Output**: A set of lines corresponding to the curve generated, which connecst the set of user-selected points.  
  
**Notes**: 
 - This is supposed to work on a string of consecutive points (not just two points) so it is better than the previous method in that way. 
 - If you want to connect different segments using two different modes (geodesic, gyrus, and sulcus), you will need to use the filter twice on two different extracted point selections and then combine them later one with the "Append Geometry" filter. However, you can use the python macros to chain some of these commands together and be more efficient.
</details>

<details> 
 <summary> <strong> SurfaceCut-ImplicitSelectionLoop (Archived) </strong> </summary>
This plugin cuts the brain surface along the curves generated by the surface tracking filter. It uses vtk filter Implicit Selection Loop along with the clip filter to extract the inner region of the loop. Documentation for implicit selection loop can be found here: https://vtk.org/doc/nightly/html/classvtkImplicitSelectionLoop.html, and this is an example using it: https://lorensen.github.io/VTKExamples/site/Cxx/PolyData/ImplicitSelectionLoop/  

  **Input**: 
  - *pipeline browser*: 
    - Brain surface that is being processed (vtkPolyData); usually rendered from a BYU file
    - Output of the Surface Tracker filter (vtkPolyData) that consists of a closed loop of points  
    
  **Output**: A vtkPolyData that is the extracted region.
  
**Notes**: 
 - This method currently doesn't work too well. It relies on the calculation of the implicit function value of each point. There is a clip function which when set to 0, will clip out the positive region (cells outside the loop will have positive implicit function values). However, I suspect this doesn't work too well because the surface itself is irregular and this may only work well for conical or spherical-shaped objects.
 - The results of this filter are variable, but it often doesn't capture the entire region you specified or captures too much.
</details>

<details> 
 <summary> <strong> SurfaceCut-ConnectedComponents (Archived) </strong> </summary>
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
 </details>
 
<details> 
 <summary> <strong> SurfaceCut-BoundaryFill </strong> </summary>
This plugin is similar in the functionality to the two previous SurfaceCut plugins, but it uses a different algorithm to extract the region inside the loop. It uses the Boundary Fill algorithm adopted from this page: https://www.geeksforgeeks.org/boundary-fill-algorithm/. It is similar to how the 'fill' command in MS Paint works. We first build an adjacency list to keep track of each vertex's neighbors. Next, we start from a user-selected vertex that is inside the desired region and recursively visit neighbors until the boundary (defined by surface tracker loop) is reached. Lastly, we remove the cells that have at least one vertex that wasn't visited, and only keep the cells where all three vertices were visited.

  **Input**: 
  - *pipeline browser*: 
    - Brain surface that is being processed (vtkPolyData); usually rendered from a BYU file
    - Output of the Surface Tracker filter (vtkPolyData) that consists of a closed loop of points
    - An inside point (vtkUnstructuredGrid): obtained by selecting a point inside of the desired region (using the interactive select points on tool) and extracting this selection  
    
  **Output**: A vtkPolyData that corresponds to the extracted region.
  
**Notes**:
 - This method is the most effective surface cut algorithm that I've experimented with so far.
 - It is effective in clipping the desired region most of the time and completes within a second.
 
</details>

## Compiling and Loading Custom Plugins ##
Before using one of these plugins, you must compile them into a dynamic library (.so, .dylib, or .dll).
All of the above plugins can be built using these commands (starting from the plugin folder).

```
$ mkdir build
$ cd build
$ ccmake ..
```

The ccmake interface is an iterative process in which you set the settings and run configure (c key), repeating until all values are set. Then, you can generate (g key) the make files.

The ccmake interface will indicate that you need to set ParaView_DIR (path to paraview build directory) and Qt_DIR (path to Qt library). For example, the file paths specified on bofur to build the plugins originally are as follows:
 1. ParaView_DIR: /export/bofur/akolasny/paraview/paraview_build
 2. Qt_DIR: /usr/local/qt/Qt-5.11.2/5.11.2/gcc_64/lib/cmake/Qt5
 
 Once you generate the makefiles without error, you can proceed with the following command which will build the plugin in a loadable format.
 ```
 $ make
 ```
 
These commands only need to be executed the first time you compile a plugin. Otherwise, after making modifications to the code you can compile again by just running "make" in the build folder. 

**To load the plugin:**
From the Paraview menu, select Tools => Manage Plugins... => Load New...

Version 5.6.0:
Navigate into the build folder for plugin and select the .so or .dylib file, and it should be good to go. You can access the filter from the Filters Dropdown menu.

Version 5.8.0:
Starting from the build folder, the loadable plugin file is located under this file path: 
```
/lib/paraview-5.8/plugins/PluginName/PluginName.so
```
------------------------------------------------------------------------------------------------------------------------------
Last Edited by: Connie He (2020 June 1)
