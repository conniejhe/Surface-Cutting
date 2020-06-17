# Surface-Cutting

## Table of Contents ##
- [Surface-Cutting](#surface-cutting)
  * [Introduction](#introduction)
  * [Dependencies and Version Compatibility](#dependencies-and-version-compatibility)
  * [Filters Overview](#filters-overview)
  * [Compiling and Loading Custom Plugins](#compiling-and-loading-custom-plugins)


## Introduction ##
This repository contains the source code for the surface cutting modules developed for Paraview. The following filters have been developed specifically for use by the Center of Imaging Science. All plugins were developed using Paraview's Plugin Development framework, which is described in greater detail here: [Plugins (v5.6.0)](https://www.paraview.org/Wiki/ParaView/Plugin_HowTo) and [Plugins (latest version)](https://kitware.github.io/paraview-docs/nightly/cxx/PluginHowto.html).

Developing and using custom plugins for Paraview requires a Paraview version that is compiled and built from source. Details on how to do that are here: [Build and Install](https://www.paraview.org/Wiki/ParaView:Build_And_Install). Make sure to turn the 'PARAVIEW_USE_QT' and 'PARAVIEW_USE_PYTHON' variables on in CMake.

CIS users can access Paraview, with the plugins autoloaded, using this command:
```
$ paraview-connie
```
The surface cutting project consists of two different parts (and hence two main plugins):
  1. **Surface Tracking**: This filter allows the user to select points on the brain surface and connects them with a continuous curve to form a closed loop.
  2. **Surface Cutting**: Given a closed loop of points (obtained using the surface tracking plugin), this plugin extracts the user-selected region formed by the closed loop as a separate PolyData object.
  
## Dependencies and Version Compatibility ##
Currently, the plugins have been developed to be compatible with versions 5.6.0 and 5.8.0 of Paraview and for MacOS and Linux (e.g. Ubuntu) distributions. Windows is not supported as of now. The 'master' branch contains a single plugin (SurfaceTrackerCut) that combines all four filters into an all-in-one plugin package to be loaded into Paraview version 5.8.0. The four filters are separated into individual plugins on branch 'v5.8.0-separate' (compatible with Paraview version 5.8.0) as well as on branch 'v5.6.0' (compatible with Paraview version 5.6.0). Branch 'v5.6.0-original' is a copy of what the github repository looked like before I started updating the plugins to be compatible with latest version of Paraview. 

The Curvature and Surface Tracking (Text and Manual) filters require LAPACK and BLAS as dependencies, as they utilize LAPACK subroutines to perform linear algebra computations. They can be installed using homebrew. Thus, the CMakeLists.txt files should be modified to reflect the appropriate LAPACK and BLAS installation directories for your machine. The directories to modify are within the target_link_libraries command in the CMake files. For version 5.8.0, the CMakeLists.txt file for a plugin that must be modified is under the following directory: PluginName/Plugin/Filter/CMakeLists.txt. For version 5.6.0, there is only one CMakeLists.txt file per plugin and that is the file in which the file should be modified.

**Note:** You should comment out line 4 in the CMakeLists.txt file if you are not using a MacOS machine.

## Filters Overview ##
The files required for the plugin and file hierarchy are dependent on which version of Paraview you use.
Each plugin consists of (at minimum) four different files:
  1. Source code **(.cxx file)** for the filter 
  2. Source code header **(.h file)**
  3. Server Manager XML configuration **(.xml file)** which determines what properties are displayed in the property panel and the required input connections (along with what datatype they should be), essentially configures the plugin front-end
  4. **CMakeLists.txt**: file contains a set of directives and instructions describing the project's source files and targets

<details> 
 <summary> <strong> Curvature1 </strong> </summary>
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
 <summary> <strong> Surface Tracker Text </strong> </summary>
This filter connects two user-inputted points and has three different modes which determines how the points are connected:
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
    - set neighborhood depth (int)
    - set voxel dimensions [dx, dy, dz] (double)
    
  **Output**: A set of lines corresponding to the curve generated to connect the two points.
</details>

<details> 
 <summary> <strong> Surface Tracker Manual </strong> </summary>
This filter has similar functionality to the Surface Tracker Text filter. The difference is that instead of inputting the start and end vertices, the user can freely select points along the surface and this filter will connect all of these vertices with curves (the three modes of geodesic, gyrus, and sulcus exist as well).

  **Input**: 
  - *pipeline browser*: 
    - Brain surface that is being processed (vtkPolyData); usually rendered from a BYU file
    - User-selected vertices (vtkUnstructuredGrid): use the "interactive Select Points On" toolbar and select points on the surface, extract these points using the "Extract Selection" filter. This will return a vtkUnstructuredGrid that is visible in the pipeline which can be used as an input into this filter.
  - *property panel*: 
    - set line type (geodesic, gyrus, and sulcus)  
    
  **Output**: A set of lines corresponding to the curve generated, which connects the set of user-selected points.  
  
**Notes**: 
 - This is supposed to work on a string of consecutive points (not just two points) so it is better than the previous method in that way. 
 - If you want to connect different segments using two different modes (geodesic, gyrus, and sulcus), you will need to use the filter twice on two different extracted point selections and then combine them later one with the "Append Geometry" filter. However, you can use the python macros to chain some of these commands together and be more efficient.
</details>

<details> 
 <summary> <strong> Surface Cut </strong> </summary>
This filter cuts the brain surface along the curves generated by the surface tracking filter. It uses the Boundary Fill algorithm adopted from this page: https://www.geeksforgeeks.org/boundary-fill-algorithm/. It is similar to how the 'fill' command in Microsoft Paint works. We first build an adjacency list to keep track of each vertex's neighbors. Next, we start from a user-selected vertex that is inside the desired region and recursively visit neighbors until the boundary (defined by surface tracker loop) is reached. Lastly, we remove the cells that have at least one vertex that wasn't visited, and only keep the cells where all three vertices were visited.

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
Before compiling the plugins, see the [Dependencies and Version Compatibility](#dependencies-and-version-compatibility) section for required dependencies and modifications to the CMakeLists.txt file that should be made.

Before using one of these plugins, you must compile them into a dynamic library (.so, .dylib, or .dll).
All of the above plugins can be built using these terminal commands, starting from the plugin folder (e.g. /Surface-Cutting/PLUGIN-NAME).

```
$ mkdir build
$ cd build
$ ccmake ..
```

The ccmake interface is an iterative process in which you set the settings and run configure (c key), repeating until all values are set. Then, you can generate (g key) the make files.

The ccmake interface may indicate that you need to set ParaView_DIR (path to paraview build directory) and Qt_DIR (path to Qt library). For example, the file paths specified on bofur to build the plugins originally are as follows:
 1. **ParaView_DIR:** /export/bofur/akolasny/paraview/paraview_build
 2. **Qt_DIR:** /usr/local/qt/Qt-5.11.2/5.11.2/gcc_64/lib/cmake/Qt5
 
 Once you generate the makefiles without error, you can proceed with the following command which will build the plugin in a loadable format.
 ```
 $ make
 ```
 
These commands only need to be executed the first time you compile a plugin. Otherwise, after making modifications to the code you can compile again by just running "make" in the build folder. 

**To load the plugin:**

From the Paraview menu, select Tools => Manage Plugins... => Load New...

**Version 5.6.0:**
Navigate into the build folder for plugin and select the .so or .dylib file, and it should be good to go. You can access the filter from the Filters Dropdown menu.

**Version 5.8.0:**
Starting from the build folder, the loadable plugin file is located under this file path: 
```
/lib/paraview-5.8/plugins/PluginName/PluginName.so
```

**Note:** The name of the filter that will appear in the Filters menu matches the name of each section in the [Filters Overview](#filters-overview) section.

------------------------------------------------------------------------------------------------------------------------------
Last Edited by: Connie He (2020 June 11)
