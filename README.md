# Surface-Cutting

## Introduction ##
This repository contains the source code for the surface cutting plugin project for Paraview. All plugins are developed using Paraview's Plugin Development framework, which is described in greater detail here: https://www.paraview.org/Wiki/ParaView/Plugin_HowTo#Adding_plugins_to_ParaView_source

Developing custom plugins for Paraview requires a Paraview version that is compiled and built from source. Details on how to do that are here: https://www.paraview.org/Wiki/ParaView:Build_And_Install  

The surface cutting project consists of two different parts (and hence two main plugins):
  1. **Surface Tracking**: This filter allows the user to select points on the brain surface and connects them with a continuous curve to form a closed loop.
  2. **Surface Cutting**: Given a closed loop of points (obtained using the surface tracking plugin), this plugin extracts the user-selected region formed by the closed loop as a separate PolyData object. Currently, there are two different implementations of this feature which are described in greater detail below.

## Plugin Overview ##

### SurfaceTracker-TextEntry ###
This plugin connects two user inputted points and has three different modes which determines how the points are connected:
  1. **Geodesic**: uses Dijkstra's algorithm - simple shortest path dynamic programming approach
  2. **Gyrus**: uses Dijkstra's algorithm with a modified cost function to ensure that the generated curve is along the Gyri of the brain
  3. **Sulcus**: similar idea as Gyrus but for sulci  
  
  **Input**: 
  - pipeline browser: Brain surface that is being processed (vtkPolyData); usually rendered from a BYU file  
  - property panel: 
    - set start vertex (double)
    - set end vertex (double)
    - set line type (geodesic, gyrus, and sulcus)  
    
  **Output**: A vtkPolyData that is the brain surface with curve drawn. 
  
[vtkDijkstraGraphGeodesicPath1.h](SurfaceTracker-TextEntry/vtkDijkstraGraphGeodesicPath1.h): Header file for Surface Tracker Plugin (text entry mode)  
[vtkDijkstraGraphGeodesicPath1.cxx](SurfaceTracker-TextEntry/vtkDijkstraGraphGeodesicPath1.cxx): Source file for Surface Tracker Plugin (text entry mode)
[vtkDijkstraGraphGeodesicPath1.xml](SurfaceTracker-TextEntry/vtkDijkstraGraphGeodesicPath1.xml): XML configuration (determines what the property panel looks like and what the user can change)

### SurfaceTracker-ManualSelection ###

### SurfaceCut-ImplicitSelectionLoop ###

### SurfaceCut-ConnectedComponents ###

## Compiling Custom Plugins ##

All of the above plugins can be built using these commands:
