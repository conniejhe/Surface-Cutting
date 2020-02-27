Below is a step-by-step procedure demonstrating how to use the surface tracking and cutting features to extract a custom region of a brain surface. Compile and load the 'SurfaceTracker-ManualSelection' and 'SurfaceCut-BoundaryFill' plugins beforehand.

This example uses the byu file named 'example_surface.byu'.

1. Load the 'example_surface.byu' file using the BYU Reader.
2. Using the 'Interactive Select Points Tool', select two points along the sulcus. Use the 'Extract Selection' tool to save these two points.
  
![Extracted Selection](./ExamplePictures/ExtractSelection1.png)
3. Select the 'Surface Tracker Manual' filter from the Filters dropdown menu and select 'example_surface.byu' as the Input and 'ExtractSelection1' as the Selection. Switch the line type to sulcus.
  
![Sulcus Line](./ExamplePictures/SurfaceTracker1.png)
4. Using the 'Interactive Select Points Tool', deselect the first point (using Shift key) and select a point at the end of the gyrus. Use the 'Extract Selection' tool to save these two points.
  
![Extracted Selection](./ExamplePictures/ExtractSelection2.png)
5. Select the 'Surface Tracker Manual' filter from the Filters dropdown menu and select 'example_surface.byu' as the Input and 'ExtractSelection2' as the Selection. Switch the line type to gyrus.
  
![Gyrus Line](./ExamplePictures/SurfaceTracker2.png)
6. Using the 'Interactive Select Points Tool', deselect the second point and select the first point originally selected (to form a closed loop of vertices). Use the 'Extract Selection' tool to save these two points.
  
![Extracted Selection](./ExamplePictures/ExtractSelection3.png)
7. Select the 'Surface Tracker Manual' filter from the Filters dropdown menu and select 'example_surface.byu' as the Input and 'ExtractSelection3' as the Selection. Switch the line type to geodesic.
  
![Geodesic Line](./ExamplePictures/SurfaceTracker3.png)
8. Make sure that the three 'SurfaceTrackerManual' objects are selected simultaneously. Now that you have the individual three segments that comprise the loop, combine them using the 'Append Geometry' filter in the Filters dropdown menu.
  
![Combined Lines](./ExamplePictures/AppendGeometry.png)
9. Select a point inside the desired region which will act as the starting point for the surface cut algorithm. Extract it using the 'Extract Selection' tool. Make sure all the other points are deselected before this.
  
![Inside Point](./ExamplePictures/InsidePoint.png)
10. Select the 'Surface Cut' filter from the Filters dropdown menu and choose 'example_surface.byu' as Input, 'AppendGeometry1' as Path, and 'ExtractSelection4' as Inside Point.
  
![Surface Cut](./ExamplePictures/SurfaceCut.png)

