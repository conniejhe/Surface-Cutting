Below is a step-by-step procedure demonstrating how to use the surface tracking and cutting features to extract a custom region of a brain surface.

This example uses the byu file named 'example_surface.byu'.

1. Load the 'example_surface.byu' file using the BYU Reader.
2. Using the 'Interactive Select Points Tool', select two points along the sulcus. Use the 'Extract Selection' tool to save these two points.
3. Select the 'Surface Tracker Manual' filter from the Filters dropdown menu and select 'example_surface.byu' as the Input and 'ExtractSelection1' as the Selection. Switch the line type to sulcus.
4. Using the 'Interactive Select Points Tool', deselect the first point (using Shift key) and select a point at the end of the gyrus. Use the 'Extract Selection' tool to save these two points.
5. Select the 'Surface Tracker Manual' filter from the Filters dropdown menu and select 'example_surface.byu' as the Input and 'ExtractSelection2' as the Selection. Switch the line type to gyrus.
6. Using the 'Interactive Select Points Tool', deselect the second point and select the first point originally selected (to form a closed loop of vertices). Use the 'Extract Selection' tool to save these two points.
7. Select the 'Surface Tracker Manual' filter from the Filters dropdown menu and select 'example_surface.byu' as the Input and 'ExtractSelection3' as the Selection. Switch the line type to geodesic.
8. Now that you have the individual three segments that comprise the loop, combine them using the 'Append Geometry' filter in the Filters dropdown menu. Make sure that the three 'SurfaceTrackerManual' objects are selected simultaneously before performing this step.
9. Select a point inside the desired region which will act as the starting point for the surface cut algorithm.

