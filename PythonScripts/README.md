The general workflow for surface tracking and cutting is described in [Example.md](../Examples/Example.md) under the Examples folder. However, there are several scripts that can be used to speed up the process and reduce the number of mouse-clicks required.

## Surface Tracking
To use any of these scripts, you should set the active source as the brain surface (selected and highlighted in the pipeline). Select three points using the 'Interactive Select Points On' tool. Make sure the correct points are selected (you know if they're selected if they are highlighted in magenta). Without extracting the points, run this script.

More specific descriptions for each script and when they should be used are below:

### SurfaceTracker-One.py
This script will call the 'Surface Tracker Manual' filter on all of these points, generating a single line connecting them. The default line type will be geodesic, but you can change it in the properties panel. You should use this script if you just want to generate lines connecting two points at a time as you proceed with surface tracking. The incremental nature of this approach allows you to change the points more easily as you go if you don't like the line segment that was drawn. Note that there is a bug where after running this script, if you try to select new points on the surface it will generate a warning and also delete the existing selection. Therefore, you will need to reselect the desired points again. 

### SurfaceTracker-Three.py
This script will call the 'Surface Tracker Manual' filter three times, to form separate segments between the 1st and 2nd point, 2nd and 3rd point, and 3rd back to the 1st point. The default line types for these are sulcus, gyrus, and geodesic, respectively, but they are modifiable using the properties panel. This is useful if you want to generate three different segments at once instead of one by one.

### SurfaceTracker-Multiple.py
This script will call the 'Surface Tracker Manual' filter multiple times, to form separate segments between each consecutive pair of points. The line type for all of these is geodesic by default but can be changed. This is useful for more complex brain surfaces, where you may want to generate more than three segments to form the boundary of the surface cut.

## Surface Cutting
### SurfaceCut.py
To use this script, make sure you have three 'SurfaceTrackerManual' objects, with the suffixes '1', '2', and '3'. This script relies on the proper naming of these segments, so it is very important that the naming is correct (which it should be by default). The active source should be the brain surface. Surface Cut requires a selected point on the inside of the desired region to indicate which portion the user would like to keep. Therefore, make sure you have an inside point selected as the only point selected - do not extract it, just leave it selected. Now, you can run the script which will first append the 'SurfaceTrackerManual' objects into one object and then run the 'Surface Cut' filter using the appended line segments and the selected inside point. This only works with three line segments, so if you want to create more during tracking you should use the regular workflow to cut the surface (as described in [Example.md](../Examples/Example.md)). 
