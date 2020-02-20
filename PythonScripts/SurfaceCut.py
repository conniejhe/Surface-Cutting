# trace generated using paraview version 5.6.0-RC2-140-g920f9a1
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
surfaceTrackerManual1 = FindSource('SurfaceTrackerManual1')

# set active source
SetActiveSource(surfaceTrackerManual1)

# find source
surfaceTrackerManual2 = FindSource('SurfaceTrackerManual2')

# find source
surfaceTrackerManual3 = FindSource('SurfaceTrackerManual3')

# create a new 'Append Geometry'
appendGeometry1 = AppendGeometry(Input=[surfaceTrackerManual1, surfaceTrackerManual2, surfaceTrackerManual3])

# find source
cON44_lh_PT_surf_child_05byu = FindSource('CON44_lh_PT_surf_child_0.5.byu')

# find source
extractSelection1 = FindSource('ExtractSelection1')

# find source
extractSelection2 = FindSource('ExtractSelection2')

# find source
extractSelection3 = FindSource('ExtractSelection3')

# find source
extractSelection4 = FindSource('ExtractSelection4')

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1562, 731]

# show data in view
appendGeometry1Display = Show(appendGeometry1, renderView1)

# trace defaults for the display properties.
appendGeometry1Display.Representation = 'Surface'

# hide data in view
Hide(surfaceTrackerManual1, renderView1)

# hide data in view
Hide(surfaceTrackerManual3, renderView1)

# hide data in view
Hide(surfaceTrackerManual2, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(cON44_lh_PT_surf_child_05byu)

# create a new 'Surface Cut'
surfaceCut1 = SurfaceCut(Input=cON44_lh_PT_surf_child_05byu,
    Path=appendGeometry1,
    InsidePoint=extractSelection4)

# show data in view
surfaceCut1Display = Show(surfaceCut1, renderView1)

# trace defaults for the display properties.
surfaceCut1Display.Representation = 'Surface'

# hide data in view
Hide(cON44_lh_PT_surf_child_05byu, renderView1)

# hide data in view
Hide(extractSelection4, renderView1)

# hide data in view
Hide(appendGeometry1, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(cON44_lh_PT_surf_child_05byu)

# hide data in view
Hide(surfaceCut1, renderView1)