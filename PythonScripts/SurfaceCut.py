#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
proxy = GetActiveSource()

if proxy is None:
   print("Proxy is None")

active_selection = proxy.GetSelectionInput(proxy.Port)

if active_selection is None:
   print("No selection is active")

#initialize array for selected points
length = len(active_selection.IDs)

id1 = active_selection.IDs[1]

query1 = "(id == " + str(id1) + ") "
selection = SelectPoints(query = query1)
extractSelection4 = ExtractSelection(Input=proxy,
    Selection=selection)

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

s1 = FindSource('SurfaceTrackerManual1')

segments = []
count = 2

while s1 is not None:
    segments.append(s1)
    nextSeg = 'SurfaceTrackerManual' + str(count)
    s1 = FindSource(nextSeg)
    count += 1

# create a new 'Append Geometry'
appendGeometry1 = AppendGeometry(Input=segments)

# show data in view
appendGeometry1Display = Show(appendGeometry1, renderView1)

# trace defaults for the display properties.
appendGeometry1Display.Representation = 'Surface'

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(proxy)

# create a new 'Surface Cut'
surfaceCut1 = SurfaceCut(Input=proxy,
    Path=appendGeometry1,
    InsidePoint=extractSelection4)

# show data in view
surfaceCut1Display = Show(surfaceCut1, renderView1)

# trace defaults for the display properties.
surfaceCut1Display.Representation = 'Surface'

# hide data in view
Hide(proxy)

# hide data in view
Hide(extractSelection4, renderView1)

# hide data in view
Hide(appendGeometry1, renderView1)

# update the view to ensure updated data information
renderView1.Update()
