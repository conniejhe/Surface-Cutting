#### import the simple module from the paraview
from paraview.simple import *
import sys
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Get Active Source
proxy = GetActiveSource()

if proxy is None:
   print("Proxy is None");

active_selection = proxy.GetSelectionInput(proxy.Port)

if active_selection is None:
   print("No selection is active");

#initialize array for selected points
length = len(active_selection.IDs)

id1 = active_selection.IDs[1]

id2 = active_selection.IDs[3]
# create a new 'Extract Selection'

query = ""

for i in range(1, length, 2):
    query = query + "(id == " + str(active_selection.IDs[i]) + ")"
    if i != (length - 1):
        query = query + " |"

selection = SelectPoints(query)
extractSelection1 = ExtractSelection(Input=proxy,
    Selection=selection)

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1613, 731]

# show data in view
extractSelection1Display = Show(extractSelection1, renderView1)

# trace defaults for the display properties.
extractSelection1Display.Representation = 'Surface'

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(proxy)

# show data in view
proxyDisplay = Show(proxy, renderView1)

# trace defaults for the display properties.
proxyDisplay.Representation = 'Surface'

# create a new 'Surface Tracker Manual'
surfaceTrackerManual1 = SurfaceTrackerManual(Input=proxy,
    Selection=extractSelection1)

# show data in view
surfaceTrackerManual1Display = Show(surfaceTrackerManual1, renderView1)

# trace defaults for the display properties.
surfaceTrackerManual1Display.Representation = 'Surface'

# hide data in view
Hide(extractSelection1, renderView1)

SetActiveSource(proxy)

# update the view to ensure updated data information
renderView1.Update()
