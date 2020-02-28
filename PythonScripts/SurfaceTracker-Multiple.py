#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Get Active Source
proxy = GetActiveSource()

if proxy is None:
   print "Proxy is None"   

active_selection = proxy.GetSelectionInput(proxy.Port)

if active_selection is None:
   print "No selection is active"

#initialize array for selected points
length = len(active_selection.IDs)

for i in range(1, length, 2):
    id1 = active_selection.IDs[i]
    id2 = active_selection.IDs[(i+2) % length]
    query1 = "(id == " + str(id1) + ") | " + "(id == " + str(id2) + ") "
    selection = SelectPoints(query = query1)
    extractSelection = ExtractSelection(Input=proxy,
        Selection=selection)

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')

    # show data in view
    extractSelectionDisplay = Show(extractSelection, renderView1)

    # trace defaults for the display properties.
    extractSelectionDisplay.Representation = 'Surface'

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
        Selection=extractSelection)

    # show data in view
    surfaceTrackerManual1Display = Show(surfaceTrackerManual1, renderView1)

    # trace defaults for the display properties.
    surfaceTrackerManual1Display.Representation = 'Surface'

    # hide data in view
    Hide(extractSelection, renderView1)

SetActiveSource(proxy)

# update the view to ensure updated data information
renderView1.Update()
