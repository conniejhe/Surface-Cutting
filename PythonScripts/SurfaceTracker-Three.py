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

id1 = active_selection.IDs[1]

id2 = active_selection.IDs[3]

id3 = active_selection.IDs[5]

# create a new 'Extract Selection'

query1 = "(id == " + str(id1) + ") | " + "(id == " + str(id2) + ") "
selection = SelectPoints(query = query1)
extractSelection1 = ExtractSelection(Input=proxy,
    Selection=selection)

query2 = "(id == " + str(id2) + ") | " + "(id == " + str(id3) + ") "
selection2 = SelectPoints(query = query2)
extractSelection2 = ExtractSelection(Input=proxy,
    Selection=selection2)

query3 = "(id == " + str(id3) + ") | " + "(id == " + str(id1) + ") "
selection3 = SelectPoints(query = query3)
extractSelection3 = ExtractSelection(Input=proxy,
    Selection=selection3)

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
extractSelection1Display = Show(extractSelection1, renderView1)

# trace defaults for the display properties.
extractSelection1Display.Representation = 'Surface'

# show data in view
extractSelection2Display = Show(extractSelection2, renderView1)

# trace defaults for the display properties.
extractSelection2Display.Representation = 'Surface'

# show data in view
extractSelection3Display = Show(extractSelection3, renderView1)

# trace defaults for the display properties.
extractSelection3Display.Representation = 'Surface'

# hide data in view
Hide(proxy, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(proxy)

# show data in view
sphere1Display = Show(proxy, renderView1)

# trace defaults for the display properties.
sphere1Display.Representation = 'Surface'

# create a new 'Surface Tracker Manual'
surfaceTrackerManual1 = SurfaceTrackerManual(Input=proxy,
    Selection=extractSelection1)

surfaceTrackerManual1.LineType = "Sulcus"

# show data in view
surfaceTrackerManual1Display = Show(surfaceTrackerManual1, renderView1)

# trace defaults for the display properties.
surfaceTrackerManual1Display.Representation = 'Surface'

# create a new 'Surface Tracker Manual'
surfaceTrackerManual2 = SurfaceTrackerManual(Input=proxy,
    Selection=extractSelection2)

surfaceTrackerManual2.LineType = "Gyrus"

# show data in view
surfaceTrackerManual2Display = Show(surfaceTrackerManual2, renderView1)

# trace defaults for the display properties.
surfaceTrackerManual2Display.Representation = 'Surface'

# create a new 'Surface Tracker Manual'
surfaceTrackerManual3 = SurfaceTrackerManual(Input=proxy,
    Selection=extractSelection3)

# show data in view
surfaceTrackerManual3Display = Show(surfaceTrackerManual3, renderView1)

# trace defaults for the display properties.
surfaceTrackerManual3Display.Representation = 'Surface'

# hide data in view
Hide(proxy, renderView1)

# hide data in view
Hide(extractSelection1, renderView1)

Hide(extractSelection2, renderView1)

Hide(extractSelection3, renderView1)

# update the view to ensure updated data information
renderView1.Update()
