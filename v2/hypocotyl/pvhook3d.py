# ParaView state for the 3D apical hook shell: cells colored by auxin,
# 3/4 view, time annotation in hours. Run with pvpython after generating
# VTK output into ./vtk3d.
import glob
import os
import sys

from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()

base = os.path.dirname(os.path.abspath(__file__))
T_END = float(sys.argv[1]) if len(sys.argv) > 1 else 10.0
cells = XMLUnstructuredGridReader(registrationName='Hook shell (auxin)',
                                  FileName=sorted(glob.glob(base + "/vtk3d/VTK_cells*.vtu")))
view = GetActiveViewOrCreate('RenderView')
view.OrientationAxesVisibility = 1

cd = Show(cells, view)
cd.Representation = 'Surface With Edges'
ColorBy(cd, ('CELLS', 'cell variable 4'))
lut = GetColorTransferFunction('cellvariable4')
lut.ApplyPreset('Viridis (matplotlib)', True)
lut.RescaleTransferFunction(0.0, 1.0)
cd.SetScalarBarVisibility(view, True)
bar = GetScalarBar(lut, view)
bar.Title = 'Auxin (relative)'
bar.ComponentTitle = ''
bar.WindowLocation = 'Lower Right Corner'

scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
steps = list(GetTimeKeeper().TimestepValues)
annotate = AnnotateTimeFilter(registrationName='Time', Input=cells)
annotate.Scale = T_END / max(1.0, float(len(steps) - 1))
annotate.Format = 'hours after illumination = {time:.1f}'
ad = Show(annotate, view)
ad.FontSize = 20
ad.Bold = 1
ad.WindowLocation = 'Upper Center'

scene.AnimationTime = steps[-1]
view.ResetCamera()
Render()
scene.AnimationTime = steps[0]
Render()
SaveState(base + '/hook3d_opening.pvsm')
view.ViewSize = [1100, 900]
camera = GetActiveCamera()
camera.Azimuth(25)
camera.Elevation(20)
for i, name in ((0, 't0'), (len(steps) // 3, 'tA'), (2 * len(steps) // 3, 'tB'),
                (len(steps) - 1, 'tEnd')):
    scene.AnimationTime = steps[i]
    Render()
    SaveScreenshot(base + f"/hook3d_{name}.png", view)
print('state written')
