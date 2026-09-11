import glob, os
from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
base = os.path.dirname(os.path.abspath(__file__))
cells = XMLUnstructuredGridReader(registrationName='Cells (auxin)',
                                  FileName=sorted(glob.glob(base+"/vtk/VTK_cells*.vtu")))
walls = XMLUnstructuredGridReader(registrationName='Walls (tension)',
                                  FileName=sorted(glob.glob(base+"/vtk/VTK_walls*.vtu")))
view = GetActiveViewOrCreate('RenderView'); view.InteractionMode='2D'
view.OrientationAxesVisibility = 0
cd = Show(cells, view); ColorBy(cd, ('CELLS','cell variable 4'))
lut = GetColorTransferFunction('cellvariable4'); lut.ApplyPreset('Viridis (matplotlib)', True)
lut.RescaleTransferFunction(0.0, 1.0)
cd.SetScalarBarVisibility(view, True)
bar = GetScalarBar(lut, view); bar.Title='Auxin (relative)'; bar.ComponentTitle=''
bar.WindowLocation = 'Lower Left Corner'
wd = Show(walls, view); ColorBy(wd, ('CELLS','wall variable 0'))
wlut = GetColorTransferFunction('wallvariable0'); wlut.ApplyPreset('Black-Body Radiation', True)
wlut.RescaleTransferFunction(0.0, 400.0)
wd.SetScalarBarVisibility(view, True)
wbar = GetScalarBar(wlut, view); wbar.Title='Wall tension'; wbar.ComponentTitle=''
wbar.WindowLocation = 'Lower Right Corner'
annotate = AnnotateTimeFilter(registrationName='Time', Input=cells)
annotate.Scale = 16.0 / 48.0  # frames -> hours
annotate.Format = 'hours after illumination = {time:.1f}'
ad = Show(annotate, view); ad.FontSize = 20; ad.Bold = 1; ad.WindowLocation = 'Upper Center'
scene = GetAnimationScene(); scene.UpdateAnimationUsingDataTimeSteps()
steps = list(GetTimeKeeper().TimestepValues)
scene.AnimationTime = steps[-1]; view.ResetCamera(); Render()
scene.AnimationTime = steps[0]; Render()
SaveState(base+'/hook_opening.pvsm')
view.ViewSize=[1000,850]
for i, name in ((0,'t0'), (len(steps)//3,'tA'), (2*len(steps)//3,'tB'), (len(steps)-1,'tEnd')):
    scene.AnimationTime = steps[i]; Render()
    SaveScreenshot(base+f"/hookfinal_{name}.png", view)
print('state written')
