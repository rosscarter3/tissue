"""ParaView state for the biological apical hook run.

Cell variable layout written by the simulation:
    0,1,2  maximal principal stress direction  -> VTK "cell vector"
    3      stress anisotropy a = 1 - s2/s1     -> VTK "cell vector length"
    4      sigma1      5  auxin      6  apoplastic acidification
    7      cell type   8  theta (deg)

Because the stress direction lands on VTK's native cell-vector slot, the
predicted CMT axis can be glyphed directly - no Calculator needed.

Usage:  pvpython paraview_state.py [t_end_hours] [vtk_dir]
"""
import glob
import os
import sys

from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()

base = os.path.dirname(os.path.abspath(__file__))
T_END = float(sys.argv[1]) if len(sys.argv) > 1 else 12.0
VTK = sys.argv[2] if len(sys.argv) > 2 else "vtk"
files = sorted(glob.glob(os.path.join(base, VTK, "VTK_cells*.vtu")))
assert files, f"no VTK_cells*.vtu in {VTK}"

AUXIN, ACID = "cell variable 5", "cell variable 6"
ANISO, SDIR = "cell vector length", "cell vector"

cells = XMLUnstructuredGridReader(registrationName="Hook epidermis", FileName=files)

view = GetActiveViewOrCreate("RenderView")
view.OrientationAxesVisibility = 1
view.Background = [0.12, 0.12, 0.14]

# Surface coloured by auxin: the gate that holds the hook shut in the dark.
shell = Show(cells, view)
shell.Representation = "Surface With Edges"
shell.EdgeColor = [0.25, 0.25, 0.3]
ColorBy(shell, ("CELLS", AUXIN))
lut = GetColorTransferFunction(AUXIN.replace(" ", ""))
lut.ApplyPreset("Viridis (matplotlib)", True)
lut.RescaleTransferFunction(0.0, 1.0)
shell.SetScalarBarVisibility(view, True)
bar = GetScalarBar(lut, view)
bar.Title = "Auxin (relative)"
bar.ComponentTitle = ""
bar.WindowLocation = "Lower Left Corner"
bar.ScalarBarLength = 0.3

# Predicted CMT axis: line glyph along the maximal principal stress direction,
# length scaled by the stress anisotropy.
centers = CellCenters(registrationName="Cell centres", Input=cells)
glyph = Glyph(registrationName="CMT axes", Input=centers, GlyphType="Line")
glyph.OrientationArray = ["POINTS", SDIR]
glyph.ScaleArray = ["POINTS", ANISO]
glyph.ScaleFactor = 42.0
glyph.GlyphMode = "All Points"
gd = Show(glyph, view)
ColorBy(gd, None)
gd.AmbientColor = [1.0, 1.0, 1.0]
gd.DiffuseColor = [1.0, 1.0, 1.0]
gd.LineWidth = 3.0

scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
steps = list(GetTimeKeeper().TimestepValues)
annotate = AnnotateTimeFilter(registrationName="Time", Input=cells)
annotate.Scale = T_END / max(1.0, float(len(steps) - 1))
annotate.Format = "hours after illumination = {time:.1f}"
ad = Show(annotate, view)
ad.FontSize = 20
ad.Bold = 1
ad.WindowLocation = "Upper Center"

scene.AnimationTime = steps[-1]
view.ResetCamera()
cam = GetActiveCamera()
cam.Azimuth(20)
cam.Elevation(18)
cam.Dolly(0.85)
scene.AnimationTime = steps[0]
Render()

SaveState(os.path.join(base, "hook.pvsm"))
print("state written: hook.pvsm")

view.ViewSize = [1200, 950]
os.makedirs(os.path.join(base, "figures"), exist_ok=True)
for i, name in ((0, "t0"), (len(steps) // 4, "tA"), (len(steps) // 2, "tB"),
                (len(steps) - 1, "tEnd")):
    scene.AnimationTime = steps[i]
    Render()
    SaveScreenshot(os.path.join(base, "figures", f"render_{name}.png"), view)
print("renders written to figures/")
