"""Build a ParaView state showing the DYNAMIC MATERIAL ANISOTROPY of the 3D
apical hook, cell by cell (Walia, Carter et al. 2024, Eq. 3 / Figs 3-4).

The shell surface is coloured by the stress anisotropy a = 1 - s2/s1, and a
line glyph at each cell centre points along that cell's maximal principal
stress direction — the axis the fibre (cellulose/CMT) stiffness is
redistributed onto. Glyph length scales with a, so an isotropic cell shows a
dot and a strongly anisotropic cell a long bar.

Reading the picture: in the closed hook the bars wrap AROUND the tube
(circumferential principal stress, matching dark CMT arrays); as the hook
opens the outer-side bars rotate to run ALONG the tube (longitudinal), while
the inner side stays circumferential — the inner/outer switch asymmetry the
paper reports.

Cell variable layout written by the simulation:
    4  auxin          5  cell type        6  hook angle (deg)
    7  stress anisotropy a                8,9,10  principal stress direction
    11 sigma1 (maximal principal stress)

Usage:  pvpython pvhook3d_aniso.py [t_end_hours]
"""
import glob
import os
import sys

from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()

base = os.path.dirname(os.path.abspath(__file__))
T_END = float(sys.argv[1]) if len(sys.argv) > 1 else 10.0
files = sorted(glob.glob(base + "/vtk3d/VTK_cells*.vtu"))
assert files, "no VTK_cells*.vtu frames in ./vtk3d"

ANISO = "cell variable 7"
DIRX, DIRY, DIRZ = "cell variable 8", "cell variable 9", "cell variable 10"
SIGMA1 = "cell variable 11"
AUXIN = "cell variable 4"

cells = XMLUnstructuredGridReader(registrationName="Hook shell", FileName=files)
cells.CellArrayStatus = [ANISO, DIRX, DIRY, DIRZ, SIGMA1, AUXIN,
                         "cell variable 5", "cell variable 6"]

view = GetActiveViewOrCreate("RenderView")
view.OrientationAxesVisibility = 1
view.Background = [0.12, 0.12, 0.14]

# --- shell surface, coloured by stress anisotropy ----------------------------
shell = Show(cells, view)
shell.Representation = "Surface With Edges"
shell.EdgeColor = [0.25, 0.25, 0.3]
ColorBy(shell, ("CELLS", ANISO))
aLut = GetColorTransferFunction(ANISO.replace(" ", ""))
aLut.ApplyPreset("Inferno (matplotlib)", True)
aLut.RescaleTransferFunction(0.0, 1.0)
shell.SetScalarBarVisibility(view, True)
aBar = GetScalarBar(aLut, view)
aBar.Title = "Stress anisotropy  a = 1 - s2/s1"
aBar.ComponentTitle = ""
aBar.WindowLocation = "Lower Right Corner"
aBar.ScalarBarLength = 0.33

# --- principal stress direction as a line glyph per cell ---------------------
# Cell centres carry the cell data as point data, which the Calculator can
# assemble into the direction vector the Glyph filter orients by.
centers = CellCenters(registrationName="Cell centres", Input=cells)
direction = Calculator(registrationName="Principal stress axis", Input=centers)
direction.AttributeType = "Point Data"
direction.ResultArrayName = "stressAxis"
direction.Function = (f'("{DIRX}")*iHat + ("{DIRY}")*jHat + ("{DIRZ}")*kHat')

glyph = Glyph(registrationName="Anisotropy axes", Input=direction,
              GlyphType="Line")
glyph.OrientationArray = ["POINTS", "stressAxis"]
glyph.ScaleArray = ["POINTS", ANISO]   # length proportional to anisotropy
glyph.ScaleFactor = 55.0               # um per unit anisotropy
glyph.GlyphMode = "All Points"
glyphDisplay = Show(glyph, view)
ColorBy(glyphDisplay, None)
glyphDisplay.AmbientColor = [1.0, 1.0, 1.0]
glyphDisplay.DiffuseColor = [1.0, 1.0, 1.0]
glyphDisplay.LineWidth = 3.0

# --- time annotation ----------------------------------------------------------
scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
steps = list(GetTimeKeeper().TimestepValues)
annotate = AnnotateTimeFilter(registrationName="Time", Input=cells)
annotate.Scale = T_END / max(1.0, float(len(steps) - 1))
annotate.Format = "hours after illumination = {time:.1f}"
aDisp = Show(annotate, view)
aDisp.FontSize = 20
aDisp.Bold = 1
aDisp.WindowLocation = "Upper Center"

# Frame the fully opened shape so the animation never leaves the view.
scene.AnimationTime = steps[-1]
view.ResetCamera()
camera = GetActiveCamera()
camera.Azimuth(20)
camera.Elevation(18)
camera.Dolly(0.85)
scene.AnimationTime = steps[0]
Render()

SaveState(base + "/hook3d_anisotropy.pvsm")
print("state written:", base + "/hook3d_anisotropy.pvsm")

view.ViewSize = [1200, 950]
for i, name in ((0, "t0"), (len(steps) // 4, "tA"), (len(steps) // 2, "tB"),
                (len(steps) - 1, "tEnd")):
    scene.AnimationTime = steps[i]
    Render()
    SaveScreenshot(base + f"/aniso_{name}.png", view)
print("screenshots written")
