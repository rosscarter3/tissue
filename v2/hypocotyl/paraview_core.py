"""ParaView animation of the hook opening with the CMT switch visible.

Colours the epidermis by the predicted CMT orientation and draws a line glyph
at each cell centre along the maximal principal stress direction - the axis the
model predicts cortical microtubules follow. The outer flank's glyphs rotate
from circumferential to longitudinal partway through; the inner flank's do not
until much later.

Cell variable layout: 0-2 stress direction (VTK "cell vector"), 3 anisotropy
("cell vector length"), 4 sigma1, 5 auxin, 6 acid.

Usage:
  pvpython paraview_core.py [t_end_hours] [vtk_dir] [--movie]
"""
import glob, os, sys
from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()
base = os.path.dirname(os.path.abspath(__file__))
T_END = float(sys.argv[1]) if len(sys.argv) > 1 else 12.0
VTK = sys.argv[2] if len(sys.argv) > 2 else "vtk_core"
MOVIE = "--movie" in sys.argv

files = sorted(glob.glob(os.path.join(base, VTK, "VTK_cells*.vtu")))
assert files, f"no VTK_cells*.vtu in {VTK}"
print(f"{len(files)} frames from {VTK}")

cells = XMLUnstructuredGridReader(registrationName="Epidermis", FileName=files)
view = GetActiveViewOrCreate("RenderView")
view.OrientationAxesVisibility = 0
view.Background = [0.10, 0.10, 0.12]
view.UseColorPaletteForBackground = 0

# Colour by stress anisotropy and let the glyphs carry the direction. Angle to
# the *local* organ axis would be the more direct readout, but the organ is
# curved, so an angle to any global axis would be wrong - animate_switch.py
# computes the local-axis angle properly and colours by that instead.
shell = Show(cells, view)
shell.Representation = "Surface With Edges"
shell.EdgeColor = [0.20, 0.20, 0.24]
ColorBy(shell, ("CELLS", "cell vector length"))
lut = GetColorTransferFunction("cellvectorlength")
lut.ApplyPreset("Inferno (matplotlib)", True)
lut.RescaleTransferFunction(0.0, 0.9)
shell.SetScalarBarVisibility(view, True)
bar = GetScalarBar(lut, view)
bar.Title = "stress anisotropy"
bar.ComponentTitle = ""
bar.WindowLocation = "Lower Left Corner"
bar.ScalarBarLength = 0.28
bar.TitleColor = bar.LabelColor = [1, 1, 1]

# Predicted CMT axes.
centres = CellCenters(registrationName="Cell centres", Input=cells)
glyph = Glyph(registrationName="CMT axes", Input=centres, GlyphType="Line")
glyph.OrientationArray = ["POINTS", "cell vector"]
glyph.ScaleArray = ["POINTS", "cell vector length"]
glyph.ScaleFactor = 46.0
glyph.GlyphMode = "All Points"
gd = Show(glyph, view)
ColorBy(gd, None)
gd.AmbientColor = gd.DiffuseColor = [1.0, 1.0, 1.0]
gd.LineWidth = 3.0

scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
steps = list(GetTimeKeeper().TimestepValues)
clock = AnnotateTimeFilter(registrationName="Clock", Input=cells)
clock.Scale = T_END / max(1.0, float(len(steps) - 1))
clock.Format = "{time:.1f} h after illumination"
cd = Show(clock, view)
cd.FontSize = 22
cd.Bold = 1
cd.Color = [1, 1, 1]
cd.WindowLocation = "Upper Center"

scene.AnimationTime = steps[-1]
view.ResetCamera()
cam = GetActiveCamera()
cam.Azimuth(25); cam.Elevation(15); cam.Dolly(0.9)
scene.AnimationTime = steps[0]
Render()

SaveState(os.path.join(base, "hook_core.pvsm"))
print("state written: hook_core.pvsm")

view.ViewSize = [1280, 1000]
figs = os.path.join(base, "figures")
os.makedirs(figs, exist_ok=True)
if MOVIE:
    SaveAnimation(os.path.join(figs, "hook_opening.png"), view,
                  ImageResolution=[1280, 1000], FrameWindow=[0, len(steps) - 1])
    print("frames written: figures/hook_opening.*.png")
else:
    for frac, name in ((0.0, "t0"), (0.25, "tA"), (0.5, "tB"), (1.0, "tEnd")):
        scene.AnimationTime = steps[int(frac * (len(steps) - 1))]
        Render()
        SaveScreenshot(os.path.join(figs, f"core_{name}.png"), view)
    print("stills written to figures/")
