A python library to read Tissue output.


It depends on the pyvista library to read meshes with information such as the standard vtk output of Tissue.
```
out = output.TissueOutput.from_dir(OUTPUT_DIR, everyN=10)
```

You can then access the cells and wall meshes along with associated variables:
```
cell_meshes = out.tcells
wall_meshes = out.twalls
```

Cell information can also be exported as a Pandas dataframe for further analysis/visualisation:
```
df = out.to_pandas_cells()
```

It can also read ply files including 'extended'-ply files with geometry+variables meshes:
```
out = output.TissueOutput.from_dir_ply(OUTPUT_DIR)
```



