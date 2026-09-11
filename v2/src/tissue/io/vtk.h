//
// VTK unstructured-grid output (.vtu + .pvd collections), format-compatible
// with the legacy VTUostream/PVD_file writers so existing ParaView pipelines
// keep working. Cells are written shrunk towards their centroid (relative
// wall thickness 0.1) and walls as quads filling the gap.
//
#ifndef TISSUE2_IO_VTK_H
#define TISSUE2_IO_VTK_H

#include <string>

#include "tissue/core/tissue.h"

namespace tissue::io {

// Writes the full .pvd collection referencing n timestep files per part.
void writeFullPvd(const std::string &pvdFile, const std::string &cellFile,
                  const std::string &wallFile, size_t n);

// Writes VTK_cells<count>.vtu / VTK_walls<count>.vtu (6-digit zero padded
// counter before the extension). singleWall selects the legacy write_walls2
// (plain wall variables) vs write_walls3 (paired two-sided wall variables).
void writeVtu(const Tissue &T, const std::string &cellFile,
              const std::string &wallFile, size_t count, bool singleWall);

} // namespace tissue::io

#endif
