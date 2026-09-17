#ifndef TISSUE2_IO_FILE_UTILS_H
#define TISSUE2_IO_FILE_UTILS_H

#include <memory>
#include <sstream>
#include <string>

namespace tissue::io {

// Opens a file and strips '#'-to-end-of-line comments (anywhere on a line),
// returning a stream over the remaining whitespace-delimited tokens — the
// legacy myFiles::openFile behavior used for model, init and solver files.
// Returns nullptr if the file cannot be opened.
std::unique_ptr<std::istringstream> openCommentFiltered(const std::string &file);

} // namespace tissue::io

#endif
