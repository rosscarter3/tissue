#include "tissue/io/file_utils.h"

#include <fstream>

namespace tissue::io {

std::unique_ptr<std::istringstream> openCommentFiltered(const std::string &file) {
  std::ifstream in(file);
  if (!in)
    return nullptr;
  std::ostringstream out;
  std::string line;
  while (std::getline(in, line)) {
    auto pos = line.find('#');
    if (pos != std::string::npos)
      line.erase(pos);
    out << line << '\n';
  }
  return std::make_unique<std::istringstream>(out.str());
}

} // namespace tissue::io
