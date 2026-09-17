//
// Minimal command-line/config registry mirroring the legacy myConfig:
// options registered with an argument count, parsed from argv (with optional
// ~/.tissue config file), remaining tokens become positional arguments.
//
#ifndef TISSUE2_IO_CONFIG_H
#define TISSUE2_IO_CONFIG_H

#include <string>

namespace tissue::config {

void registerOption(const std::string &name, int numArguments);
// Parses argv and (if present) the config file; unknown -options are fatal.
void init(int argc, char *argv[], const std::string &configFile);

// Value of an option's i-th argument; empty if unset.
std::string getValue(const std::string &name, size_t i);
bool getBooleanValue(const std::string &name);

// Positional (non-option) arguments; argv(0) is the program name.
int argc();
std::string argv(int i);

} // namespace tissue::config

#endif
