#include "tissue/io/config.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

namespace tissue::config {

namespace {
struct Option {
  int numArguments = 0;
  bool set = false;
  std::vector<std::string> values;
};
std::map<std::string, Option> &options() {
  static std::map<std::string, Option> opts;
  return opts;
}
std::vector<std::string> &positionals() {
  static std::vector<std::string> args;
  return args;
}

void parseTokens(const std::vector<std::string> &tokens) {
  for (size_t i = 0; i < tokens.size(); ++i) {
    const std::string &tok = tokens[i];
    if (!tok.empty() && tok[0] == '-' && tok.size() > 1 &&
        !(std::isdigit(static_cast<unsigned char>(tok[1])) || tok[1] == '.')) {
      std::string name = tok.substr(1);
      auto it = options().find(name);
      if (it == options().end()) {
        std::cerr << "Unknown option '" << tok << "'. Try -help." << std::endl;
        std::exit(EXIT_FAILURE);
      }
      Option &opt = it->second;
      opt.set = true;
      opt.values.clear();
      for (int a = 0; a < opt.numArguments; ++a) {
        if (++i >= tokens.size()) {
          std::cerr << "Option -" << name << " expects " << opt.numArguments
                    << " argument(s)." << std::endl;
          std::exit(EXIT_FAILURE);
        }
        opt.values.push_back(tokens[i]);
      }
    } else {
      positionals().push_back(tok);
    }
  }
}
} // namespace

void registerOption(const std::string &name, int numArguments) {
  options()[name].numArguments = numArguments;
}

void init(int argc, char *argv[], const std::string &configFile) {
  positionals().push_back(argv[0]);
  // Config file first (command line overrides).
  std::ifstream cfg(configFile);
  if (cfg) {
    std::vector<std::string> tokens;
    std::string tok;
    while (cfg >> tok)
      tokens.push_back(tok);
    parseTokens(tokens);
  }
  std::vector<std::string> tokens(argv + 1, argv + argc);
  parseTokens(tokens);
}

std::string getValue(const std::string &name, size_t i) {
  auto it = options().find(name);
  if (it == options().end() || !it->second.set ||
      i >= it->second.values.size())
    return "";
  return it->second.values[i];
}

bool getBooleanValue(const std::string &name) {
  auto it = options().find(name);
  return it != options().end() && it->second.set;
}

int argc() { return static_cast<int>(positionals().size()); }

std::string argv(int i) { return positionals()[static_cast<size_t>(i)]; }

} // namespace tissue::config
