//
// Self-registering string-keyed factory. Replaces the legacy ~1000-line
// if/else chains in baseReaction.cc, baseCompartmentChange.cc and
// baseSolver.cc: each concrete class registers itself (with any legacy alias
// names) at static-init time.
//
#ifndef TISSUE2_CORE_REGISTRY_H
#define TISSUE2_CORE_REGISTRY_H

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace tissue {

template <typename Base, typename... CtorArgs> class Registry {
public:
  using Creator = std::function<std::unique_ptr<Base>(CtorArgs...)>;

  static Registry &instance() {
    static Registry r;
    return r;
  }

  void add(std::initializer_list<std::string> names, Creator creator) {
    for (const auto &n : names)
      creators_[n] = creator;
  }

  bool contains(const std::string &name) const {
    return creators_.count(name) != 0;
  }

  std::unique_ptr<Base> create(const std::string &name, CtorArgs... args) const {
    auto it = creators_.find(name);
    if (it == creators_.end())
      return nullptr;
    return it->second(std::forward<CtorArgs>(args)...);
  }

  std::vector<std::string> names() const {
    std::vector<std::string> out;
    out.reserve(creators_.size());
    for (const auto &kv : creators_)
      out.push_back(kv.first);
    return out;
  }

private:
  std::map<std::string, Creator> creators_;
};

} // namespace tissue

#endif
