#ifndef SETTINGS_HEADER_DEFINED
#define SETTINGS_HEADER_DEFINED

#include <string>
#include <map>
#include <vector>

namespace cc4s {

  enum SettingType
    { INT     // int_64
    , NAT     // size_t
    , REAL    // double
    , BOOL
    , STRING
    , TENSOR
    , MAP
    };

  using SettingName = std::string;

  struct Setting {
    const SettingName name;
    const SettingName description;
    const SettingType type;
  };

  using Settings = std::map<SettingName, Setting>;

  std::vector<SettingName> parseDotNotation(const std::string &l);

}

#define SETTING(NAME, DESCRIPTION, TYPE) \
  { NAME, { NAME, DESCRIPTION, SettingType::TYPE }}

#endif
