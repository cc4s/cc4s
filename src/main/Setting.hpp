/* Copyright 2021 cc4s.org
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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
