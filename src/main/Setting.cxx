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

#include <Setting.hpp>
#include <regex>

using namespace cc4s;


std::vector<SettingName> cc4s::parseDotNotation(const std::string &l) {
  const std::regex rx("([^.]+)");
  std::smatch match;
  std::string lcopy(l);
  std::vector<std::string> result;
  while(std::regex_search(lcopy, match, rx)) {
    result.push_back(match[1].str());
    lcopy = match.suffix();
  }
  return result;
}
