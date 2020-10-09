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
