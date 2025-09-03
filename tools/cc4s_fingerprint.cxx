#include <array>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <functional>

#include "picosha2.h"

#ifndef FINGERPRINT_DEFAULT_PARAM_CHAR
#define FINGERPRINT_DEFAULT_PARAM_CHAR '0'
#endif

// ---------- Hash function ----------
std::string sha256sum(const std::string& data) {
  std::vector<unsigned char> hash(picosha2::k_digest_size);
  picosha2::hash256(data.begin(), data.end(), hash.begin(), hash.end());

  return picosha2::bytes_to_hex_string(hash.begin(), hash.end());
}

// ---------- Hash wrappers ----------
std::string default_hash_func(const std::string& data) { return sha256sum(data); }
std::string new_hash_func(const std::string& data)     { return sha256sum(data); }

// ---------- Parameters ----------
struct FingerprintParams {
  size_t chunks;
  size_t bchunk;
  size_t payload_len;
  std::function<std::string(const std::string&)> hash_func;
};

FingerprintParams get_params(char param_char) {
  switch (param_char) {
    case '0': return { 12, 64, 10, default_hash_func };
    case '1': return { 16, 512, 12, new_hash_func };
    default:
      std::cerr << "Warning: Unknown parameter version " << param_char << ". Using default.\n";
      return { 12, 64, 8, default_hash_func };
  }
}

constexpr char BASE36[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";

// ---------- Utilities ----------
size_t filesize(const std::string& filename) {
  std::ifstream f(filename, std::ios::binary | std::ios::ate);
  return f ? static_cast<size_t>(f.tellg()) : 0;
}

std::vector<char> read_chunk(const std::string& filename, size_t offset, size_t len) {
  std::vector<char> buf(len);
  std::ifstream f(filename, std::ios::binary);
  if (!f) return {};
  f.seekg(offset);
  f.read(buf.data(), len);
  buf.resize(f.gcount());
  return buf;
}

std::string checksum_base36(const std::string& s) {
  int sum = 0;
  for (char ch : s) {
    char u = std::toupper(static_cast<unsigned char>(ch));
    int v = 0;
    if (u >= '0' && u <= '9') v = u - '0';
    else if (u >= 'A' && u <= 'Z') v = u - 'A' + 10;
    sum = (sum + v) % 36;
  }
  return std::string(1, BASE36[sum]);
}

std::string base64_alnum(const std::string& s) {
  std::string out;
  for (char c : s) {
    if (c == '+') out += 'A';
    else if (c == '/') out += 'Z';
    else if (c != '=') out += c;
  }
  return out;
}

// ---------- Fingerprint ----------
std::string fingerprint_file(const std::string& filename, char param_char) {
  auto p = get_params(param_char);
  size_t size = filesize(filename);
  if (size == 0) {
    std::cerr << "Warning: File size is zero. Hash could not be built for " << filename << "\n";
    return "";
  }

  // collect data
  // 1.) get the size of the file
  std::string data;
  std::ostringstream oss;
  oss << std::hex << std::setw(16) << std::setfill('0') << size;
  data += oss.str();

  // 2.) get the different chunks of data from the file
  if (size > 0) {
    size_t maxpos = (size > p.bchunk ? (size - p.bchunk) : 0);
    for (size_t i = 0; i < p.chunks; ++i) {
      size_t pos = 0;
      if (p.chunks > 1) {
          pos = i * maxpos / (p.chunks - 1);
      }
      std::vector<char> chunk = read_chunk(filename, pos, p.bchunk);
      data.append(chunk.data(), chunk.size());
    }
  }

  std::string hash = p.hash_func(data);

  std::string payload = hash.substr(0, p.payload_len);
  std::string check = checksum_base36(payload);

  return std::string(1, param_char) + payload + check;
}

// ---------- Main ----------
int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " FILE\n";
    return 1;
  }
  char param_char = FINGERPRINT_DEFAULT_PARAM_CHAR;
  std::cout << fingerprint_file(argv[1], param_char) << "\n";
  return 0;
}
