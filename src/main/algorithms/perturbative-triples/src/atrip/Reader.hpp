#ifndef READER_HPP_
#define READER_HPP_

#include <string>
#include <atrip/Atrip.hpp>
#include <atrip/RankMap.hpp>
#include <atrip/Complex.hpp>

namespace atrip {

template <typename F>
Sources<F> diskReader(const std::string &file_path,
                      std::vector<size_t> tensor_dimension,
                      std::vector<size_t> slice_mapping,
                      bool rowMajor);

#if defined(HAVE_CTF)
/* Note: in the current version the ctf tensor data can only be
 *       deleted when calling directly the ctfReader (see below
 *       the implementation of the reader). This should be only
 *       a minor constraint because when we call atrip as a stand
 *       alone program it is not wise to use the ctf reader at all
 */
template <typename F>
Sources<F> ctfReader_fallback(CTF::Tensor<F>& tensor,
                              std::vector<size_t> tensor_dimension,
                              std::vector<size_t> slice_mapping,
                              bool delete_tensor_data = false);

template <typename F>
Sources<F> ctfReader(CTF::Tensor<F>& tensor,
                     std::vector<size_t> tensor_dimension,
                     std::vector<size_t> slice_mapping,
                     bool delete_tensor_data = false);

/*
template <typename F>
Sources<F> vertexReader(std::vector<size_t> tensor_dimension,
                        std::vector<size_t> slice_mapping,
                        CTF::Tensor<F>* hhVertex,
                        CTF::Tensor<F>* phVertex,
                        CTF::Tensor<F>* hpVertex,
                        CTF::Tensor<F>* ppVertex);
*/
#endif

template <typename F>
Sources<F> reader(void* tensor_,
                  std::vector<size_t> tensor_dimension,
                  std::vector<size_t> slice_mapping,
                  const std::string &file_path = "",
                  bool rowMajor = false) {
  if (tensor_ == nullptr) {
    assert(!file_path.empty());
    LOG(0, "Atrip") << "Loading tensor from disk. File path: " << file_path << std::endl;
    return diskReader<F>(file_path,
                         tensor_dimension,
                         slice_mapping,
                         rowMajor);
  }
#if defined(HAVE_CTF)
  auto *tensor = reinterpret_cast<CTF::Tensor<F>*>(tensor_);
  if (tensor == nullptr) {
    throw std::invalid_argument("Invalid tensor pointer, failed dynamic_cast");
  }
  LOG(0, "Atrip") << "Loading tensor from CTF. Tensor name: " << tensor->name << std::endl;
#if defined(CTF_SWITCH_REDISTRIBUTION)
  if (Atrip::useSwitchRedistribution) {
    return ctfReader<F>(*tensor,
                        tensor_dimension,
                        slice_mapping);
  }
#endif
  return ctfReader_fallback<F>(*tensor,
                               tensor_dimension,
                               slice_mapping);


#endif
  if (!Atrip::rank) std::cout << "CTF not available. Abort!" << std::endl;
  assert(0);
  return Sources<F>{};
}


// read file data into a vector
template <typename F>
std::vector<F> read_all(std::vector<size_t> lengths,
                        std::string const &file_path,
                        MPI_Comm comm);


} // namespace atrip

#endif
