#include <atrip/Slice.hpp>

namespace atrip {

template <typename F>
void Slice<F>::mark_ready() noexcept {
  info.state = Ready;
  info.recycling = Blank;
}

template <typename F>
bool Slice<F>::is_unwrapped() const noexcept {
  return info.state == Ready || info.state == SelfSufficient;
}

template <typename F>
bool Slice<F>::is_unwrappable() const noexcept {
  return is_unwrapped() || info.state == Recycled || info.state == Dispatched;
}

template <typename F>
inline bool Slice<F>::is_directly_fetchable() const noexcept {
  return info.state == Ready || info.state == Dispatched;
}

template <typename F>
void Slice<F>::free() noexcept {
  info.tuple = {0, 0};
  info.type = Blank;
  info.state = Acceptor;
  info.from = {0, 0};
  info.recycling = Blank;
  data = DataNullPtr;
}

template <typename F>
PartialTuple Slice<F>::subtuple_by_slice(ABCTuple abc, Type slice_type) {
  switch (slice_type) {
  case AB: return {abc[0], abc[1]};
  case BC: return {abc[1], abc[2]};
  case AC: return {abc[0], abc[2]};
  case CB: return {abc[2], abc[1]};
  case BA: return {abc[1], abc[0]};
  case CA: return {abc[2], abc[0]};
  case A: return {abc[0], 0};
  case B: return {abc[1], 0};
  case C: return {abc[2], 0};
  default: throw "Switch statement not exhaustive!";
  }
}

template <typename F>
std::vector<Slice<F> *>
Slice<F>::has_recycled_referencing_to_it(std::vector<Slice<F>> &slices,
                                         Info const &info) {
  std::vector<Slice<F> *> result;

  for (auto &s : slices)
    if (s.info.recycling == info.type && s.info.tuple == info.tuple
        && s.info.state == Recycled)
      result.push_back(&s);

  return result;
}

template <typename F>
Slice<F> &Slice<F>::find_one_by_type(std::vector<Slice<F>> &slices,
                                     Slice<F>::Type type) {
  const auto slice_it =
      std::find_if(slices.begin(), slices.end(), [&type](Slice<F> const &s) {
        return type == s.info.type;
      });
  WITH_CRAZY_DEBUG
  WITH_RANK << "\t__ looking for " << type << "\n";
  if (slice_it == slices.end())
    throw std::domain_error("Slice one by type not found!");
  return *slice_it;
}

template <typename F>
Slice<F> &Slice<F>::find_recycled_source(std::vector<Slice<F>> &slices,
                                         Slice<F>::Info info) {
  const auto slice_it =
      std::find_if(slices.begin(), slices.end(), [&info](Slice<F> const &s) {
        return info.recycling == s.info.type && info.tuple == s.info.tuple
            && State::Recycled != s.info.state;
      });

  WITH_CRAZY_DEBUG
  WITH_RANK << "__slice__:find: recycling source of " << pretty_print(info)
            << "\n";
  if (slice_it == slices.end())
    throw std::domain_error("Recycled source not found: " + pretty_print(info)
                            + " rank: " + pretty_print(Atrip::rank));
  WITH_RANK << "__slice__:find: " << pretty_print(slice_it->info) << "\n";
  return *slice_it;
}

template <typename F>
Slice<F> &Slice<F>::find_type_abc(std::vector<Slice<F>> &slices,
                                  Slice<F>::Type type,
                                  ABCTuple const &abc) {
  const auto tuple = Slice<F>::subtuple_by_slice(abc, type);
  const auto slice_it =
      std::find_if(slices.begin(),
                   slices.end(),
                   [&type, &tuple](Slice<F> const &s) {
                     return type == s.info.type && tuple == s.info.tuple;
                   });
  WITH_CRAZY_DEBUG
  WITH_RANK << "__slice__:find:" << type << " and tuple " << pretty_print(tuple)
            << "\n";
  if (slice_it == slices.end())
    throw std::domain_error("Slice by type not found: " + pretty_print(tuple)
                            + ", " + std::to_string(type)
                            + " rank: " + std::to_string(Atrip::rank));
  return *slice_it;
}

template <typename F>
Slice<F> &Slice<F>::find_by_info(std::vector<Slice<F>> &slices,
                                 Slice<F>::Info const &info) {
  const auto slice_it =
      std::find_if(slices.begin(), slices.end(), [&info](Slice<F> const &s) {
        // TODO: maybe implement comparison in Info struct
        return info.type == s.info.type && info.state == s.info.state
            && info.tuple == s.info.tuple && info.from.rank == s.info.from.rank
            && info.from.source == s.info.from.source;
      });
  WITH_CRAZY_DEBUG
  WITH_RANK << "__slice__:find:looking for " << pretty_print(info) << "\n";
  if (slice_it == slices.end())
    throw std::domain_error("Slice by info not found: " + pretty_print(info));
  return *slice_it;
}

template <typename F>
inline bool Slice<F>::is_free() const noexcept {
  return info.tuple == PartialTuple{0, 0} && info.type == Blank
      && info.state == Acceptor && info.from.rank == 0 && info.from.source == 0
      && info.recycling == Blank && data == DataNullPtr;
}

template <typename F>
inline bool Slice<F>::is_recyclable() const noexcept {
  return (info.state == Dispatched || info.state == Ready
          || info.state == Fetch)
      && has_valid_data_pointer();
}

template <typename F>
inline bool Slice<F>::has_valid_data_pointer() const noexcept {
  return data != DataNullPtr && info.state != Acceptor && info.type != Blank;
}

template <typename F>
void Slice<F>::unwrap_and_mark_ready() {
  if (info.state == Ready) return;
  if (info.state != Dispatched)
    throw std::domain_error("Can't unwrap a non-ready, non-dispatched slice!");
  mark_ready();
  MPI_Status status;
#ifdef HAVE_OCD
  WITH_RANK << "__slice__:mpi: waiting "
            << "\n";
#endif
  const int error_code = MPI_Wait(&request, &status);

  // FIXME: it appears not to work to free
  // this request, investigate if this is necessary or not
  // const auto _mpi_request_free = MPI_Request_free(&request);

  // if (MPI_SUCCESS != _mpi_request_free)
  // throw "Atrip: Error freeing MPI request";

  if (error_code != MPI_SUCCESS) throw "Atrip: Unexpected error MPI ERROR";

#if defined(HAVE_CUDA) && !defined(ATRIP_SOURCES_IN_GPU)
  // copy the retrieved mpi data to the device
  WITH_CHRONO("cuda:memcpy",
              _CHECK_CUDA_SUCCESS(
                  "copying mpi data to device",
                  cuMemcpyHtoD(data, (void *)mpi_data, sizeof(F) * size));)
  std::free(mpi_data);
#endif

#ifdef HAVE_OCD
  char error_string[MPI_MAX_ERROR_STRING];
  int error_size;
  MPI_Error_string(error_code, error_string, &error_size);

  WITH_RANK << "__slice__:mpi: status "
            << "{ .source=" << status.MPI_SOURCE << ", .tag=" << status.MPI_TAG
            << ", .error=" << status.MPI_ERROR << ", .errCode=" << error_code
            << ", .err=" << error_string << " }"
            << "\n";
#endif
}

template <typename F>
Slice<F>::Slice(size_t size_)
    : info({})
    , data(DataNullPtr)
#if defined(HAVE_CUDA) && !defined(ATRIP_SOURCES_IN_GPU)
    , mpi_data(nullptr)
#endif
    , size(size_) {
}

// utility functions

template <typename F>
std::string type_to_string(typename Slice<F>::Type t) {
  switch (t) {
  case Slice<F>::AB: return "AB";
  case Slice<F>::BC: return "BC";
  case Slice<F>::AC: return "AC";
  case Slice<F>::CB: return "CB";
  case Slice<F>::BA: return "BA";
  case Slice<F>::CA: return "CA";
  case Slice<F>::A: return "A";
  case Slice<F>::B: return "B";
  case Slice<F>::C: return "C";
  case Slice<F>::Blank: return "Blank";
  default: throw "Switch statement not exhaustive!";
  }
}

template <typename F>
std::string name_to_string(typename Slice<F>::Name t) {
  switch (t) {
  case Slice<F>::TA: return "TA";
  case Slice<F>::VIJKA: return "VIJKA";
  case Slice<F>::VABCI: return "VABCI";
  case Slice<F>::TABIJ: return "TABIJ";
  case Slice<F>::VABIJ: return "VABIJ";
  case Slice<F>::JIJKA: return "JIJKA";
  case Slice<F>::JABCI: return "JABCI";
  default: throw "Switch statement not exhaustive!";
  }
}

template <typename F>
size_t name_to_size(typename Slice<F>::Name t, size_t No, size_t Nv) {
  switch (t) {
  case Slice<F>::TA: return Nv * No * No;
  case Slice<F>::JIJKA:
  case Slice<F>::VIJKA: return No * No * No;
  case Slice<F>::JABCI:
  case Slice<F>::VABCI: return Nv * No;
  case Slice<F>::TABIJ: return No * No;
  case Slice<F>::VABIJ: return No * No;
  default: throw "Switch statement not exhaustive!";
  }
}

template <typename F>
std::string state_to_string(typename Slice<F>::State t) {
  switch (t) {
  case Slice<F>::Fetch: return "Fetch";
  case Slice<F>::Dispatched: return "Dispatched";
  case Slice<F>::Ready: return "Ready";
  case Slice<F>::SelfSufficient: return "SelfSufficient";
  case Slice<F>::Recycled: return "Recycled";
  case Slice<F>::Acceptor: return "Acceptor";
  default: throw "Switch statement not exhaustive!";
  }
}

template std::string state_to_string<double>(typename Slice<double>::State t);
template std::string state_to_string<Complex>(typename Slice<Complex>::State t);
template size_t
name_to_size<double>(typename Slice<double>::Name t, size_t No, size_t Nv);
template size_t
name_to_size<Complex>(typename Slice<Complex>::Name t, size_t No, size_t Nv);
template std::string name_to_string<double>(typename Slice<double>::Name t);
template std::string name_to_string<Complex>(typename Slice<Complex>::Name t);
template std::string type_to_string<double>(typename Slice<double>::Type t);
template std::string type_to_string<Complex>(typename Slice<Complex>::Type t);

template class Slice<double>;
template class Slice<Complex>;

} // namespace atrip
