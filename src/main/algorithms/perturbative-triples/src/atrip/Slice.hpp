// Copyright 2022 Alejandro Gallo
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// [[file:~/cuda/atrip/atrip.org::*Prolog][Prolog:1]]
#pragma once
#include <iostream>
#include <algorithm>
#include <vector>

#include <atrip/mpi.hpp>
#include <atrip/Types.hpp>
#include <atrip/Tuples.hpp>
#include <atrip/Utils.hpp>
#include <atrip/CUDA.hpp>

namespace atrip {

template <typename F = double>
class Slice {
public:
  // Prolog:1 ends here

  // [[file:~/cuda/atrip/atrip.org::*Location][Location:1]]
  struct Location {
    size_t rank;
    size_t source;
  };
  // Location:1 ends here

  // [[file:~/cuda/atrip/atrip.org::*Type][Type:1]]
  enum Type {
    A = 10,
    B,
    C,

    // Two-parameter slices
    AB = 20,
    BC,
    AC,

    // for abci and the doubles
    CB,
    BA,
    CA,

    // The non-typed slice
    Blank = 404
  };
  // Type:1 ends here

  // [[file:~/cuda/atrip/atrip.org::*State][State:1]]
  enum State {
    Fetch = 0,
    Dispatched = 2,
    Ready = 1,
    SelfSufficient = 911,
    Recycled = 123,
    Acceptor = 405
  };
  // State:1 ends here

  // [[file:~/cuda/atrip/atrip.org::*The%20Info%20structure][The Info
  // structure:1]]
  struct Info {
    // which part of a,b,c the slice holds
    PartialTuple tuple;
    // The type of slice for the user to retrieve the correct one
    Type type;
    // What is the state of the slice
    State state;
    // Where the slice is to be retrieved
    Location from;
    // If the data are actually to be found in this other slice
    Type recycling;

    Info()
        : tuple{0, 0}
        , type{Blank}
        , state{Acceptor}
        , from{0, 0}
        , recycling{Blank} {}
  };

  using Ty_x_Tu = std::pair<Type, PartialTuple>;
  // The Info structure:1 ends here

  // [[file:~/cuda/atrip/atrip.org::*Name][Name:1]]
  enum Name {
    TA = 100,
    VIJKA = 101,
    VABCI = 200,
    TABIJ = 201,
    VABIJ = 202,
    // ct tensor names
    JIJKA = 300,
    JABCI = 301
  };
  // Name:1 ends here

  // [[file:~/cuda/atrip/atrip.org::*Database][Database:1]]
  struct LocalDatabaseElement {
    Slice<F>::Name name;
    Slice<F>::Info info;
  };
  // Database:1 ends here

  // [[file:~/cuda/atrip/atrip.org::*Database][Database:2]]
  using LocalDatabase = std::vector<LocalDatabaseElement>;
  using Database = LocalDatabase;
  // Database:2 ends here

  // [[file:~/cuda/atrip/atrip.org::*MPI%20Types][MPI Types:1]]
  struct mpi {

    static MPI_Datatype vector(size_t n, MPI_Datatype const &DT) {
      MPI_Datatype dt;
      MPI_Type_vector(n, 1, 1, DT, &dt);
      MPI_Type_commit(&dt);
      return dt;
    }

    static MPI_Datatype slice_location() {
      constexpr int n = 2;
      // create a slice_location to measure in the current architecture
      // the packing of the struct
      Slice<F>::Location measure;
      MPI_Datatype dt;
      const std::vector<int> lengths(n, 1);
      const MPI_Datatype types[n] = {usize_dt(), usize_dt()};

      static_assert(sizeof(Slice<F>::Location) == 2 * sizeof(size_t),
                    "The Location packing is wrong in your compiler");

      // measure the displacements in the struct
      size_t j = 0;
      MPI_Aint base_address, displacements[n];
      MPI_Get_address(&measure, &base_address);
      MPI_Get_address(&measure.rank, &displacements[j++]);
      MPI_Get_address(&measure.source, &displacements[j++]);
      for (size_t i = 0; i < n; i++)
        displacements[i] = MPI_Aint_diff(displacements[i], base_address);

      MPI_Type_create_struct(n, lengths.data(), displacements, types, &dt);
      MPI_Type_commit(&dt);
      return dt;
    }

    static MPI_Datatype usize_dt() { return MPI_UINT64_T; }

// airmler: The following lines contain a memory leak.
//         However, they are currently not used.
/*
    static MPI_Datatype slice_info() {
      constexpr int n = 5;
      MPI_Datatype dt;
      Slice<F>::Info measure;
      const std::vector<int> lengths(n, 1);
      const MPI_Datatype types[n] = {
          vector(2, usize_dt()),
          vector(sizeof(enum Type), MPI_CHAR),
          vector(sizeof(enum State), MPI_CHAR),
          slice_location(),
          vector(sizeof(enum Type), MPI_CHAR)
          // TODO: Why this does not work on intel mpi?
          //, MPI_UINT64_T //
      };

      static_assert(sizeof(enum Type) == 4, "Enum type not 4 bytes long");
      static_assert(sizeof(enum State) == 4, "Enum State not 4 bytes long");
      static_assert(sizeof(enum Name) == 4, "Enum Name not 4 bytes long");

      // create the displacements from the info measurement struct
      size_t j = 0;
      MPI_Aint base_address, displacements[n];
      MPI_Get_address(&measure, &base_address);
      MPI_Get_address(&measure.tuple[0], &displacements[j++]);
      MPI_Get_address(&measure.type, &displacements[j++]);
      MPI_Get_address(&measure.state, &displacements[j++]);
      MPI_Get_address(&measure.from, &displacements[j++]);
      MPI_Get_address(&measure.recycling, &displacements[j++]);
      for (size_t i = 0; i < n; i++)
        displacements[i] = MPI_Aint_diff(displacements[i], base_address);

      MPI_Type_create_struct(n, lengths.data(), displacements, types, &dt);
      MPI_Type_commit(&dt);
      return dt;
    }
*/
    static MPI_Datatype local_database_element() {
      constexpr int n = 2;
      MPI_Datatype dt;
      LocalDatabaseElement measure;
      const std::vector<int> lengths(n, 1);
// airmler: remove memory leaks, unused anyways
//      const MPI_Datatype types[n] = {vector(sizeof(enum Name), MPI_CHAR),
//                                     slice_info()};

      // measure the displacements in the struct
      size_t j = 0;
      MPI_Aint base_address, displacements[n];
      MPI_Get_address(&measure, &base_address);
      MPI_Get_address(&measure.name, &displacements[j++]);
      MPI_Get_address(&measure.info, &displacements[j++]);
      for (size_t i = 0; i < n; i++)
        displacements[i] = MPI_Aint_diff(displacements[i], base_address);

      static_assert(sizeof(LocalDatabaseElement) == sizeof(measure),
                    "Measure has bad size");

//      MPI_Type_create_struct(n, lengths.data(), displacements, types, &dt);
//      MPI_Type_commit(&dt);
      return vector(sizeof(LocalDatabaseElement), MPI_CHAR);
      // TODO: write tests in order to know if this works
//      return dt;
    }
  };
  // MPI Types:1 ends here

  // [[file:~/cuda/atrip/atrip.org::*Static%20utilities][Static utilities:1]]
  static PartialTuple subtuple_by_slice(ABCTuple abc, Type slice_type);
  // Static utilities:1 ends here

  // [[file:~/cuda/atrip/atrip.org::*Static%20utilities][Static utilities:2]]
  static std::vector<Slice<F> *>
  has_recycled_referencing_to_it(std::vector<Slice<F>> &slices,
                                 Info const &info);
  // Static utilities:2 ends here

  // [[file:~/cuda/atrip/atrip.org::*Static%20utilities][Static utilities:3]]
  static Slice<F> &find_one_by_type(std::vector<Slice<F>> &slices,
                                    Slice<F>::Type type);
  // Static utilities:3 ends here

  // [[file:~/cuda/atrip/atrip.org::*Static%20utilities][Static utilities:4]]
  static Slice<F> &find_recycled_source(std::vector<Slice<F>> &slices,
                                        Slice<F>::Info info);
  // Static utilities:4 ends here

  // [[file:~/cuda/atrip/atrip.org::*Static%20utilities][Static utilities:5]]
  static Slice<F> &find_type_abc(std::vector<Slice<F>> &slices,
                                 Slice<F>::Type type,
                                 ABCTuple const &abc);
  // Static utilities:5 ends here

  // [[file:~/cuda/atrip/atrip.org::*Static%20utilities][Static utilities:6]]
  static Slice<F> &find_by_info(std::vector<Slice<F>> &slices,
                                Slice<F>::Info const &info);
  // Static utilities:6 ends here

  // [[file:~/cuda/atrip/atrip.org::*Attributes][Attributes:1]]
  Info info;
  // Attributes:1 ends here

  // [[file:~/cuda/atrip/atrip.org::*Attributes][Attributes:2]]
  DataPtr<F> data;
#if defined(HAVE_CUDA) && !defined(ATRIP_SOURCES_IN_GPU)
  F *mpi_data;
#endif
  // Attributes:2 ends here

  // [[file:~/cuda/atrip/atrip.org::*Attributes][Attributes:3]]
  MPI_Request request;
  // Attributes:3 ends here

  // [[file:~/cuda/atrip/atrip.org::*Attributes][Attributes:4]]
  const size_t size;
  // Attributes:4 ends here

  // [[file:~/cuda/atrip/atrip.org::*Member%20functions][Member functions:1]]
  void mark_ready() noexcept;
  // Member functions:1 ends here

  // [[file:~/cuda/atrip/atrip.org::*Member%20functions][Member functions:2]]
  bool is_unwrapped() const noexcept;
  // Member functions:2 ends here

  // [[file:~/cuda/atrip/atrip.org::*Member%20functions][Member functions:3]]
  bool is_unwrappable() const noexcept;

  inline bool is_directly_fetchable() const noexcept;

  void free() noexcept;

  bool is_free() const noexcept;
  // Member functions:3 ends here

  // [[file:~/cuda/atrip/atrip.org::*Member%20functions][Member functions:4]]
  bool is_recyclable() const noexcept;
  // Member functions:4 ends here

  // [[file:~/cuda/atrip/atrip.org::*Member%20functions][Member functions:5]]
  inline bool has_valid_data_pointer() const noexcept;
  // Member functions:5 ends here

  // [[file:~/cuda/atrip/atrip.org::*Member%20functions][Member functions:6]]
  void unwrap_and_mark_ready();

  // CONSTRUCTOR
  Slice(size_t size_);

}; // struct Slice
// Epilog:1 ends here

// [[file:~/cuda/atrip/atrip.org::*Debug][Debug:1]]
template <typename F = double>
std::ostream &operator<<(std::ostream &out,
                         typename Slice<F>::Location const &v) {
  // TODO: remove me
  out << "{.r(" << v.rank << "), .s(" << v.source << ")};";
  return out;
}

template <typename F = double>
std::ostream &operator<<(std::ostream &out, typename Slice<F>::Info const &i) {
  out << "«" << i.type << ", " << i.state << "»"
      << " ⊙ {" << i.from.rank << ", " << i.from.source << "}"
      << " ∴ {" << i.tuple[0] << ", " << i.tuple[1] << "}"
      << " ♲" << i.recycling;
  return out;
}

template <typename F>
std::string type_to_string(typename Slice<F>::Type t);

template <typename F>
std::string name_to_string(typename Slice<F>::Name t);

template <typename F>
size_t name_to_size(typename Slice<F>::Name t, size_t No, size_t Nv);

template <typename F>
std::string state_to_string(typename Slice<F>::State t);

template <typename F = double>
std::ostream &operator<<(std::ostream &out, typename Slice<F>::State const &i) {
  out << state_to_string<F>(i);
  return out;
}

template <typename F = double>
std::ostream &operator<<(std::ostream &out, typename Slice<F>::Name const &i) {
  out << name_to_string<F>(i);
  return out;
}

template <typename F = double>
std::ostream &operator<<(std::ostream &out, typename Slice<F>::Type const &i) {
  out << type_to_string<F>(i);
  return out;
}

} // namespace atrip
// Debug:1 ends here
