/*Copyright (c) 2018, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef LAPACK_MATRIX_DEFINED
#define LAPACK_MATRIX_DEFINED

#include <vector>
#include <sstream>
#include <ctf.hpp>

namespace cc4s {
  template <typename F=double>
  class LapackMatrix {
  public:
    /**
     * \brief Copies the content of the Lapack matrix
     **/
    LapackMatrix(
      const LapackMatrix<F> &A
    ):
      rows(A.rows), columns(A.columns), values(A.values)
    {
    }

    /**
     * \brief Construct a zero nxm Lapack.
     **/
    LapackMatrix(
      const int rows_, const int columns_
    ):
      rows(rows_), columns(columns_), values(rows_*columns_)
    {
    }

    /**
     * \brief Constructs an LapackMatrix from a CTF tensor on all ranks
     **/
    LapackMatrix(
      CTF::Tensor<F> &ctfA
    ):
      rows(ctfA.lens[0]), columns(ctfA.lens[1])
    {
      values.resize(ctfA.lens[0]*ctfA.lens[1]);
      ctfA.read_all(values.data(), true);
    }

    LapackMatrix<F> &operator =(const LapackMatrix<F> &A) {
      rows = A.rows;
      columns = A.columns;
      values = A.values;
      return *this;
    }

    const F &operator ()(const  int i, const int j) const {
      return values[i+j*rows];
    }

    F &operator ()(const int i, const int j) {
      return values[i+j*rows];
    }

    int getRows() const { return rows; }
    int getColumns() const { return columns; }

    /**
     * \brief Returns the pointer to the column major data.
     **/
    const F *getValues() const {
      return values.data();
    }

    /**
     * \brief Returns the pointer to the mutable column major data.
     **/
    F *getValues() {
      return values.data();
    }

    /**
     * \brief Writes the data of this Lapack matrix to the CTF tensor.
     **/
    void write(CTF::Tensor<F> &ctfA) const {
      if (ctfA.lens[0] != rows || ctfA.lens[1] != columns) {
        std::stringstream stream;
        stream << "Tensor is not of correct shape to receive ("
          << rows << "x" << columns << ") matrix: ";
        std::string join("");
        for (int d(0); d < ctfA.order; ++d) {
          stream << join << ctfA.lens[d];
          join = "x";
        }
        throw new EXCEPTION(stream.str());
      }
      int64_t size(values.size());
      int64_t localSize(ctfA.wrld->rank == 0 ? size : 0);
      std::vector<int64_t> indices(localSize);
      for (int64_t i(0); i < localSize; ++i) { indices[i] = i; }
      ctfA.write(localSize, indices.data(), values.data());
    }

  protected:
    /**
     * \brief Number of rows and columns.
     */
    int rows, columns;

    /**
     * \brief The column major data.
     */
    std::vector<F> values;
  };
}

#endif

