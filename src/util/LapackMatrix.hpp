/*Copyright (c) 2017, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef LAPACK_MATRIX_DEFINED
#define LAPACK_MATRIX_DEFINED

#include <vector>

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

