#ifndef SYMMETRIZE_DEFINED
#define SYMMETRIZE_DEFINED
  /**
   * \brief Calculate the sign of a permutation of strings, e.g.
   * sign("abcd", "bacd") = -1
   */
  int permutationSign(const std::string &original, const std::string &permuted);

  /**
   * \brief Apply a permutation operator and antisymmetrize accordingly, e.g.
   *  antiSymmetrize(X, "abcd", "abdc") is replaced by
   *  -1 X["abcd"] = X["abcd"] + sign(abcd -> abdc) X["abdc"]
   *
   */
  template <typename F>
  inline void antiSymmetrize(std::string indices, std::string permuted,
      CTF::Tensor<F> &t, F prefactor=1) {
    double sign(permutationSign(indices, permuted));
    t[indices.c_str()] += prefactor * sign * t[permuted.c_str()];
  }
  /**
   * \brief Apply a permutation operator and antisymmetrize accordingly, e.g.
   *  antiSymmetrize(X, "abcd", "abdc") is replaced by
   *  -1 X["abcd"] = X["abcd"] + sign(abcd -> abdc) X["abdc"]
   *
   */
  template <typename F>
  inline void symmetrize(std::string indices, std::string permuted,
      CTF::Tensor<F> &t, F prefactor=1) {
    t[indices.c_str()] += prefactor * t[permuted.c_str()];
  }
  template <typename F>
  inline void checkAntisymmetry(CTF::Tensor<F> &t){
    CTF::Tensor<F> testResultUp(t);
    CTF::Tensor<F> testResultDown(t);
    F normValue;
    testResultUp["abij"] += testResultUp["baij"];
    testResultDown["abij"] += testResultDown["abji"];
    normValue = testResultUp.norm1();
    normValue += testResultDown.norm1();
    if (normValue >= 1e-3) {
      t.print();
      LOG(0, "AntisymmetryCheck") << t.get_name()
        << ": zero tensor norm " << normValue << std::endl;
      exit(1);
    }
  }
#endif
