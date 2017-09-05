/*Copyright (c) 2017, Andreas Grueneis, Felix Hummel and Alejandro Gallo, all
 * rights reserved.*/
#ifndef MP2_EOM_DAVIDSON_DEFINED
#define MP2_EOM_DAVIDSON_DEFINED

#include <algorithms/Algorithm.hpp>
#include <math/FockVector.hpp>
#include <vector>

namespace cc4s {

  /**
   * \brief Implements the diagonal preconditionar for the davidson method
   * \tparam F It is the field variable to be used, in general it will be
   * complex
   * \tparam V The type of vectors to be used. In our case we will be spanning
   * the space in Singles and doubles excitations operators, so it will be
   * a FockVector.
   */
  template <typename F = complex, typename V = FockVector<F> >
  class Mp2PreConditioner
  {
  public:
    Mp2PreConditioner (
    ){
    };
    ~Mp2PreConditioner (){
    };
    std::vector<V> getInitialBasis(int eigenVectorsCount) const;
    V getCorrection(const complex eigenValue, const V &residuum) const;


  private:
    /* data */
  };

  class Mp2EquationOfMotionDavidson: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(Mp2EquationOfMotionDavidson);
    Mp2EquationOfMotionDavidson(
      std::vector<Argument> const &argumentList
    );
    virtual ~Mp2EquationOfMotionDavidson();

    virtual void run();

  protected:
    static constexpr int DEFAULT_MAX_ITERATIONS = 16;

    template <typename F = double>
    void getCanonicalPerturbationBasis(
        CTF::Tensor<F> &Tai, CTF::Tensor<F> &Tabij, int64_t i
    );

  };
}

#endif

