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
      CTF::Tensor<F> &Fij,
      CTF::Tensor<F> &Fab,
      CTF::Tensor<F> &Tabij,
      CTF::Tensor<F> &Vabcd,
      CTF::Tensor<F> &Viajb,
      CTF::Tensor<F> &Vijab,
      CTF::Tensor<F> &Vijkl
    ): Fij(Fij),
      Fab(Fab),
      Tabij(Tabij),
      Vabcd(Vabcd),
      Viajb(Viajb),
      Vijab(Vijab),
      Vijkl(Vijkl)
    { };
    ~Mp2PreConditioner (){
    };
    std::vector<V> getInitialBasis(int eigenVectorsCount) const;
    V getCorrection(const complex eigenValue, const V &residuum) const;

  private:

    CTF::Tensor<F> Fij, Fab, Tabij, Vabcd, Viajb, Vijab, Vijkl;

    F computeDiagonalElement(
      CTF::Tensor<F> &Rai,
      CTF::Tensor<F> &Lia,
      CTF::Tensor<F> &Rabij,
      CTF::Tensor<F> &Lijab
    )
    {
      CTF::Scalar<F> energy;
      energy[""]  = ( - 1.0  ) * Rai["aj"] * Fij["jk"] * Lia["ka"];
      energy[""] += ( + 1.0  ) * Rai["aj"] * Fab["ca"] * Lia["jc"];
      energy[""] += ( - 1.0  ) * Rai["aj"] * Viajb["jcla"] * Lia["lc"];
      energy[""] += ( + 1.0  ) * Rai["aj"] * Tabij["cdmn"] * Vijab["njda"] * Lia["mc"];
      energy[""] += ( + 0.5  ) * Rai["aj"] * Tabij["cdmn"] * Vijab["njcd"] * Lia["ma"];
      energy[""] += ( + 0.5  ) * Rai["aj"] * Tabij["cdmn"] * Vijab["mnda"] * Lia["jc"];
      energy[""] += ( - 1.0  ) * Rabij["abkl"] * Fij["km"] * Lijab["mlab"];
      energy[""] += ( + 1.0  ) * Rabij["abkl"] * Fij["lm"] * Lijab["mkab"];
      energy[""] += ( - 1.0  ) * Rabij["abkl"] * Fab["eb"] * Lijab["klea"];
      energy[""] += ( + 1.0  ) * Rabij["abkl"] * Fab["ea"] * Lijab["kleb"];
      energy[""] += ( - 0.5  ) * Rabij["abkl"] * Vijkl["klmn"] * Lijab["nmab"];
      energy[""] += ( + 1.0  ) * Rabij["abkl"] * Viajb["kenb"] * Lijab["nlea"];
      energy[""] += ( - 1.0  ) * Rabij["abkl"] * Viajb["kena"] * Lijab["nleb"];
      energy[""] += ( - 1.0  ) * Rabij["abkl"] * Viajb["lenb"] * Lijab["nkea"];
      energy[""] += ( + 1.0  ) * Rabij["abkl"] * Viajb["lena"] * Lijab["nkeb"];
      energy[""] += ( - 0.5  ) * Rabij["abkl"] * Vabcd["efab"] * Lijab["klfe"];
      energy[""] += ( + 0.5  ) * Rabij["abkl"] * Tabij["efop"] * Vijab["klfb"] * Lijab["poea"];
      energy[""] += ( - 0.5  ) * Rabij["abkl"] * Tabij["efop"] * Vijab["klfa"] * Lijab["poeb"];
      energy[""] += ( - 0.25  ) * Rabij["abkl"] * Tabij["efop"] * Vijab["klef"] * Lijab["poab"];
      energy[""] += ( - 0.5  ) * Rabij["abkl"] * Tabij["efop"] * Vijab["pkab"] * Lijab["olfe"];
      energy[""] += ( + 0.5  ) * Rabij["abkl"] * Tabij["efop"] * Vijab["plab"] * Lijab["okfe"];
      energy[""] += ( - 1.0  ) * Rabij["abkl"] * Tabij["efop"] * Vijab["pkfb"] * Lijab["olea"];
      energy[""] += ( + 1.0  ) * Rabij["abkl"] * Tabij["efop"] * Vijab["pkfa"] * Lijab["oleb"];
      energy[""] += ( + 1.0  ) * Rabij["abkl"] * Tabij["efop"] * Vijab["plfb"] * Lijab["okea"];
      energy[""] += ( - 1.0  ) * Rabij["abkl"] * Tabij["efop"] * Vijab["plfa"] * Lijab["okeb"];
      energy[""] += ( + 0.5  ) * Rabij["abkl"] * Tabij["efop"] * Vijab["pkef"] * Lijab["olab"];
      energy[""] += ( - 0.5  ) * Rabij["abkl"] * Tabij["efop"] * Vijab["plef"] * Lijab["okab"];
      energy[""] += ( - 0.25  ) * Rabij["abkl"] * Tabij["efop"] * Vijab["opab"] * Lijab["klfe"];
      energy[""] += ( - 0.5  ) * Rabij["abkl"] * Tabij["efop"] * Vijab["opfb"] * Lijab["klea"];
      energy[""] += ( + 0.5  ) * Rabij["abkl"] * Tabij["efop"] * Vijab["opfa"] * Lijab["kleb"];
      return energy.get_val();
    }
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

