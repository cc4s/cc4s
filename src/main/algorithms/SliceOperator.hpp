/*Copyright (c) 2020, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef SLICE_OPERATOR_DEFINED
#define SLICE_OPERATOR_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class SliceOperator: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(SliceOperator);
    Ptr<MapNode> run(const Ptr<MapNode> &arguments) override;

  protected:
    template <typename F, typename TE>
    Ptr<MapNode> run(const Ptr<MapNode> &arguments);
    template <typename F, typename TE>
    void slice(const Ptr<Tensor<F,TE>> &tensor, const std::string &parts);

    Ptr<MapNode> slices;
    std::vector<size_t> dims;
    size_t No, Nv;
  };
}

#endif

