/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <algorithms/Algorithm.hpp>
#include <Data.hpp>
#include <tcc/Tcc.hpp>
#include <math/Real.hpp>
#include <math/Complex.hpp>
#include <util/Exception.hpp>
#include <util/Log.hpp>

#include <iostream>
#include <sstream>

using namespace cc4s;

// default constructor and destructor
Algorithm::Algorithm() {
}

Algorithm::~Algorithm() {
}

/**
 * \brief The dryRun estimates resource consumption, especially
 * memory and processor time.
 */
Ptr<MapNode> Algorithm::dryRun(const Ptr<MapNode> &arguments) {
  LOG(0, getName()) << "dry run not implemented" << std::endl;
  return arguments;
}

Ptr<AlgorithmFactory::AlgorithmMap> AlgorithmFactory::algorithmMap;

