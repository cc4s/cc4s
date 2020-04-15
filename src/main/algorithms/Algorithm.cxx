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

Ptr<AlgorithmFactory::AlgorithmMap> AlgorithmFactory::algorithmMap;

