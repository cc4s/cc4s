/* Copyright 2021 cc4s.org
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <algorithms/Algorithm.hpp>
#include <Node.hpp>
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

