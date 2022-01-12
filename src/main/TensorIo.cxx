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

#include <TensorIo.hpp>
#include <Reader.hpp>

// using namespace cc4s;

int cc4s::TensorIo::WRITE_REGISTERED =
  cc4s::Writer::registerWriteFunction("tensor", cc4s::TensorIo::write);

int cc4s::TensorIo::READ_REGISTERED =
  cc4s::Reader::registerReadFunction("tensor", cc4s::TensorIo::read);

