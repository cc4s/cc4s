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

#ifndef INTEGRATION_DEFINED
#define INTEGRATION_DEFINED

namespace cc4s {
  template<typename F, typename Real, typename Method>
  Real integrate(F f, Real a, Real b, size_t steps, Method m) {
    Real s(0);
    Real h((b-a)/steps);
    for (size_t i(0); i <= steps; ++i) {
      s += m(f, (a*(steps-i) + b*i)/steps, h);
    }
    return h*s;
  }

  class Trapezium {
  public:
    template<typename F, typename Real>
    Real operator()(F f, Real x, Real h) const {
      return (f(x) + f(x+h))/2;
    }
  };

  class Simpson {
  public:
    template<typename F, typename Real>
    Real operator()(F f, Real x, Real h) const {
      return (f(x) + 4*f(x+h/2) + f(x+h))/6;
    }
  };
}
#endif
