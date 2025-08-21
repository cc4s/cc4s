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

#ifndef TCC_CONTRACTION_OPERATION_DEFINED
#define TCC_CONTRACTION_OPERATION_DEFINED

#include <tcc/IndexedTensorOperation.hpp>

#include <tcc/Costs.hpp>
#include <tcc/Tensor.hpp>
#include <SharedPointer.hpp>
#include <tcc/Tcc.hpp>
#include <tcc/Tensor_flags.hpp>
#include <Timer.hpp>


#include <string>
#include <unordered_map>
namespace cc4s {
  template <typename F, typename TE> class Contraction;

  template <typename F, typename TE>
  class ContractionOperation: public IndexedTensorOperation<F,TE> {
  public:
    /**
     * \brief Creates a contraction operation contracting the results of
     * the left and the right sub-operations where the result is to be
     * stored in the specified result tensor.
     * Not intended for direct invocation. Use Tcc::compile(expression) to
     * generate operations.
     **/
    ContractionOperation(
      const Ptr<IndexedTensorOperation<F,TE>> &left_,
      const Ptr<IndexedTensorOperation<F,TE>> &right_,
      const Ptr<Tensor<F,TE>> &result_,
      const char *resultIndices_,
      Costs contractionCosts_,
      const std::string &file_, const size_t line_,
      const typename Operation<TE>::ProtectedToken &
    ):
      IndexedTensorOperation<F,TE>(
        result_, resultIndices_,
        contractionCosts_, left_->costs + right_->costs,
        file_, line_, typename Operation<TE>::ProtectedToken()
      ),
      contractionCosts(contractionCosts_),
      left(left_), right(right_)
    {
    }

    void execute() override {
      left->execute();
      right->execute();
      if (
        this->template isOlderThan<F>(left) ||
        this->template isOlderThan<F>(right)
      ) {
        std::string s_contr = (cc4s::is_lowmem_enabled())
                            ? "exexuting: lowmem contraction "
                            : "executing: contraction ";
        LOG_LOCATION(SourceLocation(this->file, this->line)) <<
          s_contr <<
          this->getName() << " <<= " <<
          this->alpha << " * " <<
          left->getName() << " * " << right->getName() << " + " <<
          this->beta << " * " << this->getName() << std::endl;
        auto startTime = Time::getCurrentRealTime();
        if (cc4s::is_lowmem_enabled()) {
          this->getResult()->getMachineTensor()->contract(
            this->alpha,
            left->getResult()->getMachineTensor(), left->getResultIndices(),
            right->getResult()->getMachineTensor(), right->getResultIndices(),
            this->beta,
            this->getResultIndices(),
            0.1
          );
        } else {
          this->getResult()->getMachineTensor()->contract(
            this->alpha,
            left->getResult()->getMachineTensor(), left->getResultIndices(),
            right->getResult()->getMachineTensor(), right->getResultIndices(),
            this->beta,
            this->getResultIndices()
          );
        }
        auto contractTime = Time::getCurrentRealTime() - startTime;
        // To represent the contraction information uniquely,
        // we remove the tensor names and eliminate arbitrary index labeling.
        // Additionally, we include the dimensions (lengths) of each index.
        std::unordered_map<char,char> idx_map;
        char next_label = 'a';
        auto map_index = [&](char c) {
          if (idx_map.find(c) == idx_map.end()) {
            idx_map[c] = next_label++;
          }
          return idx_map[c];
        };
        for (char c :  this->getResultIndices()) map_index(c);
        for (char c :  left->getResultIndices()) map_index(c);
        for (char c : right->getResultIndices()) map_index(c);
        auto remap = [&](const std::string& s) {
            std::string out;
            for (char c : s) out += idx_map[c];
            return out;
        };
        // We need the lens of a tensor in comma-separated string
        auto join = [](const std::vector<size_t>& vec) -> std::string {
          std::ostringstream oss;
          for (size_t i = 0; i < vec.size(); ++i) {
            if (i > 0) oss << ",";
            oss << vec[i];
          }
          return oss.str();
        };
        std::string line("[");
        line += remap( this->getResultIndices()) + ":";
        line += join( this->getResult()->lens) + "] = [";
        line += remap( left->getResultIndices()) + ":";
        line += join( left->getResult()->lens) + "] * [";
        line += remap(right->getResultIndices()) + ":";
        line += join(right->getResult()->lens) + "]";
        /*
         * Actually i dont need the hash function. i leave it.
        // Base62 hash of length 10.
        auto hash10 = [](const std::string& input) -> std::string {
          static const std::string chars
            = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
          uint64_t hash = std::hash<std::string>{}(input);
          std::string result;
          for (int i = 0; i < 10; ++i) {
            result += chars[hash % 62];
            hash /= 62;
          }
          return result;
        };
        */
        this->updated();
        auto _before = this->getFloatingPointOperations();
        this->accountFlops(contractionCosts);
        auto _after = this->getFloatingPointOperations();
        line = "Contraction: " + line + " | "
               + std::to_string((_after - _before) / 1e9) + " GF";
        if (!Cc4s::dryRun) {
          std::ostringstream oss;
          oss << ", " << contractTime << " s";
          line += oss.str();
        }
        LOG_LOCATION(SourceLocation(this->file, this->line)) << line << std::endl;
      } else {
          LOG_LOCATION(SourceLocation(this->file, this->line)) <<
          this->getResult()->getName() << " up-to-date with " <<
          "(" << left->getName() << ", " << right->getName() << ")" <<
          std::endl;
      }
    }

    size_t getLatestSourceVersion() override {
      return std::max(
        left->getLatestSourceVersion(),
        right->getLatestSourceVersion()
      );
    }

    operator std::string () const override {
      std::stringstream stream;
      stream << "contraction( " << this->alpha << ", " <<
        std::string(*left) << ", " << std::string(*right) << ", " <<
        this->beta << " )";
      return stream.str();
    }

  protected:
    /**
     * \brief Creates a contraction operation contracting the results of
     * the left and the right sub-operations where the result is to be
     * stored in the specified result tensor.
     **/
    static Ptr<ContractionOperation<F,TE>> create(
      const Ptr<IndexedTensorOperation<F,TE>> &left_,
      const Ptr<IndexedTensorOperation<F,TE>> &right_,
      const Ptr<Tensor<F,TE>> &result_,
      const char *resultIndices_,
      const Costs &contractionCosts,
      const Scope &scope
    ) {
      return New<ContractionOperation<F,TE>>(
        left_, right_,
        result_, resultIndices_,
        contractionCosts,
        scope.file, scope.line,
        typename Operation<TE>::ProtectedToken()
      );
    }

    Costs contractionCosts;
    Ptr<IndexedTensorOperation<F,TE>> left;
    Ptr<IndexedTensorOperation<F,TE>> right;

    friend class Contraction<F,TE>;
  };
}

#endif

