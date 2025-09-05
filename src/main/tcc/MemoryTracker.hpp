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

#ifndef TCC_MEMORY_TRACKER_DEFINED
#define TCC_MEMORY_TRACKER_DEFINED

#include <cstddef>
namespace cc4s {

class MemoryTracker {
public:
    virtual ~MemoryTracker() = default;

    // ---- Global API ----
    virtual void global_reset() = 0;
    virtual size_t current_bytes_allocated() const = 0;
    virtual size_t high_watermark_bytes() const = 0;
    virtual size_t number_of_allocations() const = 0;
    virtual size_t number_of_deallocations() const = 0;

    // ---- Local API (scoped measurements) ----
    virtual size_t reset() = 0; // local reset
    virtual size_t local_high_watermark() const = 0;
    virtual size_t local_allocations() const = 0;
    virtual size_t local_deallocations() const = 0;

    // ---- Controls ----
    virtual void set_size_threshold(size_t threshold) = 0;
};

class NullMemoryTracker : public MemoryTracker {
public:
    NullMemoryTracker() = default;
    ~NullMemoryTracker() override = default;
    void global_reset() override {}
    size_t current_bytes_allocated() const override { return 0; }
    size_t high_watermark_bytes() const override { return 0; }
    size_t number_of_allocations() const override { return 0; }
    size_t number_of_deallocations() const override { return 0; }
    size_t reset() override { return 0; }
    size_t local_high_watermark() const override { return 0; }
    size_t local_allocations() const override { return 0; }
    size_t local_deallocations() const override { return 0; }
    void set_size_threshold(size_t) override {}
};


// Forward-declare CTF-specific factory (if CTF is enabled)
#if defined(CTF_MEMORY_LOGGER)
class CtfMemoryTracker;
std::unique_ptr<MemoryTracker> createCtfMemoryTracker();
#endif

// Factory function (defined here, but CTF-specific part is conditional)
inline std::unique_ptr<MemoryTracker> createMemoryTracker() {
#if defined(CTF_MEMORY_LOGGER)
    return createCtfMemoryTracker();
#else
    return std::unique_ptr<MemoryTracker>(new NullMemoryTracker());
#endif
}


/* dont need this guy anymore
// Factory function forward-declared
std::unique_ptr<MemoryTracker> createMemoryTracker();
*/
}


#endif
