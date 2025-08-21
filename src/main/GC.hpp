#ifndef _CC4S_GC_HPP_
#define _CC4S_GC_HPP_
#include <Node.hpp>

namespace cc4s {

void collect_garbage(size_t const i, Ptr<MapNode> &steps,
                     Ptr<MapNode> &storage);

} // namespace cc4s

#endif
