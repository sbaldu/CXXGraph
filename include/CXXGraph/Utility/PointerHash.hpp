/***********************************************************/
/***      ______  ____  ______                 _         ***/
/***     / ___\ \/ /\ \/ / ___|_ __ __ _ _ __ | |__	     ***/
/***    | |    \  /  \  / |  _| '__/ _` | '_ \| '_ \	 ***/
/***    | |___ /  \  /  \ |_| | | | (_| | |_) | | | |    ***/
/***     \____/_/\_\/_/\_\____|_|  \__,_| .__/|_| |_|    ***/
/***                                    |_|			     ***/
/***********************************************************/
/***     Header-Only C++ Library for Graph			     ***/
/***	 Representation and Algorithms				     ***/
/***********************************************************/
/***     Author: ZigRazor ***/
/***	 E-Mail: zigrazor@gmail.com 				     ***/
/***********************************************************/
/***	 Collaboration: ----------- 				     ***/
/***********************************************************/
/***	 License: AGPL v3.0 ***/
/***********************************************************/

#ifndef __CXXGRAPH_POINTER_HASH__
#define __CXXGRAPH_POINTER_HASH__

#pragma once

#include <memory>
#include <string>

#include "CXXGraph/Edge/DirectedEdge.hpp"
#include "CXXGraph/Edge/DirectedWeightedEdge.hpp"
#include "CXXGraph/Edge/Edge.hpp"
#include "CXXGraph/Edge/UndirectedEdge.hpp"
#include "CXXGraph/Edge/UndirectedWeightedEdge.hpp"
#include "CXXGraph/Edge/Weighted.hpp"
#include "CXXGraph/Node/Node.hpp"

namespace CXXGraph {
template <typename T>
using shared = std::shared_ptr<T>;

// Redefine the hash functions and equality operators for shared pointers of nodes and edges
struct nodeHash {
  size_t operator()(const shared<const Node>& node) const {
    return node->getId();
  }
  size_t operator()(const shared<Node>& node) const {
    return node->getId();
  }
};

struct edgeHash {
  size_t operator()(const shared<const Edge>& edge) const {
    return (edge->getNodePair().first->getId()) ^ (edge->getNodePair().second->getId());
  }
};

bool operator==(shared<const Node> p1, shared<const Node> p2) {
  return p1->getUserId() == p2->getUserId();
}

bool operator==(shared<Node> p1, shared<Node> p2) {
  return p1->getUserId() == p2->getUserId();
}

bool operator==(shared<const Edge> p1, shared<const Edge> p2) {
  return p1->getNodePair().first->getUserId() ==
             p2->getNodePair().first->getUserId() &&
         p1->getNodePair().second->getUserId() ==
             p2->getNodePair().second->getUserId();
}
}  // namespace CXXGraph

#endif
