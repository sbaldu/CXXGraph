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

#ifndef __CXXGRAPH_BRONKERBASCH_IMPL_H__
#define __CXXGRAPH_BRONKERBASCH_IMPL_H__

#pragma once

#include <unordered_set>

#include "CXXGraph/Graph/Graph_decl.h"

namespace CXXGraph {
template <typename T>
virtual const BronKerboschResult<T> Graph<T>::BronKerbosch(
    std::unordered_set<shared<const Node<T>>>& P,
    std::unordered_set<shared<const Node<T>>>& R,
    std::unordered_set<shared<const Node<T>>>& X,
    BronKerboschResult<T>& result) const {
  if (P.size() == 0 && X.size() == 0) {
    result.push_back(R);
  } else {
	for (const auto& nodeIt : P) {
	  if (nodeIt)
	}
  }
}

template <typename T>
virtual const BronKerboschResult<T> Graph<T>::BronKerbosch() const {
  std::unordered_set<shared<const Node<T>>> P{this->getNodeSet()};
  std::unordered_set<shared<const Node<T>>> R{};
  std::unordered_set<shared<const Node<T>>> X{};

  // choose a pivot
  BronKerboschResult<T> result;
  BronKerbosch(P, R, X, result);
  return result;
}
}  // namespace CXXGraph

#endif  // __CXXGRAPH_BELLMANFORD_IMPL_H__
