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

#ifndef __CXXGRAPH_UNDIRECTEDWEIGHTEDEDGE_H__
#define __CXXGRAPH_UNDIRECTEDWEIGHTEDEDGE_H__

#pragma once

#include "UndirectedEdge.hpp"
#include "Weighted.hpp"

namespace CXXGraph {
// Smart pointers alias
template <typename T>
using unique = std::unique_ptr<T>;
template <typename T>
using shared= std::shared_ptr<T>;

using std::make_unique;
using std::make_shared;

// Foward Declaration
template <typename T, typename Data>
class DirectedWeightedEdge;

template <typename T, typename Data>
class UndirectedWeightedEdge;

// ostream operator
template <typename T, typename Data>
std::ostream &operator<<(std::ostream &o,
                         const UndirectedWeightedEdge<T, Data> &edge);

template <typename T, typename Data>
class UndirectedWeightedEdge : public UndirectedEdge<T, Data>, public Weighted {
 public:
  UndirectedWeightedEdge(const unsigned long id, const Node<T> &node1,
                         const Node<T> &node2, const double weight);
  UndirectedWeightedEdge(const unsigned long id, shared<const Node<T>> node1,
                         shared<const Node<T>> node2, const double weight);
  UndirectedWeightedEdge(
      const unsigned long id,
      const std::pair<const Node<T> *, const Node<T> *> &nodepair,
      const double weight);
  UndirectedWeightedEdge(
      const unsigned long id,
      const std::pair<shared<const Node<T>>, shared<const Node<T>>> &nodepair,
      const double weight);
  UndirectedWeightedEdge(const UndirectedEdge<T, Data> &edge, const double weight);
  UndirectedWeightedEdge(const Edge<T, Data> &edge, const double weight);
  UndirectedWeightedEdge(const UndirectedEdge<T, Data> &edge);
  UndirectedWeightedEdge(const Edge<T, Data> &edge);
  UndirectedWeightedEdge(const DirectedWeightedEdge<T, Data> &edge);
  virtual ~UndirectedWeightedEdge() = default;
  const std::optional<bool> isWeighted() const override;
  // operator
  explicit operator DirectedWeightedEdge<T, Data>() const {
    return DirectedWeightedEdge<T, Data>(Edge<T, Data>::getId(), Edge<T, Data>::getNodePair(),
                                   Weighted::getWeight());
  }

  friend std::ostream &operator<< <>(std::ostream &os,
                                     const UndirectedWeightedEdge<T, Data> &edge);
};

template <typename T, typename Data>
UndirectedWeightedEdge<T, Data>::UndirectedWeightedEdge(const unsigned long id,
                                                  const Node<T> &node1,
                                                  const Node<T> &node2,
                                                  const double weight)
    : UndirectedEdge<T, Data>(id, node1, node2), Weighted(weight) {}

template <typename T, typename Data>
UndirectedWeightedEdge<T, Data>::UndirectedWeightedEdge(const unsigned long id,
                                                  shared<const Node<T>> node1,
                                                  shared<const Node<T>> node2,
                                                  const double weight)
    : UndirectedEdge<T, Data>(id, node1, node2), Weighted(weight) {}

template <typename T, typename Data>
UndirectedWeightedEdge<T, Data>::UndirectedWeightedEdge(
    const unsigned long id,
    const std::pair<const Node<T> *, const Node<T> *> &nodepair,
    const double weight)
    : UndirectedEdge<T, Data>(id, nodepair), Weighted(weight) {}

template <typename T, typename Data>
UndirectedWeightedEdge<T, Data>::UndirectedWeightedEdge(
    const unsigned long id,
    const std::pair<shared<const Node<T>>, shared<const Node<T>>> &nodepair,
    const double weight)
    : UndirectedEdge<T, Data>(id, nodepair), Weighted(weight) {}

template <typename T, typename Data>
UndirectedWeightedEdge<T, Data>::UndirectedWeightedEdge(const UndirectedEdge<T, Data> &edge,
                                                  const double weight)
    : UndirectedEdge<T, Data>(edge), Weighted(weight) {}

template <typename T, typename Data>
UndirectedWeightedEdge<T, Data>::UndirectedWeightedEdge(const Edge<T, Data> &edge,
                                                  const double weight)
    : UndirectedEdge<T, Data>(edge), Weighted(weight) {}

template <typename T, typename Data>
UndirectedWeightedEdge<T, Data>::UndirectedWeightedEdge(const UndirectedEdge<T, Data> &edge)
    : UndirectedEdge<T, Data>(edge), Weighted() {}

template <typename T, typename Data>
UndirectedWeightedEdge<T, Data>::UndirectedWeightedEdge(const Edge<T, Data> &edge)
    : UndirectedEdge<T, Data>(edge), Weighted() {}

template <typename T, typename Data>
UndirectedWeightedEdge<T, Data>::UndirectedWeightedEdge(
    const DirectedWeightedEdge<T, Data> &edge)
    : UndirectedEdge<T, Data>(edge), Weighted(edge.getWeight()) {}

template <typename T, typename Data>
const std::optional<bool> UndirectedWeightedEdge<T, Data>::isWeighted() const {
  return true;
}

template <typename T, typename Data>
std::ostream &operator<<(std::ostream &os,
                         const UndirectedWeightedEdge<T, Data> &edge) {
  os << "((Node: " << edge.getNode1().getId() << ")) <----- |Edge: #"
     << edge.getId() << " W:" << edge.getWeight()
     << "|-----> ((Node: " << edge.getNode2().getId() << "))";
  return os;
}

}  // namespace CXXGraph

#endif  // __CXXGRAPH_UNDIRECTEDWEIGHTEDEDGE_H__
