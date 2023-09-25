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
#ifndef __CXXGRAPH_DIRECTEDWEIGHTEDEDGE_H__
#define __CXXGRAPH_DIRECTEDWEIGHTEDEDGE_H__

#pragma once

#include "DirectedEdge.hpp"
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
class UndirectedWeightedEdge;

template <typename T, typename Data>
class DirectedWeightedEdge;

// ostream operator
template <typename T, typename Data>
std::ostream &operator<<(std::ostream &o, const DirectedWeightedEdge<T, Data> &edge);

template <typename T, typename Data>
class DirectedWeightedEdge : public DirectedEdge<T, Data>, public Weighted {
 public:
  DirectedWeightedEdge(const unsigned long id, const Node<T> &node1,
                       const Node<T> &node2, const double weight);
  DirectedWeightedEdge(const unsigned long id, shared<const Node<T>> node1,
                       shared<const Node<T>> node2, const double weight);
  DirectedWeightedEdge(
      const unsigned long id,
      const std::pair<const Node<T> *, const Node<T> *> &nodepair,
      const double weight);
  DirectedWeightedEdge(
      const unsigned long id,
      const std::pair<shared<const Node<T>>, shared<const Node<T>>> &nodepair,
      const double weight);
  DirectedWeightedEdge(const DirectedEdge<T, Data> &edge, const double weight);
  DirectedWeightedEdge(const Edge<T, Data> &edge, const double weight);
  DirectedWeightedEdge(const DirectedEdge<T, Data> &edge);
  DirectedWeightedEdge(const Edge<T, Data> &edge);
  DirectedWeightedEdge(const UndirectedWeightedEdge<T, Data> &edge);
  virtual ~DirectedWeightedEdge() = default;
  const std::optional<bool> isWeighted() const override;
  // operator
  explicit operator UndirectedWeightedEdge<T, Data>() const {
    return UndirectedWeightedEdge<T, Data>(Edge<T, Data>::getId(), Edge<T>::getNodePair(),
                                     Weighted::getWeight());
  }

  friend std::ostream &operator<< <>(std::ostream &os,
                                     const DirectedWeightedEdge<T, Data> &edge);
};

template <typename T, typename Data>
DirectedWeightedEdge<T, Data>::DirectedWeightedEdge(const unsigned long id,
                                              const Node<T> &node1,
                                              const Node<T> &node2,
                                              const double weight)
    : DirectedEdge<T, Data>(id, node1, node2), Weighted(weight) {}

template <typename T, typename Data>
DirectedWeightedEdge<T, Data>::DirectedWeightedEdge(const unsigned long id,
                                              shared<const Node<T>> node1,
                                              shared<const Node<T>> node2,
                                              const double weight)
    : DirectedEdge<T, Data>(id, node1, node2), Weighted(weight) {}

template <typename T, typename Data>
DirectedWeightedEdge<T, Data>::DirectedWeightedEdge(
    const unsigned long id,
    const std::pair<const Node<T> *, const Node<T> *> &nodepair,
    const double weight)
    : DirectedEdge<T, Data>(id, nodepair), Weighted(weight) {}

template <typename T, typename Data>
DirectedWeightedEdge<T, Data>::DirectedWeightedEdge(
    const unsigned long id,
    const std::pair<shared<const Node<T>>, shared<const Node<T>>> &nodepair,
    const double weight)
    : DirectedEdge<T, Data>(id, nodepair), Weighted(weight) {}

template <typename T, typename Data>
DirectedWeightedEdge<T, Data>::DirectedWeightedEdge(const DirectedEdge<T, Data> &edge,
                                              const double weight)
    : DirectedEdge<T, Data>(edge), Weighted(weight) {}

template <typename T, typename Data>
DirectedWeightedEdge<T, Data>::DirectedWeightedEdge(const Edge<T, Data> &edge,
                                              const double weight)
    : DirectedEdge<T, Data>(edge), Weighted(weight) {}

template <typename T, typename Data>
DirectedWeightedEdge<T, Data>::DirectedWeightedEdge(const DirectedEdge<T, Data> &edge)
    : DirectedEdge<T, Data>(edge), Weighted() {}

template <typename T, typename Data>
DirectedWeightedEdge<T, Data>::DirectedWeightedEdge(const Edge<T, Data> &edge)
    : DirectedEdge<T, Data>(edge), Weighted() {}

template <typename T, typename Data>
DirectedWeightedEdge<T, Data>::DirectedWeightedEdge(
    const UndirectedWeightedEdge<T, Data> &edge)
    : DirectedEdge<T, Data>(edge), Weighted(edge.getWeight()) {}

template <typename T, typename Data>
const std::optional<bool> DirectedWeightedEdge<T, Data>::isWeighted() const {
  return true;
}

template <typename T, typename Data>
std::ostream &operator<<(std::ostream &os,
                         const DirectedWeightedEdge<T, Data> &edge) {
  os << "((Node: " << edge.getFrom().getId() << ")) +----- |Edge: #"
     << edge.getId() << " W:" << edge.getWeight()
     << "|-----> ((Node: " << edge.getTo().getId() << "))";
  return os;
}

}  // namespace CXXGraph

#endif  // __CXXGRAPH_DIRECTEDWEIGHTEDEDGE_H__
