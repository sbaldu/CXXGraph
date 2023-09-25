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

#ifndef __CXXGRAPH_UNDIRECTEDEDGE_H__
#define __CXXGRAPH_UNDIRECTEDEDGE_H__

#pragma once

#include "Edge.hpp"

namespace CXXGraph {
// Smart pointers alias
template <typename T>
using unique = std::unique_ptr<T>;
template <typename T>
using shared= std::shared_ptr<T>;

using std::make_unique;
using std::make_shared;

template <typename T, typename Data>
class UndirectedEdge;

// ostream operator
template <typename T, typename Data>
std::ostream &operator<<(std::ostream &o, const UndirectedEdge<T, Data> &edge);

template <typename T, typename Data>
class UndirectedEdge : public Edge<T, Data> {
 public:
  UndirectedEdge(const unsigned long id, const Node<T> &node1,
                 const Node<T> &node2);
  UndirectedEdge(const unsigned long id, shared<const Node<T>> node1,
                 shared<const Node<T>> node2);
  UndirectedEdge(const unsigned long id,
                 const std::pair<const Node<T> *, const Node<T> *> &nodepair);
  UndirectedEdge(const unsigned long id,
                 const std::pair<shared<const Node<T>>, shared<const Node<T>>> &nodepair);
  UndirectedEdge(const Edge<T, Data> &edge);
  virtual ~UndirectedEdge() = default;
  const Node<T> &getNode1() const;
  const Node<T> &getNode2() const;
  const std::optional<bool> isDirected() const override;
  const std::optional<bool> isWeighted() const override;
  // operator
  explicit operator DirectedEdge<T, Data>() const {
    return DirectedEdge<T, Data>(Edge<T, Data>::getId(), Edge<T, Data>::getNodePair());
  }

  friend std::ostream &operator<< <>(std::ostream &os,
                                     const UndirectedEdge<T, Data> &edge);
};

template <typename T, typename Data>
UndirectedEdge<T, Data>::UndirectedEdge(const unsigned long id, const Node<T> &node1,
                                  const Node<T> &node2)
    : Edge<T, Data>(id, node1, node2) {}

template <typename T, typename Data>
UndirectedEdge<T, Data>::UndirectedEdge(const unsigned long id, shared<const Node<T>> node1,
                                  shared<const Node<T>> node2)
    : Edge<T, Data>(id, node1, node2) {}

template <typename T, typename Data>
UndirectedEdge<T, Data>::UndirectedEdge(
    const unsigned long id,
    const std::pair<const Node<T> *, const Node<T> *> &nodepair)
    : Edge<T, Data>(id, nodepair) {}

template <typename T, typename Data>
UndirectedEdge<T, Data>::UndirectedEdge(
    const unsigned long id,
    const std::pair<shared<const Node<T>>, shared<const Node<T>>> &nodepair)
    : Edge<T, Data>(id, nodepair) {}

template <typename T, typename Data>
UndirectedEdge<T, Data>::UndirectedEdge(const Edge<T, Data> &edge)
    : UndirectedEdge(edge.getId(), *(edge.getNodePair().first),
                     *(edge.getNodePair().second)) {}

template <typename T, typename Data>
const Node<T> &UndirectedEdge<T, Data>::getNode1() const {
  return *(Edge<T, Data>::getNodePair().first);
}

template <typename T, typename Data>
const Node<T> &UndirectedEdge<T, Data>::getNode2() const {
  return *(Edge<T, Data>::getNodePair().second);
}

template <typename T, typename Data>
const std::optional<bool> UndirectedEdge<T, Data>::isDirected() const {
  return false;
}

template <typename T, typename Data>
const std::optional<bool> UndirectedEdge<T, Data>::isWeighted() const {
  return false;
}

template <typename T, typename Data>
std::ostream &operator<<(std::ostream &os, const UndirectedEdge<T, Data> &edge) {
  os << "((Node: " << edge.getNode1().getId() << ")) <----- |Edge: #"
     << edge.getId() << "|-----> ((Node: " << edge.getNode2().getId() << "))";
  return os;
}
}  // namespace CXXGraph

#endif  // __CXXGRAPH_UNDIRECTEDEDGE_H__
