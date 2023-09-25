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

#ifndef __CXXGRAPH_DIRECTEDEDGE_H__
#define __CXXGRAPH_DIRECTEDEDGE_H__

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

template <typename T, typename Data>
class DirectedEdge;
// ostream operator
template <typename T, typename Data>
std::ostream &operator<<(std::ostream &o, const DirectedEdge<T, Data> &edge);
template <typename T, typename Data>
class DirectedEdge : public Edge<T, Data> {
 public:
  DirectedEdge(const unsigned long id, const Node<T> &node1,
               const Node<T> &node2);
  DirectedEdge(const unsigned long id, shared<const Node<T>> node1,
               shared<const Node<T>> node2);
  DirectedEdge(const unsigned long id,
               const std::pair<const Node<T> *, const Node<T> *> &nodepair);
  DirectedEdge(const unsigned long id,
               const std::pair<shared<const Node<T>>, shared<const Node<T>>> &nodepair);
  DirectedEdge(const unsigned long id, const Node<T> &node1,
               const Node<T> &node2, Data data);
  DirectedEdge(const unsigned long id, shared<const Node<T>> node1,
               shared<const Node<T>> node2, Data data);
  DirectedEdge(const unsigned long id,
               const std::pair<const Node<T> *, const Node<T> *> &nodepair, Data data);
  DirectedEdge(const unsigned long id,
               const std::pair<shared<const Node<T>>, shared<const Node<T>>> &nodepair, Data data);
  DirectedEdge(const Edge<T, Data> &edge);
  virtual ~DirectedEdge() = default;
  const Node<T> &getFrom() const;
  const Node<T> &getTo() const;
  const std::optional<bool> isDirected() const override;
  const std::optional<bool> isWeighted() const override;
  // operator
  explicit operator UndirectedEdge<T, Data>() const {
    return UndirectedEdge<T, Data>(Edge<T, Data>::getId(), Edge<T, Data>::getNodePair());
  }

  friend std::ostream &operator<< <>(std::ostream &os,
                                     const DirectedEdge<T, Data> &edge);
};

template <typename T, typename Data>
DirectedEdge<T, Data>::DirectedEdge(const unsigned long id, const Node<T> &node1,
                              const Node<T> &node2)
    : Edge<T, Data>(id, node1, node2) {}

template <typename T, typename Data>
DirectedEdge<T, Data>::DirectedEdge(const unsigned long id, shared<const Node<T>> node1,
			 shared<const Node<T>> node2) : Edge<T, Data>(id, node1, node2) {}

template <typename T, typename Data>
DirectedEdge<T, Data>::DirectedEdge(
    const unsigned long id,
    const std::pair<const Node<T> *, const Node<T> *> &nodepair)
    : Edge<T, Data>(id, nodepair) {}

template <typename T, typename Data>
DirectedEdge<T, Data>::DirectedEdge(
    const unsigned long id,
    const std::pair<shared<const Node<T>>, shared<const Node<T>>> &nodepair)
    : Edge<T, Data>(id, nodepair) {}

template <typename T, typename Data>
DirectedEdge<T, Data>::DirectedEdge(const Edge<T, Data> &edge)
    : DirectedEdge(edge.getId(), *(edge.getNodePair().first),
                   *(edge.getNodePair().second)) {}

template <typename T, typename Data>
DirectedEdge<T, Data>::DirectedEdge(const unsigned long id, const Node<T> &node1,
                              const Node<T> &node2, Data data)
    : Edge<T, Data>(id, node1, node2, data) {}

template <typename T, typename Data>
DirectedEdge<T, Data>::DirectedEdge(const unsigned long id, shared<const Node<T>> node1,
			 shared<const Node<T>> node2, Data data) : Edge<T, Data>(id, node1, node2, data) {}

template <typename T, typename Data>
DirectedEdge<T, Data>::DirectedEdge(
    const unsigned long id,
    const std::pair<const Node<T> *, const Node<T> *> &nodepair, Data data)
    : Edge<T, Data>(id, nodepair, data) {}

template <typename T, typename Data>
DirectedEdge<T, Data>::DirectedEdge(
    const unsigned long id,
    const std::pair<shared<const Node<T>>, shared<const Node<T>>> &nodepair, Data data)
    : Edge<T, Data>(id, nodepair, data) {}

template <typename T, typename Data>
const Node<T> &DirectedEdge<T, Data>::getFrom() const {
  return *(Edge<T, Data>::getNodePair().first);
}

template <typename T, typename Data>
const Node<T> &DirectedEdge<T, Data>::getTo() const {
  return *(Edge<T, Data>::getNodePair().second);
}

template <typename T, typename Data>
const std::optional<bool> DirectedEdge<T, Data>::isDirected() const {
  return true;
}

template <typename T, typename Data>
const std::optional<bool> DirectedEdge<T, Data>::isWeighted() const {
  return false;
}

template <typename T, typename Data>
std::ostream &operator<<(std::ostream &os, const DirectedEdge<T, Data> &edge) {
  os << "((Node: " << edge.getFrom().getId() << ")) +----- |Edge: #"
     << edge.getId() << "|-----> ((Node: " << edge.getTo().getId() << "))";
  return os;
}
}  // namespace CXXGraph

#endif  // __CXXGRAPH_DIRECTEDEDGE_H__
