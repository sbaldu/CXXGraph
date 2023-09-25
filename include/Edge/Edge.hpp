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

#ifndef __CXXGRAPH_EDGE_H__
#define __CXXGRAPH_EDGE_H__

#pragma once

#include <memory>
#include <optional>
#include <utility>

#include "Node/Node.hpp"

namespace CXXGraph {
// Smart pointers alias
template <typename T>
using unique = std::unique_ptr<T>;
template <typename T>
using shared = std::shared_ptr<T>;

using std::make_shared;
using std::make_unique;

template <typename T, typename Data>
class Edge;
// ostream operator
template <typename T, typename Data>
std::ostream &operator<<(std::ostream &o, const Edge<T, Data> &edge);

template <typename T, typename Data>
class Edge {
 private:
  unsigned long long id = 0;
  Data data;
  std::pair<shared<const Node<T>>, shared<const Node<T>>> nodePair;

 public:
  Edge(const unsigned long long id, const Node<T> &node1, const Node<T> &node2);
  Edge(const unsigned long long id, shared<const Node<T>> node1,
       shared<const Node<T>> node2);
  Edge(const unsigned long long id,
       const std::pair<const Node<T> *, const Node<T> *> &nodepair);
  Edge(const unsigned long long id,
       const std::pair<shared<const Node<T>>, shared<const Node<T>>> &nodepair);
  Edge(const unsigned long long id, const Node<T> &node1, const Node<T> &node2,
       Data data);
  Edge(const unsigned long long id, shared<const Node<T>> node1,
       shared<const Node<T>> node2, Data data);
  Edge(const unsigned long long id,
       const std::pair<const Node<T> *, const Node<T> *> &nodepair, Data data);
  Edge(const unsigned long long id,
       const std::pair<shared<const Node<T>>, shared<const Node<T>>> &nodepair,
       Data data);
  virtual ~Edge() = default;
  void setFirstNode(shared<const Node<T>> node);
  void setSecondNode(shared<const Node<T>> node);
  const unsigned long long &getId() const;
  const std::pair<shared<const Node<T>>, shared<const Node<T>>> &getNodePair()
      const;
  shared<const Node<T>> getOtherNode(shared<const Node<T>> node) const;
  virtual const std::optional<bool> isDirected() const;
  virtual const std::optional<bool> isWeighted() const;
  // operator
  virtual bool operator==(const Edge<T, Data> &b) const;
  bool operator<(const Edge<T, Data> &b) const;
  // operator DirectedEdge<T>() const { return DirectedEdge<T>(id, nodePair); }
  // operator UndirectedEdge<T>() const { return UndirectedEdge<T>(id,
  // nodePair); }

  friend std::ostream &operator<< <>(std::ostream &os,
                                     const Edge<T, Data> &edge);
};

template <typename T, typename Data>
Edge<T, Data>::Edge(const unsigned long long id, const Node<T> &node1,
                    const Node<T> &node2) {
  this->nodePair.first = make_shared<const Node<T>>(node1);
  this->nodePair.second = make_shared<const Node<T>>(node2);
  this->id = id;
}

template <typename T, typename Data>
Edge<T, Data>::Edge(const unsigned long long id, shared<const Node<T>> node1,
                    shared<const Node<T>> node2) {
  this->nodePair.first = node1;
  this->nodePair.second = node2;
  this->id = id;
}

template <typename T, typename Data>
Edge<T, Data>::Edge(
    const unsigned long long id,
    const std::pair<const Node<T> *, const Node<T> *> &nodepair) {
  this->nodePair.first = make_shared<const Node<T>>(*(nodepair.first));
  this->nodePair.second = make_shared<const Node<T>>(*(nodepair.second));
  this->id = id;
}

template <typename T, typename Data>
Edge<T, Data>::Edge(
    const unsigned long long id,
    const std::pair<shared<const Node<T>>, shared<const Node<T>>> &nodepair)
    : nodePair(nodepair) {
  this->id = id;
}

template <typename T, typename Data>
Edge<T, Data>::Edge(const unsigned long long id, const Node<T> &node1,
                    const Node<T> &node2, Data data)
    : data{std::move(data)} {
  this->nodePair.first = make_shared<const Node<T>>(node1);
  this->nodePair.second = make_shared<const Node<T>>(node2);
  this->id = id;
}

template <typename T, typename Data>
Edge<T, Data>::Edge(const unsigned long long id, shared<const Node<T>> node1,
                    shared<const Node<T>> node2, Data data)
    : data{std::move(data)} {
  this->nodePair.first = node1;
  this->nodePair.second = node2;
  this->id = id;
}

template <typename T, typename Data>
Edge<T, Data>::Edge(const unsigned long long id,
                    const std::pair<const Node<T> *, const Node<T> *> &nodepair,
                    Data data)
    : data{std::move(data)} {
  this->nodePair.first = make_shared<const Node<T>>(*(nodepair.first));
  this->nodePair.second = make_shared<const Node<T>>(*(nodepair.second));
  this->id = id;
}

template <typename T, typename Data>
Edge<T, Data>::Edge(
    const unsigned long long id,
    const std::pair<shared<const Node<T>>, shared<const Node<T>>> &nodepair,
    Data data)
    : nodePair(nodepair), data{std::move(data)} {
  this->id = id;
}

template <typename T, typename Data>
void Edge<T, Data>::setFirstNode(shared<const Node<T>> node) {
  /* this->nodePair = std::make_pair(node, this->nodePair.second); */
  this->nodePair.first = node;
}

template <typename T, typename Data>
void Edge<T, Data>::setSecondNode(shared<const Node<T>> node) {
  /* this->nodePair = std::make_pair(this->nodePair.first, node); */
  this->nodePair.second = node;
}

template <typename T, typename Data>
const unsigned long long &Edge<T, Data>::getId() const {
  return id;
}

template <typename T, typename Data>
const std::pair<shared<const Node<T>>, shared<const Node<T>>>
    &Edge<T, Data>::getNodePair() const {
  return nodePair;
}

template <typename T, typename Data>
shared<const Node<T>> Edge<T, Data>::getOtherNode(
    shared<const Node<T>> node) const {
  if (this->getNodePair().first == node) {
    return this->getNodePair().second;
  } else {
    return this->getNodePair().first;
  }
}

template <typename T, typename Data>
const std::optional<bool> Edge<T, Data>::isDirected() const {
  return std::nullopt;
}

template <typename T, typename Data>
const std::optional<bool> Edge<T, Data>::isWeighted() const {
  return std::nullopt;
}

template <typename T, typename Data>
bool Edge<T, Data>::operator==(const Edge<T, Data> &b) const {
  return (this->id == b.id && this->nodePair == b.nodePair);
}

template <typename T, typename Data>
bool Edge<T, Data>::operator<(const Edge<T, Data> &b) const {
  return (this->id < b.id);
}

template <typename T, typename Data>
std::ostream &operator<<(std::ostream &os, const Edge<T, Data> &edge) {
  os << "((Node: " << edge.nodePair.first->getId()
     << ")) ?----- |Edge: " << edge.id
     << "|-----? ((Node: " << edge.nodePair.second->getId() << "))";
  return os;
}
}  // namespace CXXGraph

#endif  // __CXXGRAPH_EDGE_H__
