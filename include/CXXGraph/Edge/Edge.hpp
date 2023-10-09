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

#include "CXXGraph/Node/Node.hpp"
#include "CXXGraph/Node/DataNode.hpp"
#include "CXXGraph/Utility/id_t.hpp"

namespace CXXGraph {
// Smart pointers alias
template <typename T>
using unique = std::unique_ptr<T>;
template <typename T>
using shared = std::shared_ptr<T>;

using std::make_unique;
using std::make_shared;

class Edge;
// ostream operator
std::ostream &operator<<(std::ostream &o, const Edge &edge);
class Edge {
 private:
  CXXGraph::id_t id = 0;
  std::pair<shared<const Node>, shared<const Node>> nodePair;

 public:
  Edge(const CXXGraph::id_t id, const Node &node1, const Node &node2);
  Edge(const CXXGraph::id_t id, shared<const Node> node1, shared<const Node> node2);
  Edge(const CXXGraph::id_t id,
       const std::pair<const Node *, const Node *> &nodepair);
  Edge(const CXXGraph::id_t id,
       const std::pair<shared<const Node>, shared<const Node>> &nodepair);
  virtual ~Edge() = default;
  void setFirstNode(shared<const Node> node);
  void setSecondNode(shared<const Node> node);
  const unsigned long long getId() const;
  const std::pair<shared<const Node>, shared<const Node>> &getNodePair() const;
  shared<const Node> getOtherNode(shared<const Node> node) const;
  virtual const std::optional<bool> isDirected() const;
  virtual const std::optional<bool> isWeighted() const;
  // operator
  virtual bool operator==(const Edge &b) const;
  bool operator<(const Edge &b) const;
  // operator DirectedEdge() const { return DirectedEdge(id, nodePair); }
  // operator UndirectedEdge() const { return UndirectedEdge(id,
  // nodePair); }

  friend std::ostream &operator<<(std::ostream &os, const Edge &edge);
};

Edge::Edge(const CXXGraph::id_t id, const Node &node1,
              const Node &node2) {
  this->nodePair.first = make_shared<const Node>(node1);
  this->nodePair.second = make_shared<const Node>(node2);
  this->id = id;
}

Edge::Edge(const CXXGraph::id_t id, shared<const Node> node1, shared<const Node> node2) {
  this->nodePair.first = node1;
  this->nodePair.second = node2;
  this->id = id;
}

Edge::Edge(const CXXGraph::id_t id,
              const std::pair<const Node *, const Node *> &nodepair) {
  this->nodePair.first = make_shared<const Node>(*(nodepair.first));
  this->nodePair.second = make_shared<const Node>(*(nodepair.second));
  this->id = id;
}

Edge::Edge(const CXXGraph::id_t id,
              const std::pair<shared<const Node>, shared<const Node>> &nodepair)
    : nodePair(nodepair) {
  this->id = id;
}

void Edge::setFirstNode(shared<const Node> node) {
  /* this->nodePair = std::make_pair(node, this->nodePair.second); */
  this->nodePair.first = node;
}

void Edge::setSecondNode(shared<const Node> node) {
  /* this->nodePair = std::make_pair(this->nodePair.first, node); */
  this->nodePair.second = node;
}

const unsigned long long Edge::getId() const {
  return id;
}

const std::pair<shared<const Node>, shared<const Node>> &Edge::getNodePair()
    const {
  return nodePair;
}

shared<const Node> Edge::getOtherNode(shared<const Node> node) const {
  if (this->getNodePair().first == node) {
    return this->getNodePair().second;
  } else {
    return this->getNodePair().first;
  }
}

const std::optional<bool> Edge::isDirected() const {
  return std::nullopt;
}

const std::optional<bool> Edge::isWeighted() const {
  return std::nullopt;
}

bool Edge::operator==(const Edge &b) const {
  return (this->id == b.id && this->nodePair == b.nodePair);
}

bool Edge::operator<(const Edge &b) const {
  return (this->id < b.id);
}

std::ostream &operator<<(std::ostream &os, const Edge &edge) {
  os << "((Node: " << edge.nodePair.first->getId()
     << ")) ?----- |Edge: " << edge.id
     << "|-----? ((Node: " << edge.nodePair.second->getId() << "))";
  return os;
}
}  // namespace CXXGraph

#endif  // __CXXGRAPH_EDGE_H__
