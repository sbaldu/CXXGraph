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
using shared = std::shared_ptr<T>;

using std::make_shared;
using std::make_unique;

class UndirectedEdge;

class DirectedEdge;
// ostream operator
std::ostream &operator<<(std::ostream &o, const DirectedEdge &edge);
class DirectedEdge : public Edge {
 public:
  DirectedEdge(const CXXGraph::id_t id, const Node &node1, const Node &node2);
  DirectedEdge(const CXXGraph::id_t id, shared<const Node> node1,
               shared<const Node> node2);
  DirectedEdge(const CXXGraph::id_t id,
               const std::pair<const Node *, const Node *> &nodepair);
  DirectedEdge(
      const CXXGraph::id_t id,
      const std::pair<shared<const Node>, shared<const Node>> &nodepair);
  DirectedEdge(const Edge &edge);
  virtual ~DirectedEdge() = default;
  const Node &getFrom() const;
  const Node &getTo() const;
  const std::optional<bool> isDirected() const override;
  const std::optional<bool> isWeighted() const override;
  // operator
  explicit operator UndirectedEdge() const {
    return UndirectedEdge(Edge::getId(), Edge::getNodePair());
  }

  friend std::ostream &operator<<(std::ostream &os, const DirectedEdge &edge);
};

DirectedEdge::DirectedEdge(const CXXGraph::id_t id, const Node &node1,
                           const Node &node2)
    : Edge(id, node1, node2) {}

DirectedEdge::DirectedEdge(const CXXGraph::id_t id, shared<const Node> node1,
                           shared<const Node> node2)
    : Edge(id, node1, node2) {}

DirectedEdge::DirectedEdge(
    const CXXGraph::id_t id,
    const std::pair<const Node *, const Node *> &nodepair)
    : Edge(id, nodepair) {}

DirectedEdge::DirectedEdge(
    const CXXGraph::id_t id,
    const std::pair<shared<const Node>, shared<const Node>> &nodepair)
    : Edge(id, nodepair) {}

DirectedEdge::DirectedEdge(const Edge &edge)
    : DirectedEdge(edge.getId(), *(edge.getNodePair().first),
                   *(edge.getNodePair().second)) {}

const Node &DirectedEdge::getFrom() const {
  return *(Edge::getNodePair().first);
}

const Node &DirectedEdge::getTo() const {
  return *(Edge::getNodePair().second);
}

const std::optional<bool> DirectedEdge::isDirected() const { return true; }

const std::optional<bool> DirectedEdge::isWeighted() const { return false; }

std::ostream &operator<<(std::ostream &os, const DirectedEdge &edge) {
  os << "((Node: " << edge.getFrom().getId() << ")) +----- |Edge: #"
     << edge.getId() << "|-----> ((Node: " << edge.getTo().getId() << "))";
  return os;
}
}  // namespace CXXGraph

#endif  // __CXXGRAPH_DIRECTEDEDGE_H__
