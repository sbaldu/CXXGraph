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
using shared = std::shared_ptr<T>;

using std::make_shared;
using std::make_unique;

class UndirectedEdge;

// ostream operator
std::ostream &operator<<(std::ostream &o, const UndirectedEdge &edge);

class UndirectedEdge : public Edge {
 public:
  UndirectedEdge(const CXXGraph::id_t id, const Node &node1, const Node &node2);
  UndirectedEdge(const CXXGraph::id_t id, shared<const Node> node1,
                 shared<const Node> node2);
  UndirectedEdge(const CXXGraph::id_t id,
                 const std::pair<const Node *, const Node *> &nodepair);
  UndirectedEdge(
      const CXXGraph::id_t id,
      const std::pair<shared<const Node>, shared<const Node>> &nodepair);
  UndirectedEdge(const Edge &edge);
  virtual ~UndirectedEdge() = default;
  const Node &getNode1() const;
  const Node &getNode2() const;
  const std::optional<bool> isDirected() const override;
  const std::optional<bool> isWeighted() const override;
  // operator
  explicit operator DirectedEdge() const {
    return DirectedEdge(Edge::getId(), Edge::getNodePair());
  }

  friend std::ostream &operator<<(std::ostream &os, const UndirectedEdge &edge);
};

UndirectedEdge::UndirectedEdge(const CXXGraph::id_t id, const Node &node1,
                               const Node &node2)
    : Edge(id, node1, node2) {}

UndirectedEdge::UndirectedEdge(const CXXGraph::id_t id,
                               shared<const Node> node1,
                               shared<const Node> node2)
    : Edge(id, node1, node2) {}

UndirectedEdge::UndirectedEdge(
    const CXXGraph::id_t id,
    const std::pair<const Node *, const Node *> &nodepair)
    : Edge(id, nodepair) {}

UndirectedEdge::UndirectedEdge(
    const CXXGraph::id_t id,
    const std::pair<shared<const Node>, shared<const Node>> &nodepair)
    : Edge(id, nodepair) {}

UndirectedEdge::UndirectedEdge(const Edge &edge)
    : UndirectedEdge(edge.getId(), *(edge.getNodePair().first),
                     *(edge.getNodePair().second)) {}

const Node &UndirectedEdge::getNode1() const {
  return *(Edge::getNodePair().first);
}

const Node &UndirectedEdge::getNode2() const {
  return *(Edge::getNodePair().second);
}

const std::optional<bool> UndirectedEdge::isDirected() const { return false; }

const std::optional<bool> UndirectedEdge::isWeighted() const { return false; }

std::ostream &operator<<(std::ostream &os, const UndirectedEdge &edge) {
  os << "((Node: " << edge.getNode1().getId() << ")) <----- |Edge: #"
     << edge.getId() << "|-----> ((Node: " << edge.getNode2().getId() << "))";
  return os;
}
}  // namespace CXXGraph

#endif  // __CXXGRAPH_UNDIRECTEDEDGE_H__
