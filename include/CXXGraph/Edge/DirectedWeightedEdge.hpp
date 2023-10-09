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
using shared = std::shared_ptr<T>;

using std::make_shared;
using std::make_unique;

// Foward Declaration
class UndirectedWeightedEdge;

class DirectedWeightedEdge;

// ostream operator
std::ostream &operator<<(std::ostream &o, const DirectedWeightedEdge &edge);

class DirectedWeightedEdge : public DirectedEdge, public Weighted {
 public:
  DirectedWeightedEdge(const CXXGraph::id_t id, const Node &node1,
                       const Node &node2, const double weight);
  DirectedWeightedEdge(const CXXGraph::id_t id, shared<const Node> node1,
                       shared<const Node> node2, const double weight);
  DirectedWeightedEdge(const CXXGraph::id_t id,
                       const std::pair<const Node *, const Node *> &nodepair,
                       const double weight);
  DirectedWeightedEdge(
      const CXXGraph::id_t id,
      const std::pair<shared<const Node>, shared<const Node>> &nodepair,
      const double weight);
  DirectedWeightedEdge(const DirectedEdge &edge, const double weight);
  DirectedWeightedEdge(const Edge &edge, const double weight);
  DirectedWeightedEdge(const DirectedEdge &edge);
  DirectedWeightedEdge(const Edge &edge);
  DirectedWeightedEdge(const UndirectedWeightedEdge &edge);
  virtual ~DirectedWeightedEdge() = default;
  const std::optional<bool> isWeighted() const override;
  // operator
  explicit operator UndirectedWeightedEdge() const {
    return UndirectedWeightedEdge(Edge::getId(), Edge::getNodePair(),
                                  Weighted::getWeight());
  }

  friend std::ostream &operator<<(std::ostream &os,
                                  const DirectedWeightedEdge &edge);
};

DirectedWeightedEdge::DirectedWeightedEdge(const CXXGraph::id_t id,
                                           const Node &node1, const Node &node2,
                                           const double weight)
    : DirectedEdge(id, node1, node2), Weighted(weight) {}

DirectedWeightedEdge::DirectedWeightedEdge(const CXXGraph::id_t id,
                                           shared<const Node> node1,
                                           shared<const Node> node2,
                                           const double weight)
    : DirectedEdge(id, node1, node2), Weighted(weight) {}

DirectedWeightedEdge::DirectedWeightedEdge(
    const CXXGraph::id_t id,
    const std::pair<const Node *, const Node *> &nodepair, const double weight)
    : DirectedEdge(id, nodepair), Weighted(weight) {}

DirectedWeightedEdge::DirectedWeightedEdge(
    const CXXGraph::id_t id,
    const std::pair<shared<const Node>, shared<const Node>> &nodepair,
    const double weight)
    : DirectedEdge(id, nodepair), Weighted(weight) {}

DirectedWeightedEdge::DirectedWeightedEdge(const DirectedEdge &edge,
                                           const double weight)
    : DirectedEdge(edge), Weighted(weight) {}

DirectedWeightedEdge::DirectedWeightedEdge(const Edge &edge,
                                           const double weight)
    : DirectedEdge(edge), Weighted(weight) {}

DirectedWeightedEdge::DirectedWeightedEdge(const DirectedEdge &edge)
    : DirectedEdge(edge), Weighted() {}

DirectedWeightedEdge::DirectedWeightedEdge(const Edge &edge)
    : DirectedEdge(edge), Weighted() {}

DirectedWeightedEdge::DirectedWeightedEdge(const UndirectedWeightedEdge &edge)
    : DirectedEdge(edge), Weighted(edge.getWeight()) {}

const std::optional<bool> DirectedWeightedEdge::isWeighted() const {
  return true;
}

std::ostream &operator<<(std::ostream &os, const DirectedWeightedEdge &edge) {
  os << "((Node: " << edge.getFrom().getId() << ")) +----- |Edge: #"
     << edge.getId() << " W:" << edge.getWeight()
     << "|-----> ((Node: " << edge.getTo().getId() << "))";
  return os;
}

}  // namespace CXXGraph

#endif  // __CXXGRAPH_DIRECTEDWEIGHTEDEDGE_H__
