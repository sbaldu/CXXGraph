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
using shared = std::shared_ptr<T>;

using std::make_shared;
using std::make_unique;

// Foward Declaration
class DirectedWeightedEdge;

class UndirectedWeightedEdge;

// ostream operator
std::ostream &operator<<(std::ostream &o, const UndirectedWeightedEdge &edge);

class UndirectedWeightedEdge : public UndirectedEdge, public Weighted {
 public:
  UndirectedWeightedEdge(const CXXGraph::id_t id, const Node &node1,
                         const Node &node2, const double weight);
  UndirectedWeightedEdge(const CXXGraph::id_t id, shared<const Node> node1,
                         shared<const Node> node2, const double weight);
  UndirectedWeightedEdge(const CXXGraph::id_t id,
                         const std::pair<const Node *, const Node *> &nodepair,
                         const double weight);
  UndirectedWeightedEdge(
      const CXXGraph::id_t id,
      const std::pair<shared<const Node>, shared<const Node>> &nodepair,
      const double weight);
  UndirectedWeightedEdge(const UndirectedEdge &edge, const double weight);
  UndirectedWeightedEdge(const Edge &edge, const double weight);
  UndirectedWeightedEdge(const UndirectedEdge &edge);
  UndirectedWeightedEdge(const Edge &edge);
  UndirectedWeightedEdge(const DirectedWeightedEdge &edge);
  virtual ~UndirectedWeightedEdge() = default;
  const std::optional<bool> isWeighted() const override;
  // operator
  explicit operator DirectedWeightedEdge() const {
    return DirectedWeightedEdge(Edge::getId(), Edge<T>::getNodePair(),
                                Weighted::getWeight());
  }

  friend std::ostream &operator<<(std::ostream &os,
                                  const UndirectedWeightedEdge &edge);
};

UndirectedWeightedEdge::UndirectedWeightedEdge(const CXXGraph::id_t id,
                                               const Node &node1,
                                               const Node &node2,
                                               const double weight)
    : UndirectedEdge(id, node1, node2), Weighted(weight) {}

UndirectedWeightedEdge::UndirectedWeightedEdge(const CXXGraph::id_t id,
                                               shared<const Node> node1,
                                               shared<const Node> node2,
                                               const double weight)
    : UndirectedEdge(id, node1, node2), Weighted(weight) {}

UndirectedWeightedEdge::UndirectedWeightedEdge(
    const CXXGraph::id_t id,
    const std::pair<const Node *, const Node *> &nodepair, const double weight)
    : UndirectedEdge(id, nodepair), Weighted(weight) {}

UndirectedWeightedEdge::UndirectedWeightedEdge(
    const CXXGraph::id_t id,
    const std::pair<shared<const Node>, shared<const Node>> &nodepair,
    const double weight)
    : UndirectedEdge(id, nodepair), Weighted(weight) {}

UndirectedWeightedEdge::UndirectedWeightedEdge(const UndirectedEdge &edge,
                                               const double weight)
    : UndirectedEdge(edge), Weighted(weight) {}

UndirectedWeightedEdge::UndirectedWeightedEdge(const Edge &edge,
                                               const double weight)
    : UndirectedEdge(edge), Weighted(weight) {}

UndirectedWeightedEdge::UndirectedWeightedEdge(const UndirectedEdge &edge)
    : UndirectedEdge(edge), Weighted() {}

UndirectedWeightedEdge::UndirectedWeightedEdge(const Edge &edge)
    : UndirectedEdge(edge), Weighted() {}

UndirectedWeightedEdge::UndirectedWeightedEdge(const DirectedWeightedEdge &edge)
    : UndirectedEdge(edge), Weighted(edge.getWeight()) {}

const std::optional<bool> UndirectedWeightedEdge::isWeighted() const {
  return true;
}

std::ostream &operator<<(std::ostream &os, const UndirectedWeightedEdge &edge) {
  os << "((Node: " << edge.getNode1().getId() << ")) <----- |Edge: #"
     << edge.getId() << " W:" << edge.getWeight()
     << "|-----> ((Node: " << edge.getNode2().getId() << "))";
  return os;
}

}  // namespace CXXGraph

#endif  // __CXXGRAPH_UNDIRECTEDWEIGHTEDEDGE_H__
