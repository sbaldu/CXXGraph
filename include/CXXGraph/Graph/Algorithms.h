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

#ifndef __CXXGRAPH_GRAPH_H__
#define __CXXGRAPH_GRAPH_H__

#include <cstdio>
#pragma once

#include <limits.h>

#include <atomic>
#include <cmath>
#include <condition_variable>
#include <cstring>
#include <deque>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "CXXGraph/Edge/DirectedEdge.hpp"
#include "CXXGraph/Edge/DirectedWeightedEdge.hpp"
#include "CXXGraph/Edge/Edge.hpp"
#include "CXXGraph/Edge/UndirectedEdge.hpp"
#include "CXXGraph/Edge/UndirectedWeightedEdge.hpp"
#include "CXXGraph/Edge/Weighted.hpp"
#include "CXXGraph/Graph/Graph.hpp"
#include "CXXGraph/Node/Node.hpp"
#include "CXXGraph/Partitioning/Partition.hpp"
#include "CXXGraph/Partitioning/PartitionAlgorithm.hpp"
#include "CXXGraph/Partitioning/Partitioner.hpp"
#include "CXXGraph/Partitioning/Utility/Globals.hpp"
#include "CXXGraph/Utility/ConstString.hpp"
#include "CXXGraph/Utility/ConstValue.hpp"
#include "CXXGraph/Utility/PointerHash.hpp"
#include "CXXGraph/Utility/Reader.hpp"
#include "CXXGraph/Utility/ThreadSafe.hpp"
#include "CXXGraph/Utility/Typedef.hpp"
#include "CXXGraph/Utility/Writer.hpp"

#ifdef WITH_COMPRESSION
#include <zlib.h>
#endif

namespace CXXGraph {
// Smart pointers alias
template <typename T>
using unique = std::unique_ptr<T>;
template <typename T>
using shared = std::shared_ptr<T>;

using std::make_shared;
using std::make_unique;

template <typename T>
using T_EdgeSet = std::unordered_set<shared<const Edge<T>>, edgeHash<T>>;

template <typename T>
using T_NodeSet = std::unordered_set<shared<const Node<T>>, nodeHash<T>>;

namespace Partitioning {
template <typename T>
class Partition;
}

template <typename T>
std::ostream &operator<<(std::ostream &o, const Graph<T> &graph);
template <typename T>
std::ostream &operator<<(std::ostream &o, const AdjacencyMatrix<T> &adj);

template <typename T>
const DijkstraResult Graph<T>::dijkstra(const Node<T> &source,
                                        const Node<T> &target) const {
  DijkstraResult result;
  auto nodeSet = Graph<T>::getNodeSet();

  auto source_node_it = std::find_if(
      nodeSet.begin(), nodeSet.end(),
      [&source](auto node) { return node->getUserId() == source.getUserId(); });
  if (source_node_it == nodeSet.end()) {
    // check if source node exist in the graph
    result.errorMessage = ERR_SOURCE_NODE_NOT_IN_GRAPH;
    return result;
  }

  auto target_node_it = std::find_if(
      nodeSet.begin(), nodeSet.end(),
      [&target](auto node) { return node->getUserId() == target.getUserId(); });
  if (target_node_it == nodeSet.end()) {
    // check if target node exist in the graph
    result.errorMessage = ERR_TARGET_NODE_NOT_IN_GRAPH;
    return result;
  }
  const auto adj = getAdjMatrix();
  // n denotes the number of vertices in graph
  auto n = adj->size();

  // setting all the distances initially to INF_DOUBLE
  std::unordered_map<shared<const Node<T>>, double, nodeHash<T>> dist;

  for (const auto &node : nodeSet) {
    dist[node] = INF_DOUBLE;
  }

  // creating a min heap using priority queue
  // first element of pair contains the distance
  // second element of pair contains the vertex
  std::priority_queue<std::pair<double, shared<const Node<T>>>,
                      std::vector<std::pair<double, shared<const Node<T>>>>,
                      std::greater<std::pair<double, shared<const Node<T>>>>>
      pq;

  // pushing the source vertex 's' with 0 distance in min heap
  pq.push(std::make_pair(0.0, *source_node_it));

  // marking the distance of source as 0
  dist[*source_node_it] = 0;

  std::unordered_map<std::string, std::string> parent;
  parent[source.getUserId()] = "";

  while (!pq.empty()) {
    // second element of pair denotes the node / vertex
    shared<const Node<T>> currentNode = pq.top().second;
    // first element of pair denotes the distance
    double currentDist = pq.top().first;

    pq.pop();

    // for all the reachable vertex from the currently exploring vertex
    // we will try to minimize the distance
    if (adj->find(currentNode) != adj->end()) {
      for (const auto &elem : adj->at(currentNode)) {
        // minimizing distances
        if (elem.second->isWeighted().has_value() &&
            elem.second->isWeighted().value()) {
          if (elem.second->isDirected().has_value() &&
              elem.second->isDirected().value()) {
            shared<const DirectedWeightedEdge<T>> dw_edge =
                std::static_pointer_cast<const DirectedWeightedEdge<T>>(
                    elem.second);
            if (dw_edge->getWeight() < 0) {
              result.errorMessage = ERR_NEGATIVE_WEIGHTED_EDGE;
              return result;
            } else if (currentDist + dw_edge->getWeight() < dist[elem.first]) {
              dist[elem.first] = currentDist + dw_edge->getWeight();
              pq.push(std::make_pair(dist[elem.first], elem.first));
              parent[elem.first.get()->getUserId()] =
                  currentNode.get()->getUserId();
            }
          } else if (elem.second->isDirected().has_value() &&
                     !elem.second->isDirected().value()) {
            shared<const UndirectedWeightedEdge<T>> udw_edge =
                std::static_pointer_cast<const UndirectedWeightedEdge<T>>(
                    elem.second);
            if (udw_edge->getWeight() < 0) {
              result.errorMessage = ERR_NEGATIVE_WEIGHTED_EDGE;
              return result;
            } else if (currentDist + udw_edge->getWeight() < dist[elem.first]) {
              dist[elem.first] = currentDist + udw_edge->getWeight();
              pq.push(std::make_pair(dist[elem.first], elem.first));
              parent[elem.first.get()->getUserId()] =
                  currentNode.get()->getUserId();
            }
          } else {
            // ERROR it shouldn't never returned ( does not exist a Node
            // Weighted and not Directed/Undirected)
            result.errorMessage = ERR_NO_DIR_OR_UNDIR_EDGE;
            return result;
          }
        } else {
          // No Weighted Edge
          result.errorMessage = ERR_NO_WEIGHTED_EDGE;
          return result;
        }
      }
    }
  }
  if (dist[*target_node_it] != INF_DOUBLE) {
    result.success = true;
    result.errorMessage = "";
    result.result = dist[*target_node_it];
    std::string id = target.getUserId();
    while (parent[id] != "") {
      result.path.push_back(id);
      id = parent[id];
    }
    result.path.push_back(source.getUserId());
    std::reverse(result.path.begin(), result.path.end());
    return result;
  }
  result.errorMessage = ERR_TARGET_NODE_NOT_REACHABLE;
  return result;
}

template <typename T>
const BellmanFordResult Graph<T>::bellmanford(const Node<T> &source,
                                              const Node<T> &target) const {
  BellmanFordResult result;
  result.success = false;
  result.errorMessage = "";
  result.result = INF_DOUBLE;
  auto nodeSet = Graph<T>::getNodeSet();
  auto source_node_it = std::find_if(
      nodeSet.begin(), nodeSet.end(),
      [&source](auto node) { return node->getUserId() == source.getUserId(); });
  if (source_node_it == nodeSet.end()) {
    // check if source node exist in the graph
    result.errorMessage = ERR_SOURCE_NODE_NOT_IN_GRAPH;
    return result;
  }
  auto target_node_it = std::find_if(
      nodeSet.begin(), nodeSet.end(),
      [&target](auto node) { return node->getUserId() == target.getUserId(); });
  if (target_node_it == nodeSet.end()) {
    // check if target node exist in the graph
    result.errorMessage = ERR_TARGET_NODE_NOT_IN_GRAPH;
    return result;
  }
  // setting all the distances initially to INF_DOUBLE
  std::unordered_map<shared<const Node<T>>, double, nodeHash<T>> dist,
      currentDist;
  // n denotes the number of vertices in graph
  auto n = nodeSet.size();
  for (const auto &elem : nodeSet) {
    dist[elem] = INF_DOUBLE;
    currentDist[elem] = INF_DOUBLE;
  }

  // marking the distance of source as 0
  dist[*source_node_it] = 0;
  // set if node distances in two consecutive
  // iterations remain the same.
  auto earlyStopping = false;
  // outer loop for vertex relaxation
  for (int i = 0; i < n - 1; ++i) {
    auto edgeSet = Graph<T>::getEdgeSet();
    // inner loop for distance updates of
    // each relaxation
    for (const auto &edge : edgeSet) {
      auto elem = edge->getNodePair();
      if (edge->isWeighted().has_value() && edge->isWeighted().value()) {
        auto edge_weight =
            (std::dynamic_pointer_cast<const Weighted>(edge))->getWeight();
        if (dist[elem.first] + edge_weight < dist[elem.second])
          dist[elem.second] = dist[elem.first] + edge_weight;
      } else {
        // No Weighted Edge
        result.errorMessage = ERR_NO_WEIGHTED_EDGE;
        return result;
      }
    }
    auto flag = true;
    for (const auto &[key, value] : dist) {
      if (currentDist[key] != value) {
        flag = false;
        break;
      }
    }
    for (const auto &[key, value] : dist) {
      currentDist[key] = value;  // update the current distance
    }
    if (flag) {
      earlyStopping = true;
      break;
    }
  }

  // check if there exists a negative cycle
  if (!earlyStopping) {
    auto edgeSet = Graph<T>::getEdgeSet();
    for (const auto &edge : edgeSet) {
      auto elem = edge->getNodePair();
      auto edge_weight =
          (std::dynamic_pointer_cast<const Weighted>(edge))->getWeight();
      if (dist[elem.first] + edge_weight < dist[elem.second]) {
        result.success = true;
        result.negativeCycle = true;
        result.errorMessage = "";
        return result;
      }
    }
  }

  if (dist[*target_node_it] != INF_DOUBLE) {
    result.success = true;
    result.errorMessage = "";
    result.negativeCycle = false;
    result.result = dist[*target_node_it];
    return result;
  }
  result.errorMessage = ERR_TARGET_NODE_NOT_REACHABLE;
  result.result = -1;
  return result;
}

/*
 * See Harry Hsu. "An algorithm for finding a minimal equivalent graph of a
 * digraph.", Journal of the ACM, 22(1):11-16, January 1975
 *
 * foreach x in graph.vertices
 *   foreach y in graph.vertices
 *     foreach z in graph.vertices
 *       delete edge xz if edges xy and yz exist
 */
template <typename T>
const Graph<T> Graph<T>::transitiveReduction() const {
  Graph<T> result(this->edgeSet);

  CXXGraph::id_t edgeId = 0;
  std::unordered_set<shared<const Node<T>>, nodeHash<T>> nodes =
      this->getNodeSet();
  for (auto x : nodes) {
    for (auto y : nodes) {
      if (this->findEdge(x, y, edgeId)) {
        for (auto z : nodes) {
          if (this->findEdge(y, z, edgeId)) {
            if (this->findEdge(x, z, edgeId)) {
              result.removeEdge(edgeId);
            }
          }
        }
      }
    }
  }

  return result;
}

template <typename T>
const FWResult Graph<T>::floydWarshall() const {
  FWResult result;
  result.success = false;
  result.errorMessage = "";
  std::unordered_map<std::pair<std::string, std::string>, double,
                     CXXGraph::pair_hash>
      pairwise_dist;
  const auto &nodeSet = Graph<T>::getNodeSet();
  // create a pairwise distance matrix with distance node distances
  // set to inf. Distance of node to itself is set as 0.
  for (const auto &elem1 : nodeSet) {
    for (const auto &elem2 : nodeSet) {
      auto key = std::make_pair(elem1->getUserId(), elem2->getUserId());
      if (elem1 != elem2)
        pairwise_dist[key] = INF_DOUBLE;
      else
        pairwise_dist[key] = 0.0;
    }
  }

  const auto &edgeSet = Graph<T>::getEdgeSet();
  // update the weights of nodes
  // connected by edges
  for (const auto &edge : edgeSet) {
    const auto &elem = edge->getNodePair();
    if (edge->isWeighted().has_value() && edge->isWeighted().value()) {
      auto edgeWeight =
          (std::dynamic_pointer_cast<const Weighted>(edge))->getWeight();
      auto key =
          std::make_pair(elem.first->getUserId(), elem.second->getUserId());
      pairwise_dist[key] = edgeWeight;
    } else {
      // if an edge exists but has no weight associated
      // with it, we return an error message
      result.errorMessage = ERR_NO_WEIGHTED_EDGE;
      return result;
    }
  }

  for (const auto &k : nodeSet) {
    // set all vertices as source one by one
    for (const auto &src : nodeSet) {
      // iterate through all vertices as destination for the
      // current source
      auto src_k = std::make_pair(src->getUserId(), k->getUserId());
      for (const auto &dst : nodeSet) {
        // If vertex k provides a shorter path than
        // src to dst, update the value of
        // pairwise_dist[src_to_dst]
        auto src_dst = std::make_pair(src->getUserId(), dst->getUserId());
        auto k_dst = std::make_pair(k->getUserId(), dst->getUserId());
        if (pairwise_dist[src_dst] >
                (pairwise_dist[src_k] + pairwise_dist[k_dst]) &&
            (pairwise_dist[k_dst] != INF_DOUBLE &&
             pairwise_dist[src_k] != INF_DOUBLE))
          pairwise_dist[src_dst] = pairwise_dist[src_k] + pairwise_dist[k_dst];
      }
    }
  }

  result.success = true;
  // presense of negative number in the diagonal indicates
  // that that the graph contains a negative cycle
  for (const auto &node : nodeSet) {
    auto diag = std::make_pair(node->getUserId(), node->getUserId());
    if (pairwise_dist[diag] < 0.) {
      result.negativeCycle = true;
      return result;
    }
  }
  result.result = std::move(pairwise_dist);
  return result;
}

template <typename T>
const MstResult Graph<T>::prim() const {
  MstResult result;
  result.success = false;
  result.errorMessage = "";
  result.mstCost = INF_DOUBLE;
  if (!isUndirectedGraph()) {
    result.errorMessage = ERR_DIR_GRAPH;
    return result;
  }
  if (!isConnectedGraph()) {
    result.errorMessage = ERR_NOT_STRONG_CONNECTED;
    return result;
  }
  auto nodeSet = Graph<T>::getNodeSet();
  auto n = nodeSet.size();
  const auto adj = Graph<T>::getAdjMatrix();

  // setting all the distances initially to INF_DOUBLE
  std::unordered_map<shared<const Node<T>>, double, nodeHash<T>> dist;
  for (const auto &elem : (*adj)) {
    dist[elem.first] = INF_DOUBLE;
  }

  // creating a min heap using priority queue
  // first element of pair contains the distance
  // second element of pair contains the vertex
  std::priority_queue<std::pair<double, shared<const Node<T>>>,
                      std::vector<std::pair<double, shared<const Node<T>>>>,
                      std::greater<std::pair<double, shared<const Node<T>>>>>
      pq;

  // pushing the source vertex 's' with 0 distance in min heap
  auto source = *(nodeSet.begin());
  pq.push(std::make_pair(0.0, source));
  result.mstCost = 0;
  std::vector<CXXGraph::id_t> doneNode;
  // mark source node as done
  // otherwise we get (0, 0) also in mst
  doneNode.push_back(source->getId());
  // stores the parent and corresponding child node
  // of the edges that are part of MST
  std::unordered_map<CXXGraph::id_t, std::string> parentNode;
  while (!pq.empty()) {
    // second element of pair denotes the node / vertex
    shared<const Node<T>> currentNode = pq.top().second;
    auto nodeId = currentNode->getId();
    if (std::find(doneNode.begin(), doneNode.end(), nodeId) == doneNode.end()) {
      auto pair = std::make_pair(parentNode[nodeId], currentNode->getUserId());
      result.mst.push_back(pair);
      result.mstCost += pq.top().first;
      doneNode.push_back(nodeId);
    }

    pq.pop();
    // for all the reachable vertex from the currently exploring vertex
    // we will try to minimize the distance
    if (adj->find(currentNode) != adj->end()) {
      for (const auto &elem : adj->at(currentNode)) {
        // minimizing distances
        if (elem.second->isWeighted().has_value() &&
            elem.second->isWeighted().value()) {
          shared<const UndirectedWeightedEdge<T>> udw_edge =
              std::static_pointer_cast<const UndirectedWeightedEdge<T>>(
                  elem.second);
          if ((udw_edge->getWeight() < dist[elem.first]) &&
              (std::find(doneNode.begin(), doneNode.end(),
                         elem.first->getId()) == doneNode.end())) {
            dist[elem.first] = udw_edge->getWeight();
            parentNode[elem.first->getId()] = currentNode->getUserId();
            pq.push(std::make_pair(dist[elem.first], elem.first));
          }
        } else {
          // No Weighted Edge
          result.errorMessage = ERR_NO_WEIGHTED_EDGE;
          return result;
        }
      }
    }
  }
  result.success = true;
  return result;
}

template <typename T>
const MstResult Graph<T>::boruvka() const {
  MstResult result;
  result.success = false;
  result.errorMessage = "";
  result.mstCost = INF_DOUBLE;
  if (!isUndirectedGraph()) {
    result.errorMessage = ERR_DIR_GRAPH;
    return result;
  }
  const auto nodeSet = Graph<T>::getNodeSet();
  const auto n = nodeSet.size();

  // Use std map for storing n subsets.
  auto subsets = make_shared<std::unordered_map<CXXGraph::id_t, Subset>>();

  // Initially there are n different trees.
  // Finally there will be one tree that will be MST
  auto numTrees = n;

  // check if all edges are weighted and store the weights
  // in a map whose keys are the edge ids and values are the edge weights
  const auto edgeSet = Graph<T>::getEdgeSet();
  std::unordered_map<CXXGraph::id_t, double> edgeWeight;
  for (const auto &edge : edgeSet) {
    if (edge->isWeighted().has_value() && edge->isWeighted().value())
      edgeWeight[edge->getId()] =
          (std::dynamic_pointer_cast<const Weighted>(edge))->getWeight();
    else {
      // No Weighted Edge
      result.errorMessage = ERR_NO_WEIGHTED_EDGE;
      return result;
    }
  }

  for (const auto &node : nodeSet) {
    Subset set{node->getId(), 0};
    (*subsets)[node->getId()] = set;
  }

  result.mstCost = 0;  // we will store the cost here
  // exit when only 1 tree i.e. mst
  while (numTrees > 1) {
    // Everytime initialize cheapest map
    // It stores index of the cheapest edge of subset.
    std::unordered_map<CXXGraph::id_t, CXXGraph::id_t> cheapest;
    for (const auto &node : nodeSet) cheapest[node->getId()] = INT_MAX;

    // Traverse through all edges and update
    // cheapest of every component
    for (const auto &edge : edgeSet) {
      auto elem = edge->getNodePair();
      auto edgeId = edge->getId();
      // Find sets of two corners of current edge
      auto set1 = Graph<T>::setFind(subsets, elem.first->getId());
      auto set2 = Graph<T>::setFind(subsets, elem.second->getId());

      // If two corners of current edge belong to
      // same set, ignore current edge
      if (set1 == set2) continue;

      // Else check if current edge is closer to previous
      // cheapest edges of set1 and set2
      if (cheapest[set1] == INT_MAX ||
          edgeWeight[cheapest[set1]] > edgeWeight[edgeId])
        cheapest[set1] = edgeId;

      if (cheapest[set2] == INT_MAX ||
          edgeWeight[cheapest[set2]] > edgeWeight[edgeId])
        cheapest[set2] = edgeId;
    }

    // iterate over all the vertices and add picked
    // cheapest edges to MST
    for (const auto &[nodeId, edgeId] : cheapest) {
      // Check if cheapest for current set exists
      if (edgeId != INT_MAX) {
        auto cheapestNode = Graph<T>::getEdge(edgeId).value()->getNodePair();
        auto set1 = Graph<T>::setFind(subsets, cheapestNode.first->getId());
        auto set2 = Graph<T>::setFind(subsets, cheapestNode.second->getId());
        if (set1 == set2) continue;
        result.mstCost += edgeWeight[edgeId];
        auto newEdgeMST = std::make_pair(cheapestNode.first->getUserId(),
                                         cheapestNode.second->getUserId());
        result.mst.push_back(newEdgeMST);
        // take union of set1 and set2 and decrease number of trees
        Graph<T>::setUnion(subsets, set1, set2);
        numTrees--;
      }
    }
  }
  result.success = true;
  return result;
}

template <typename T>
const MstResult Graph<T>::kruskal() const {
  MstResult result;
  result.success = false;
  result.errorMessage = "";
  result.mstCost = INF_DOUBLE;
  if (!isUndirectedGraph()) {
    result.errorMessage = ERR_DIR_GRAPH;
    return result;
  }
  const auto nodeSet = Graph<T>::getNodeSet();
  auto n = nodeSet.size();

  // check if all edges are weighted and store the weights
  // in a map whose keys are the edge ids and values are the edge weights
  auto edgeSet = Graph<T>::getEdgeSet();
  std::priority_queue<std::pair<double, shared<const Edge<T>>>,
                      std::vector<std::pair<double, shared<const Edge<T>>>>,
                      std::greater<std::pair<double, shared<const Edge<T>>>>>
      sortedEdges;
  for (const auto &edge : edgeSet) {
    if (edge->isWeighted().has_value() && edge->isWeighted().value()) {
      auto weight =
          (std::dynamic_pointer_cast<const Weighted>(edge))->getWeight();
      sortedEdges.push(std::make_pair(weight, edge));
    } else {
      // No Weighted Edge
      result.errorMessage = ERR_NO_WEIGHTED_EDGE;
      return result;
    }
  }

  auto subset = make_shared<std::unordered_map<CXXGraph::id_t, Subset>>();

  for (const auto &node : nodeSet) {
    Subset set{node->getId(), 0};
    (*subset)[node->getId()] = set;
  }
  result.mstCost = 0;
  while ((!sortedEdges.empty()) && (result.mst.size() < n)) {
    auto [edgeWeight, cheapestEdge] = sortedEdges.top();
    sortedEdges.pop();
    auto &[first, second] = cheapestEdge->getNodePair();
    auto set1 = Graph<T>::setFind(subset, first->getId());
    auto set2 = Graph<T>::setFind(subset, second->getId());
    if (set1 != set2) {
      result.mst.push_back(
          std::make_pair(first->getUserId(), second->getUserId()));
      result.mstCost += edgeWeight;
    }
    Graph<T>::setUnion(subset, set1, set2);
  }
  result.success = true;
  return result;
}

template <typename T>
BestFirstSearchResult<T> Graph<T>::best_first_search(
    const Node<T> &source, const Node<T> &target) const {
  BestFirstSearchResult<T> result;
  auto &nodeSet = Graph<T>::getNodeSet();
  using pq_type = std::pair<double, shared<const Node<T>>>;

  auto source_node_it = std::find_if(
      nodeSet.begin(), nodeSet.end(),
      [&source](auto node) { return node->getUserId() == source.getUserId(); });
  if (source_node_it == nodeSet.end()) {
    result.errorMessage = ERR_SOURCE_NODE_NOT_IN_GRAPH;
    return result;
  }

  auto target_node_it = std::find_if(
      nodeSet.begin(), nodeSet.end(),
      [&target](auto node) { return node->getUserId() == target.getUserId(); });
  if (target_node_it == nodeSet.end()) {
    result.errorMessage = ERR_TARGET_NODE_NOT_IN_GRAPH;
    return result;
  }

  auto adj = Graph<T>::getAdjMatrix();
  std::priority_queue<pq_type, std::vector<pq_type>, std::greater<pq_type>> pq;

  std::vector<Node<T>> visited;
  visited.push_back(source);
  pq.push(std::make_pair(0.0, *source_node_it));

  while (!pq.empty()) {
    shared<const Node<T>> currentNode = pq.top().second;
    pq.pop();
    result.nodesInBestSearchOrder.push_back(*currentNode);

    if (*currentNode == target) {
      break;
    }
    if (adj->find(currentNode) != adj->end()) {
      for (const auto &elem : adj->at(currentNode)) {
        if (elem.second->isWeighted().has_value()) {
          if (elem.second->isDirected().has_value()) {
            shared<const DirectedWeightedEdge<T>> dw_edge =
                std::static_pointer_cast<const DirectedWeightedEdge<T>>(
                    elem.second);
            if (std::find(visited.begin(), visited.end(), *(elem.first)) ==
                visited.end()) {
              visited.push_back(*(elem.first));
              pq.push(std::make_pair(dw_edge->getWeight(), elem.first));
            }
          } else {
            shared<const UndirectedWeightedEdge<T>> dw_edge =
                std::static_pointer_cast<const UndirectedWeightedEdge<T>>(
                    elem.second);
            if (std::find(visited.begin(), visited.end(), *(elem.first)) ==
                visited.end()) {
              visited.push_back(*(elem.first));
              pq.push(std::make_pair(dw_edge->getWeight(), elem.first));
            }
          }
        } else {
          result.errorMessage = ERR_NO_WEIGHTED_EDGE;
          result.nodesInBestSearchOrder.clear();
          return result;
        }
      }
    }
  }

  result.success = true;
  return result;
}

template <typename T>
const std::vector<Node<T>> Graph<T>::breadth_first_search(
    const Node<T> &start) const {
  // vector to keep track of visited nodes
  std::vector<Node<T>> visited;
  auto &nodeSet = Graph<T>::getNodeSet();
  // check is exist node in the graph
  auto start_node_it = std::find_if(
      nodeSet.begin(), nodeSet.end(),
      [&start](auto node) { return node->getUserId() == start.getUserId(); });
  if (start_node_it == nodeSet.end()) {
    return visited;
  }
  const auto adj = Graph<T>::getAdjMatrix();
  // queue that stores vertices that need to be further explored
  std::queue<shared<const Node<T>>> tracker;

  // mark the starting node as visited
  visited.push_back(start);
  tracker.push(*start_node_it);
  while (!tracker.empty()) {
    shared<const Node<T>> node = tracker.front();
    tracker.pop();
    if (adj->find(node) != adj->end()) {
      for (const auto &elem : adj->at(node)) {
        // if the node is not visited then mark it as visited
        // and push it to the queue
        if (std::find(visited.begin(), visited.end(), *(elem.first)) ==
            visited.end()) {
          visited.push_back(*(elem.first));
          tracker.push(elem.first);
        }
      }
    }
  }

  return visited;
}

template <typename T>
const std::vector<Node<T>> Graph<T>::concurrency_breadth_first_search(
    const Node<T> &start, size_t num_threads) const {
  std::vector<Node<T>> bfs_result;
  // check is exist node in the graph
  auto &nodeSet = Graph<T>::getNodeSet();
  auto start_node_it = std::find_if(
      nodeSet.begin(), nodeSet.end(),
      [&start](auto node) { return node->getUserId() == start.getUserId(); });
  if (start_node_it == nodeSet.end()) {
    return bfs_result;
  }

  std::unordered_map<shared<const Node<T>>, size_t, nodeHash<T>> node_to_index;
  for (const auto &node : nodeSet) {
    node_to_index[node] = node_to_index.size();
  }
  std::vector<size_t> visited(nodeSet.size(), 0);

  // parameter limitations
  if (num_threads <= 0) {
    std::cout << "Error: number of threads should be greater than 0"
              << std::endl;
    num_threads = 2;
  }

  const auto &adj = Graph<T>::getAdjMatrix();
  // vector that stores vertices to be visit
  std::vector<shared<const Node<T>>> level_tracker, next_level_tracker;
  level_tracker.reserve(static_cast<int>(1.0 * nodeSet.size()));
  next_level_tracker.reserve(static_cast<int>(1.0 * nodeSet.size()));

  // mark the starting node as visited
  visited[node_to_index[*start_node_it]] = 1;
  level_tracker.push_back(*start_node_it);

  // a worker is assigned a small part of tasks for each time
  // assignments of tasks in current level and updates of tasks in next
  // level are inclusive
  std::mutex tracker_mutex;
  std::mutex next_tracker_mutex;
  std::atomic<int> assigned_tasks = 0;
  int num_tasks = 1;
  // unit of task assignment, which mean assign block_size tasks to a
  // worker each time
  int block_size = 1;
  int level = 1;

  auto extract_tasks = [&level_tracker, &tracker_mutex, &assigned_tasks,
                        &num_tasks, &block_size]() -> std::pair<int, int> {
    /*
    std::lock_guard<std::mutex> tracker_guard(tracker_mutex);
    int task_block_size = std::min(num_tasks - assigned_tasks,
    block_size); std::pair<int,int> task_block{assigned_tasks,
    assigned_tasks + task_block_size}; assigned_tasks += task_block_size;
    return task_block;
    */
    int start = assigned_tasks.fetch_add(block_size);
    int end = std::min(num_tasks, start + block_size);
    return {start, end};
  };

  auto submit_result =
      [&next_level_tracker, &next_tracker_mutex](
          std::vector<shared<const Node<T>>> &submission) -> void {
    std::lock_guard<std::mutex> tracker_guard(next_tracker_mutex);
    next_level_tracker.insert(std::end(next_level_tracker),
                              std::begin(submission), std::end(submission));
  };

  // worker thread sleep until it begin to search nodes of next level
  std::mutex next_level_mutex;
  std::condition_variable next_level_cond;
  std::atomic<int> waiting_workers = 0;

  auto bfs_worker = [&]() -> void {
    // algorithm is not done
    while (!level_tracker.empty()) {
      // search for nodes in a level is not done
      std::vector<shared<const Node<T>>> local_tracker;
      while (true) {
        auto [start_index, end_index] = extract_tasks();
        if (start_index >= end_index) {
          break;
        }

        for (int i = start_index; i < end_index; ++i) {
          if (adj->count(level_tracker[i])) {
            for (const auto &elem : adj->at(level_tracker[i])) {
              int index = (int)node_to_index[elem.first];
              if (visited[index] == 0) {
                visited[index] = 1;
                local_tracker.push_back(elem.first);
              }
            }
          }
        }
      }

      // submit local result to global result
      if (!local_tracker.empty()) {
        submit_result(local_tracker);
      }

      // last worker need to do preparation for the next iteration
      int cur_level = level;
      if (num_threads == 1 + waiting_workers.fetch_add(1)) {
        swap(level_tracker, next_level_tracker);
        next_level_tracker.clear();

        // adjust block_size according to number of nodes in next level
        block_size = 4;
        if (level_tracker.size() <= num_threads * 4) {
          block_size = std::max(
              1, static_cast<int>(std::ceil(
                     static_cast<double>(level_tracker.size()) / num_threads)));
        } else if (level_tracker.size() >= num_threads * 64) {
          block_size = 16;
        }

        num_tasks = (int)level_tracker.size();
        waiting_workers = 0;
        assigned_tasks = 0;
        level = level + 1;
        next_level_cond.notify_all();
      } else {
        // not to wait if last worker reachs last statement before notify
        // all or even further
        std::unique_lock<std::mutex> next_level_lock(next_level_mutex);
        next_level_cond.wait(next_level_lock, [&level, cur_level]() {
          return level != cur_level;
        });
      }
    }
  };

  std::vector<std::thread> workers;
  for (int i = 0; i < num_threads - 1; ++i) {
    workers.emplace_back(std::thread(bfs_worker));
  }
  bfs_worker();

  for (auto &worker : workers) {
    if (worker.joinable()) {
      worker.join();
    }
  }

  for (const auto &visited_node : nodeSet) {
    if (visited[node_to_index[visited_node]] != 0) {
      bfs_result.push_back(*visited_node);
    }
  }

  return bfs_result;
}
template <typename T>
const std::vector<Node<T>> Graph<T>::depth_first_search(
    const Node<T> &start) const {
  // vector to keep track of visited nodes
  std::vector<Node<T>> visited;
  auto nodeSet = Graph<T>::getNodeSet();
  // check is exist node in the graph
  auto start_node_it = std::find_if(
      nodeSet.begin(), nodeSet.end(),
      [&start](auto node) { return node->getUserId() == start.getUserId(); });
  if (start_node_it == nodeSet.end()) {
    return visited;
  }
  const auto adj = Graph<T>::getAdjMatrix();
  std::function<void(const std::shared_ptr<AdjacencyMatrix<T>>,
                     shared<const Node<T>>, std::vector<Node<T>> &)>
      explore;
  explore = [&explore](const std::shared_ptr<AdjacencyMatrix<T>> adj,
                       shared<const Node<T>> node,
                       std::vector<Node<T>> &visited) -> void {
    visited.push_back(*node);
    if (adj->find(node) != adj->end()) {
      for (const auto &x : adj->at(node)) {
        if (std::find(visited.begin(), visited.end(), *(x.first)) ==
            visited.end()) {
          explore(adj, x.first, visited);
        }
      }
    }
  };
  explore(adj, *start_node_it, visited);

  return visited;
}

template <typename T>
bool Graph<T>::isCyclicDirectedGraphDFS() const {
  if (!isDirectedGraph()) {
    return false;
  }
  enum nodeStates : uint8_t { not_visited, in_stack, visited };
  auto nodeSet = Graph<T>::getNodeSet();
  auto adjMatrix = Graph<T>::getAdjMatrix();

  /* State of the node.
   *
   * It is a vector of "nodeStates" which represents the state node is in.
   * It can take only 3 values: "not_visited", "in_stack", and "visited".
   *
   * Initially, all nodes are in "not_visited" state.
   */
  std::unordered_map<CXXGraph::id_t, nodeStates> state;
  for (const auto &node : nodeSet) {
    state[node->getId()] = not_visited;
  }
  int stateCounter = 0;

  // Start visiting each node.
  for (const auto &node : nodeSet) {
    // If a node is not visited, only then check for presence of cycle.
    // There is no need to check for presence of cycle for a visited
    // node as it has already been checked for presence of cycle.
    if (state[node->getId()] == not_visited) {
      // Check for cycle.
      std::function<bool(const std::shared_ptr<AdjacencyMatrix<T>>,
                         std::unordered_map<CXXGraph::id_t, nodeStates> &,
                         shared<const Node<T>>)>
          isCyclicDFSHelper;
      isCyclicDFSHelper =
          [this, &isCyclicDFSHelper](
              const std::shared_ptr<AdjacencyMatrix<T>> adjMatrix,
              std::unordered_map<CXXGraph::id_t, nodeStates> &states,
              shared<const Node<T>> node) {
            // Add node "in_stack" state.
            states[node->getId()] = in_stack;

            // If the node has children, then recursively visit all
            // children of the node.
            auto const it = adjMatrix->find(node);
            if (it != adjMatrix->end()) {
              for (const auto &child : it->second) {
                // If state of child node is "not_visited", evaluate that
                // child for presence of cycle.
                auto state_of_child = states.at((std::get<0>(child))->getId());
                if (state_of_child == not_visited) {
                  if (isCyclicDFSHelper(adjMatrix, states,
                                        std::get<0>(child))) {
                    return true;
                  }
                } else if (state_of_child == in_stack) {
                  // If child node was "in_stack", then that means that
                  // there is a cycle in the graph. Return true for
                  // presence of the cycle.
                  return true;
                }
              }
            }

            // Current node has been evaluated for the presence of cycle
            // and had no cycle. Mark current node as "visited".
            states[node->getId()] = visited;
            // Return that current node didn't result in any cycles.
            return false;
          };
      if (isCyclicDFSHelper(adjMatrix, state, node)) {
        return true;
      }
    }
  }

  // All nodes have been safely traversed, that means there is no cycle in
  // the graph. Return false.
  return false;
}

template <typename T>
bool Graph<T>::containsCycle(const T_EdgeSet<T> *edgeSet) const {
  auto edgeSet_ptr = make_shared<const T_EdgeSet<T>>(*edgeSet);
  auto subset = make_shared<std::unordered_map<CXXGraph::id_t, Subset>>();
  // initialize the subset parent and rank values
  for (const auto &edge : *edgeSet_ptr) {
    auto &[first, second] = edge->getNodePair();
    std::vector<CXXGraph::id_t> nodeId(2);
    nodeId.push_back(first->getId());
    nodeId.push_back(second->getId());
    for (const auto &id : nodeId) {
      auto nodeExists = [id](const auto &it) {
        return (id == (it.second).parent);
      };

      if (std::find_if((*subset).begin(), (*subset).end(), nodeExists) ==
          (*subset).end()) {
        Subset set;
        set.parent = id;
        set.rank = 0;
        (*subset)[id] = set;
      }
    }
  }
  return Graph<T>::containsCycle(edgeSet_ptr, subset);
}

template <typename T>
bool Graph<T>::containsCycle(shared<const T_EdgeSet<T>> edgeSet) const {
  auto subset = make_shared<std::unordered_map<CXXGraph::id_t, Subset>>();
  // initialize the subset parent and rank values
  for (const auto &edge : *edgeSet) {
    auto &[first, second] = edge->getNodePair();
    std::vector<CXXGraph::id_t> nodeId(2);
    nodeId.push_back(first->getId());
    nodeId.push_back(second->getId());
    for (const auto &id : nodeId) {
      auto nodeExists = [id](const auto &it) {
        return (id == (it.second).parent);
      };

      if (std::find_if((*subset).begin(), (*subset).end(), nodeExists) ==
          (*subset).end()) {
        Subset set;
        set.parent = id;
        set.rank = 0;
        (*subset)[id] = set;
      }
    }
  }
  return Graph<T>::containsCycle(edgeSet, subset);
}

template <typename T>
bool Graph<T>::containsCycle(
    shared<const T_EdgeSet<T>> edgeSet,
    shared<std::unordered_map<CXXGraph::id_t, Subset>> subset) const {
  for (const auto &edge : *edgeSet) {
    auto &[first, second] = edge->getNodePair();
    auto set1 = Graph<T>::setFind(subset, first->getId());
    auto set2 = Graph<T>::setFind(subset, second->getId());
    if (set1 == set2) return true;
    Graph<T>::setUnion(subset, set1, set2);
  }
  return false;
}

template <typename T>
bool Graph<T>::isCyclicDirectedGraphBFS() const {
  if (!isDirectedGraph()) {
    return false;
  }
  const auto adjMatrix = Graph<T>::getAdjMatrix();
  auto nodeSet = Graph<T>::getNodeSet();

  std::unordered_map<CXXGraph::id_t, unsigned int> indegree;
  for (const auto &node : nodeSet) {
    indegree[node->getId()] = 0;
  }
  // Calculate the indegree i.e. the number of incident edges to the node.
  for (auto const &list : (*adjMatrix)) {
    auto children = list.second;
    for (auto const &child : children) {
      indegree[std::get<0>(child)->getId()]++;
    }
  }

  std::queue<shared<const Node<T>>> can_be_solved;
  for (const auto &node : nodeSet) {
    // If a node doesn't have any input edges, then that node will
    // definately not result in a cycle and can be visited safely.
    if (!indegree[node->getId()]) {
      can_be_solved.emplace(node);
    }
  }

  // Vertices that need to be traversed.
  auto remain = nodeSet.size();
  // While there are safe nodes that we can visit.
  while (!can_be_solved.empty()) {
    auto solved = can_be_solved.front();
    // Visit the node.
    can_be_solved.pop();
    // Decrease number of nodes that need to be traversed.
    remain--;

    // Visit all the children of the visited node.
    auto it = adjMatrix->find(solved);
    if (it != adjMatrix->end()) {
      for (const auto &child : it->second) {
        // Check if we can visited the node safely.
        if (--indegree[std::get<0>(child)->getId()] == 0) {
          // if node can be visited safely, then add that node to
          // the visit queue.
          can_be_solved.emplace(std::get<0>(child));
        }
      }
    }
  }

  // If there are still nodes that we can't visit, then it means that
  // there is a cycle and return true, else return false.
  return !(remain == 0);
}

template <typename T>
bool Graph<T>::isDirectedGraph() const {
  auto edgeSet = Graph<T>::getEdgeSet();
  for (const auto &edge : edgeSet) {
    if (!(edge->isDirected().has_value() && edge->isDirected().value())) {
      // Found Undirected Edge
      return false;
    }
  }
  // No Undirected Edge
  return true;
}

template <typename T>
bool Graph<T>::isUndirectedGraph() const {
  auto edgeSet = Graph<T>::getEdgeSet();
  for (const auto &edge : edgeSet) {
    if ((edge->isDirected().has_value() && edge->isDirected().value())) {
      // Found Directed Edge
      return false;
    }
  }
  // No Directed Edge
  return true;
}

template <typename T>
void Graph<T>::reverseDirectedGraph() {
  if (!isDirectedGraph()) {
    throw std::runtime_error(ERR_UNDIR_GRAPH);
  }
  auto oldEdgeSet = Graph<T>::getEdgeSet();
  for (const auto &edge : oldEdgeSet) {
    auto &[first, second] = edge->getNodePair();
    auto id = edge->getId();
    this->removeEdge(id);
    auto newEdge = std::make_shared<DirectedEdge<T>>(id, second, first);
    this->addEdge(newEdge);
  }
}

template <typename T>
bool Graph<T>::isConnectedGraph() const {
  if (!isUndirectedGraph()) {
    return false;
  } else {
    auto nodeSet = getNodeSet();
    const auto adjMatrix = getAdjMatrix();
    // created visited map
    std::unordered_map<CXXGraph::id_t, bool> visited;
    for (const auto &node : nodeSet) {
      visited[node->getId()] = false;
    }
    std::function<void(shared<const Node<T>>)> dfs_helper =
        [this, &adjMatrix, &visited,
         &dfs_helper](shared<const Node<T>> source) {
          // mark the vertex visited
          visited[source->getId()] = true;

          // travel the neighbors
          for (int i = 0; i < (*adjMatrix)[source].size(); i++) {
            shared<const Node<T>> neighbor = (*adjMatrix)[source].at(i).first;
            if (visited[neighbor->getId()] == false) {
              // make recursive call from neighbor
              dfs_helper(neighbor);
            }
          }
        };
    // call dfs_helper for the first node
    dfs_helper(*(nodeSet.begin()));

    // check if all the nodes are visited
    for (const auto &node : nodeSet) {
      if (visited[node->getId()] == false) {
        return false;
      }
    }
    return true;
  }
}

template <typename T>
bool Graph<T>::isStronglyConnectedGraph() const {
  if (!isDirectedGraph()) {
    return false;
  } else {
    auto nodeSet = getNodeSet();
    const auto adjMatrix = getAdjMatrix();
    for (const auto &start_node : nodeSet) {
      // created visited map
      std::unordered_map<CXXGraph::id_t, bool> visited;
      for (const auto &node : nodeSet) {
        visited[node->getId()] = false;
      }
      std::function<void(shared<const Node<T>>)> dfs_helper =
          [this, &adjMatrix, &visited,
           &dfs_helper](shared<const Node<T>> source) {
            // mark the vertex visited
            visited[source->getId()] = true;

            // travel the neighbors
            for (int i = 0; i < (*adjMatrix)[source].size(); i++) {
              shared<const Node<T>> neighbor = (*adjMatrix)[source].at(i).first;
              if (visited[neighbor->getId()] == false) {
                // make recursive call from neighbor
                dfs_helper(neighbor);
              }
            }
          };
      // call dfs_helper for the first node
      dfs_helper(start_node);

      // check if all the nodes are visited
      for (const auto &node : nodeSet) {
        if (visited[node->getId()] == false) {
          return false;
        }
      }
    }
    return true;
  }
}

template <typename T>
const TarjanResult<T> Graph<T>::tarjan(const unsigned int typeMask) const {
  TarjanResult<T> result;
  result.success = false;
  bool isDirected = this->isDirectedGraph();
  if (isDirected) {
    // check whether targetMask is a subset of the mask for directed graph
    unsigned int directedMask = TARJAN_FIND_SCC;
    if ((typeMask | directedMask) != directedMask) {
      result.errorMessage = ERR_DIR_GRAPH;
      return result;
    }
  } else {
    // check whether targetMask is a subset of the mask for undirected graph
    unsigned int undirectedMask = (TARJAN_FIND_CUTV | TARJAN_FIND_BRIDGE |
                                   TARJAN_FIND_VBCC | TARJAN_FIND_EBCC);
    if ((typeMask | undirectedMask) != undirectedMask) {
      result.errorMessage = ERR_UNDIR_GRAPH;
      return result;
    }
  }

  const auto &adjMatrix = getAdjMatrix();
  const auto &nodeSet = getNodeSet();
  std::unordered_map<CXXGraph::id_t, int>
      discoveryTime;  // the timestamp when a node is visited
  std::unordered_map<CXXGraph::id_t, int>
      lowestDisc;  // the lowest discovery time of all
                   // reachable nodes from current node
  int timestamp = 0;
  CXXGraph::id_t rootId = 0;
  std::stack<Node<T>> sccNodeStack;
  std::stack<Node<T>> ebccNodeStack;
  std::stack<Node<T>> vbccNodeStack;
  std::unordered_set<CXXGraph::id_t> inStack;
  std::function<void(const shared<const Node<T>>, const shared<const Edge<T>>)>
      dfs_helper = [this, typeMask, isDirected, &dfs_helper, &adjMatrix,
                    &discoveryTime, &lowestDisc, &timestamp, &rootId,
                    &sccNodeStack, &ebccNodeStack, &vbccNodeStack, &inStack,
                    &result](const shared<const Node<T>> curNode,
                             const shared<const Edge<T>> prevEdge) {
        // record the visited time of current node
        discoveryTime[curNode->getId()] = timestamp;
        lowestDisc[curNode->getId()] = timestamp;
        timestamp++;
        if (typeMask & TARJAN_FIND_SCC) {
          sccNodeStack.emplace(*curNode);
          inStack.emplace(curNode->getId());
        }
        if (typeMask & TARJAN_FIND_EBCC) {
          ebccNodeStack.emplace(*curNode);
        }
        if (typeMask & TARJAN_FIND_VBCC) {
          vbccNodeStack.emplace(*curNode);
        }
        // travel the neighbors
        int numSon = 0;
        bool nodeIsAdded =
            false;  // whether a node has been marked as a cut vertice
        if (adjMatrix->find(curNode) != adjMatrix->end()) {
          for (const auto &[neighborNode, edge] : adjMatrix->at(curNode)) {
            if (!discoveryTime.count(neighborNode->getId())) {
              dfs_helper(neighborNode, edge);
              lowestDisc[curNode->getId()] =
                  std::min(lowestDisc[curNode->getId()],
                           lowestDisc[neighborNode->getId()]);

              if (typeMask & TARJAN_FIND_BRIDGE) {
                // lowestDisc of neighbor node is larger than that of current
                // node means we can travel back to a visited node only through
                // this edge
                if (discoveryTime[curNode->getId()] <
                    lowestDisc[neighborNode->getId()]) {
                  result.bridges.emplace_back(*edge);
                }
              }

              if ((typeMask & TARJAN_FIND_CUTV) && (nodeIsAdded == false)) {
                if (curNode->getId() == rootId) {
                  numSon++;
                  // a root node is a cut vertices only when it connects at
                  // least two connected components
                  if (numSon == 2) {
                    nodeIsAdded = true;
                    result.cutVertices.emplace_back(*curNode);
                  }
                } else {
                  if (discoveryTime[curNode->getId()] <=
                      lowestDisc[neighborNode->getId()]) {
                    nodeIsAdded = true;
                    result.cutVertices.emplace_back(*curNode);
                  }
                }
              }

              if (typeMask & TARJAN_FIND_VBCC) {
                if (discoveryTime[curNode->getId()] <=
                    lowestDisc[neighborNode->getId()]) {
                  // if current node is a cut vertice or the root node, the vbcc
                  // a vertice-biconnect-component which contains the neighbor
                  // node
                  std::vector<Node<T>> vbcc;
                  while (true) {
                    // pop a top node out of stack until
                    // the neighbor node has been poped out
                    Node<T> nodeAtTop = vbccNodeStack.top();
                    vbccNodeStack.pop();
                    vbcc.emplace_back(nodeAtTop);
                    if (nodeAtTop == *neighborNode) {
                      break;
                    }
                  }
                  vbcc.emplace_back(*curNode);
                  result.verticeBiconnectedComps.emplace_back(std::move(vbcc));
                }
              }
            } else if ((edge != prevEdge) &&
                       ((isDirected == false) ||
                        (inStack.count(neighborNode->getId())))) {
              // it's not allowed to go through the previous edge back
              // for a directed graph, it's also not allowed to visit
              // a node that is not in stack
              lowestDisc[curNode->getId()] =
                  std::min(lowestDisc[curNode->getId()],
                           lowestDisc[neighborNode->getId()]);
            }
          }
        }
        // find sccs for a undirected graph is very similar with
        // find ebccs for a directed graph
        if ((typeMask & TARJAN_FIND_SCC) || (typeMask & TARJAN_FIND_EBCC)) {
          std::stack<Node<T>> &nodeStack =
              (typeMask & TARJAN_FIND_SCC) ? sccNodeStack : ebccNodeStack;
          if (discoveryTime[curNode->getId()] == lowestDisc[curNode->getId()]) {
            std::vector<Node<T>> connectedComp;
            while (true) {
              // pop a top node out of stack until
              // the current node has been poped out
              Node<T> nodeAtTop = nodeStack.top();
              nodeStack.pop();
              if (typeMask & TARJAN_FIND_SCC) {
                inStack.erase(nodeAtTop.getId());
              }
              connectedComp.emplace_back(nodeAtTop);
              if (nodeAtTop == *curNode) {
                break;
              }
            }
            // store this component in result
            (typeMask & TARJAN_FIND_SCC)
                ? result.stronglyConnectedComps.emplace_back(
                      std::move(connectedComp))
                : result.edgeBiconnectedComps.emplace_back(
                      std::move(connectedComp));
          }
        }
      };

  for (const auto &node : nodeSet) {
    if (!discoveryTime.count(node->getId())) {
      rootId = node->getId();
      dfs_helper(node, nullptr);
    }
  }

  result.success = true;
  return result;
}

template <typename T>
TopoSortResult<T> Graph<T>::topologicalSort() const {
  TopoSortResult<T> result;
  result.success = false;

  if (!isDirectedGraph()) {
    result.errorMessage = ERR_UNDIR_GRAPH;
    return result;
  } else if (isCyclicDirectedGraphBFS()) {
    result.errorMessage = ERR_CYCLIC_GRAPH;
    return result;
  } else {
    const auto &adjMatrix = getAdjMatrix();
    const auto &nodeSet = getNodeSet();
    std::unordered_map<shared<const Node<T>>, bool, nodeHash<T>> visited;

    std::function<void(shared<const Node<T>>)> postorder_helper =
        [&postorder_helper, &adjMatrix, &visited,
         &result](shared<const Node<T>> curNode) {
          visited[curNode] = true;

          if (adjMatrix->find(curNode) != adjMatrix->end()) {
            for (const auto &edge : adjMatrix->at(curNode)) {
              const auto &nextNode = edge.first;
              if (false == visited[nextNode]) {
                postorder_helper(nextNode);
              }
            }
          }

          result.nodesInTopoOrder.push_back(*curNode);
        };

    auto numNodes = adjMatrix->size();
    result.nodesInTopoOrder.reserve(numNodes);

    for (const auto &node : nodeSet) {
      if (false == visited[node]) {
        postorder_helper(node);
      }
    }

    result.success = true;
    std::reverse(result.nodesInTopoOrder.begin(),
                 result.nodesInTopoOrder.end());
    return result;
  }
}

template <typename T>
TopoSortResult<T> Graph<T>::kahn() const {
  TopoSortResult<T> result;

  if (!isDirectedGraph()) {
    result.errorMessage = ERR_UNDIR_GRAPH;
    return result;
  } else {
    const auto adjMatrix = Graph<T>::getAdjMatrix();
    const auto nodeSet = Graph<T>::getNodeSet();
    result.nodesInTopoOrder.reserve(adjMatrix->size());

    std::unordered_map<CXXGraph::id_t, unsigned int> indegree;
    for (const auto &node : nodeSet) {
      indegree[node->getId()] = 0;
    }
    for (const auto &list : *adjMatrix) {
      auto children = list.second;
      for (const auto &child : children) {
        indegree[std::get<0>(child)->getId()]++;
      }
    }

    std::queue<shared<const Node<T>>> topologicalOrder;

    for (const auto &node : nodeSet) {
      if (!indegree[node->getId()]) {
        topologicalOrder.emplace(node);
      }
    }

    size_t visited = 0;
    while (!topologicalOrder.empty()) {
      shared<const Node<T>> currentNode = topologicalOrder.front();
      topologicalOrder.pop();
      result.nodesInTopoOrder.push_back(*currentNode);

      if (adjMatrix->find(currentNode) != adjMatrix->end()) {
        for (const auto &child : adjMatrix->at(currentNode)) {
          if (--indegree[std::get<0>(child)->getId()] == 0) {
            topologicalOrder.emplace(std::get<0>(child));
          }
        }
      }
      visited++;
    }

    if (visited != nodeSet.size()) {
      result.errorMessage = ERR_CYCLIC_GRAPH;
      result.nodesInTopoOrder.clear();
      return result;
    }

    result.success = true;
    return result;
  }
}

template <typename T>
SCCResult<T> Graph<T>::kosaraju() const {
  SCCResult<T> result;
  result.success = false;

  if (!isDirectedGraph()) {
    result.errorMessage = ERR_UNDIR_GRAPH;
    return result;
  } else {
    auto nodeSet = getNodeSet();
    const auto adjMatrix = getAdjMatrix();
    // created visited map
    std::unordered_map<CXXGraph::id_t, bool> visited;
    for (const auto &node : nodeSet) {
      visited[node->getId()] = false;
    }

    std::stack<shared<const Node<T>>> st;
    std::function<void(shared<const Node<T>>)> dfs_helper =
        [this, &adjMatrix, &visited, &dfs_helper,
         &st](shared<const Node<T>> source) {
          // mark the vertex visited
          visited[source->getId()] = true;

          // travel the neighbors
          for (int i = 0; i < (*adjMatrix)[source].size(); i++) {
            shared<const Node<T>> neighbor = (*adjMatrix)[source].at(i).first;
            if (visited[neighbor->getId()] == false) {
              // make recursive call from neighbor
              dfs_helper(neighbor);
            }
          }

          st.push(source);
        };

    for (const auto &node : nodeSet) {
      if (visited[node->getId()] == false) {
        dfs_helper(node);
      }
    }

    // construct the transpose of the given graph
    AdjacencyMatrix<T> rev;
    auto addElementToAdjMatrix = [&rev](shared<const Node<T>> nodeFrom,
                                        shared<const Node<T>> nodeTo,
                                        shared<const Edge<T>> edge) {
      std::pair<shared<const Node<T>>, shared<const Edge<T>>> elem = {nodeTo,
                                                                      edge};
      rev[nodeFrom].push_back(std::move(elem));
    };
    for (const auto &edgeSetIt : edgeSet) {
      shared<const DirectedEdge<T>> d_edge =
          std::static_pointer_cast<const DirectedEdge<T>>(edgeSetIt);
      // Add the reverse edge to the reverse adjacency matrix
      addElementToAdjMatrix(d_edge->getNodePair().second,
                            d_edge->getNodePair().first, d_edge);
    }

    visited.clear();

    std::function<void(shared<const Node<T>>, std::vector<Node<T>> &)>
        dfs_helper1 =
            [this, &rev, &visited, &dfs_helper1](shared<const Node<T>> source,
                                                 std::vector<Node<T>> &comp) {
              // mark the vertex visited
              visited[source->getId()] = true;
              // Add the current vertex to the strongly connected
              // component
              comp.push_back(*source);

              // travel the neighbors
              for (int i = 0; i < rev[source].size(); i++) {
                shared<const Node<T>> neighbor = rev[source].at(i).first;
                if (visited[neighbor->getId()] == false) {
                  // make recursive call from neighbor
                  dfs_helper1(neighbor, comp);
                }
              }
            };

    while (st.size() != 0) {
      auto rem = st.top();
      st.pop();
      if (visited[rem->getId()] == false) {
        std::vector<Node<T>> comp;
        dfs_helper1(rem, comp);
        result.stronglyConnectedComps.push_back(comp);
      }
    }

    result.success = true;
    return result;
  }
}

template <typename T>
const DialResult Graph<T>::dial(const Node<T> &source, int maxWeight) const {
  DialResult result;
  result.success = false;

  const auto adj = getAdjMatrix();
  auto nodeSet = getNodeSet();

  auto source_node_it = std::find_if(
      nodeSet.begin(), nodeSet.end(),
      [&source](auto node) { return node->getUserId() == source.getUserId(); });
  if (source_node_it == nodeSet.end()) {
    // check if source node exist in the graph
    result.errorMessage = ERR_SOURCE_NODE_NOT_IN_GRAPH;
    return result;
  }
  /* With each distance, iterator to that vertex in
          its bucket is stored so that vertex can be deleted
          in O(1) at time of updation. So
          dist[i].first = distance of ith vertex from src vertex
          dits[i].second = vertex i in bucket number */
  auto V = nodeSet.size();
  std::unordered_map<shared<const Node<T>>,
                     std::pair<long, shared<const Node<T>>>, nodeHash<T>>
      dist;

  // Initialize all distances as infinite (INF)
  for (const auto &node : nodeSet) {
    dist[node].first = std::numeric_limits<long>::max();
  }

  // Create buckets B[].
  // B[i] keep vertex of distance label i
  std::vector<std::deque<shared<const Node<T>>>> B((maxWeight * V + 1));

  B[0].push_back(*source_node_it);
  dist[*source_node_it].first = 0;

  int idx = 0;
  while (true) {
    // Go sequentially through buckets till one non-empty
    // bucket is found
    while (B[idx].size() == 0u && idx < maxWeight * V) {
      idx++;
    }

    // If all buckets are empty, we are done.
    if (idx == maxWeight * V) {
      break;
    }

    // Take top vertex from bucket and pop it
    auto u = B[idx].front();
    B[idx].pop_front();

    // Process all adjacents of extracted vertex 'u' and
    // update their distanced if required.
    for (const auto &i : (*adj)[u]) {
      auto v = i.first;
      int weight = 0;
      if (i.second->isWeighted().has_value() &&
          i.second->isWeighted().value()) {
        if (i.second->isDirected().has_value() &&
            i.second->isDirected().value()) {
          shared<const DirectedWeightedEdge<T>> dw_edge =
              std::static_pointer_cast<const DirectedWeightedEdge<T>>(i.second);
          weight = (int)dw_edge->getWeight();
        } else if (i.second->isDirected().has_value() &&
                   !i.second->isDirected().value()) {
          shared<const UndirectedWeightedEdge<T>> udw_edge =
              std::static_pointer_cast<const UndirectedWeightedEdge<T>>(
                  i.second);
          weight = (int)udw_edge->getWeight();
        } else {
          // ERROR it shouldn't never returned ( does not exist a Node
          // Weighted and not Directed/Undirected)
          result.errorMessage = ERR_NO_DIR_OR_UNDIR_EDGE;
          return result;
        }
      } else {
        // No Weighted Edge
        result.errorMessage = ERR_NO_WEIGHTED_EDGE;
        return result;
      }
      auto u_i = std::find_if(
          dist.begin(), dist.end(),
          [u](std::pair<shared<const Node<T>>,
                        std::pair<long, shared<const Node<T>>>> const &it) {
            return (*u == *(it.first));
          });

      auto v_i = std::find_if(
          dist.begin(), dist.end(),
          [v](std::pair<shared<const Node<T>>,
                        std::pair<long, shared<const Node<T>>>> const &it) {
            return (*v == *(it.first));
          });
      long du = u_i->second.first;
      long dv = v_i->second.first;

      // If there is shorted path to v through u.
      if (dv > du + weight) {
        // If dv is not INF then it must be in B[dv]
        // bucket, so erase its entry using iterator
        // in O(1)
        if (dv != std::numeric_limits<long>::max()) {
          auto findIter = std::find(B[dv].begin(), B[dv].end(), dist[v].second);
          B[dv].erase(findIter);
        }

        //  updating the distance
        dist[v].first = du + weight;
        dv = dist[v].first;

        // pushing vertex v into updated distance's bucket
        B[dv].push_front(v);

        // storing updated iterator in dist[v].second
        dist[v].second = *(B[dv].begin());
      }
    }
  }
  for (const auto &dist_i : dist) {
    result.minDistanceMap[dist_i.first->getId()] = dist_i.second.first;
  }
  result.success = true;

  return result;
}

template <typename T>
double Graph<T>::fordFulkersonMaxFlow(const Node<T> &source,
                                      const Node<T> &target) const {
  if (!isDirectedGraph()) {
    return -1;
  }
  double maxFlow = 0;
  std::unordered_map<shared<const Node<T>>, shared<const Node<T>>, nodeHash<T>>
      parent;
  std::unordered_map<
      shared<const Node<T>>,
      std::unordered_map<shared<const Node<T>>, double, nodeHash<T>>,
      nodeHash<T>>
      weightMap;
  // build weight map
  auto edgeSet = this->getEdgeSet();
  for (const auto &edge : edgeSet) {
    // The Edge are all Directed at this point because is checked at the
    // start
    if (edge->isWeighted().value_or(false)) {
      shared<const DirectedWeightedEdge<T>> dw_edge =
          std::static_pointer_cast<const DirectedWeightedEdge<T>>(edge);
      weightMap[edge->getNodePair().first][edge->getNodePair().second] =
          dw_edge->getWeight();
    } else {
      weightMap[edge->getNodePair().first][edge->getNodePair().second] =
          0;  // No Weighted Edge are assumed to be 0 weigthed
    }
  }

  // Constuct iterators for source and target nodes in nodeSet
  auto nodeSet = getNodeSet();
  auto source_node_ptr = *std::find_if(
      nodeSet.begin(), nodeSet.end(),
      [&source](auto node) { return node->getUserId() == source.getUserId(); });
  auto target_node_ptr = *std::find_if(
      nodeSet.begin(), nodeSet.end(),
      [&target](auto node) { return node->getUserId() == target.getUserId(); });

  auto bfs_helper = [this, &source_node_ptr, &target_node_ptr, &parent,
                     &weightMap]() -> bool {
    std::unordered_map<shared<const Node<T>>, bool, nodeHash<T>> visited;
    std::queue<shared<const Node<T>>> queue;
    queue.push(source_node_ptr);
    visited[source_node_ptr] = true;
    parent[source_node_ptr] = nullptr;
    while (!queue.empty()) {
      auto u = queue.front();
      queue.pop();
      for (auto &v : weightMap[u]) {
        if (!visited[v.first] && v.second > 0) {
          queue.push(v.first);
          visited[v.first] = true;
          parent[v.first] = u;
        }
      }
    }

    return (visited[target_node_ptr]);
  };
  // Updating the residual values of edges
  while (bfs_helper()) {
    double pathFlow = std::numeric_limits<double>::max();
    for (auto v = target_node_ptr; v != source_node_ptr; v = parent[v]) {
      auto u = parent[v];
      pathFlow = std::min(pathFlow, weightMap[u][v]);
    }
    for (auto v = target_node_ptr; v != source_node_ptr; v = parent[v]) {
      auto u = parent[v];
      weightMap[u][v] -= pathFlow;
      weightMap[v][u] += pathFlow;
    }
    // Adding the path flows
    maxFlow += pathFlow;
  }

  return maxFlow;
}

template <typename T>
const std::vector<Node<T>> Graph<T>::graph_slicing(const Node<T> &start) const {
  std::vector<Node<T>> result;

  auto nodeSet = Graph<T>::getNodeSet();
  // check if start node in the graph
  auto start_node_it = std::find_if(
      nodeSet.begin(), nodeSet.end(),
      [&start](auto node) { return node->getUserId() == start.getUserId(); });
  if (start_node_it == nodeSet.end()) {
    return result;
  }
  std::vector<Node<T>> C = Graph<T>::depth_first_search(start);
  std::deque<shared<const Node<T>>> C1;  // complement of C i.e. nodeSet - C
  for (auto const &node : nodeSet) {
    // from the set of all nodes, remove nodes that exist in C
    if (std::find_if(C.begin(), C.end(), [node](const Node<T> nodeC) {
          return (*node == nodeC);
        }) == C.end())
      C1.push_back(node);
  }

  // For all nodes in C', apply DFS
  //  and get the list of reachable nodes and store in M
  std::vector<Node<T>> M;
  for (auto const &node : C1) {
    std::vector<Node<T>> reachableNodes = Graph<T>::depth_first_search(*node);
    M.insert(M.end(), reachableNodes.begin(), reachableNodes.end());
  }
  // removes nodes from C that are reachable from M.
  for (const auto &nodeC : C) {
    if (std::find_if(M.begin(), M.end(), [nodeC](const Node<T> nodeM) {
          return (nodeM == nodeC);
        }) == M.end())
      result.push_back(nodeC);
  }
  return result;
}
