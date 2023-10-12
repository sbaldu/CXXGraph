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

#ifndef __CXXGRAPH_PARTITION_H__
#define __CXXGRAPH_PARTITION_H__

#pragma once

#include <list>
#include <memory>
#include <unordered_set>

#include "PartitioningStats.hpp"
#include "CXXGraph/Utility/Typedef.hpp"

namespace CXXGraph {
// Smart pointers alias
template <typename T>
using unique = std::unique_ptr<T>;
template <typename T>
using shared= std::shared_ptr<T>;

using std::make_unique;
using std::make_shared;

class Graph;

using T_EdgeSet = std::unordered_set<shared<const Edge>, edgeHash>;
namespace Partitioning {
std::ostream &operator<<(std::ostream &o, const Partition &partition);

class Partition : public Graph {
 public:
  Partition();
  explicit Partition(const CXXGraph::id_t partitionId);
  explicit Partition(const T_EdgeSet &edgeSet);
  Partition(const CXXGraph::id_t partitionId, const T_EdgeSet &edgeSet);
  ~Partition() = default;
  /**
   * @brief Get the Partition ID
   *
   * @return The ID of the partition
   */
  CXXGraph::id_t getPartitionId() const;
  /**
   * @brief Set the Partition ID
   *
   * @param partitionId the ID to set
   */
  void setPartitionId(const CXXGraph::id_t partitionId);

 private:
  CXXGraph::id_t partitionId = 0;
};

/**
 * @brief Calculate and return the statistic of the Partitioned Graph
 *
 * @param partitionMap the Partition Map
 *
 * @return The Statistic of the Partioned Graph
 */
static PartitioningStats getPartitionStats(const PartitionMap &partitionMap);

/**
 * @brief Calculate the Maximum Load in a single partition (in terms of edges)
 * for the Partioned Graph
 *
 * @param partitionMap the Partition Map
 *
 * @return The value of the Maximum Load
 */
static unsigned int getMaxEdgesLoad(const PartitionMap &partitionMap);

/**
 * @brief Calculate the Minimum Load in a single partition (in terms of edges)
 * for the Partioned Graph
 *
 * @param partitionMap the Partition Map
 *
 * @return The value of the Minimum Load
 */
static unsigned int getMinEdgesLoad(const PartitionMap &partitionMap);

/**
 * @brief Calculate the Maximum Load in a single partition (in terms of nodes)
 * for the Partioned Graph
 *
 * @param partitionMap the Partition Map
 *
 * @return The value of the Maximum Load
 */
static unsigned int getMaxNodesLoad(const PartitionMap &partitionMap);

/**
 * @brief Calculate the Minimum Load in a single partition (in terms of nodes)
 * for the Partioned Graph
 *
 * @param partitionMap the Partition Map
 *
 * @return The value of the Minimum Load
 */
static unsigned int getMinNodesLoad(const PartitionMap &partitionMap);

/**
 * @brief Calculate the Number of Unique Edges in the Partitioned Graph ( this
 * value is equal to the number of edges in the Original Graph)
 *
 * @param partitionMap the Partition Map
 *
 * @return The number of Edges
 */
static unsigned int getNumberOfEdges(const PartitionMap &partitionMap);

/**
 * @brief Calculate the Number of Unique Nodes in the Partitioned Graph ( this
 * value is equal to the number of nodes in the Original Graph)
 *
 * @param partitionMap the Partition Map
 *
 * @return The number of Nodes
 */
static unsigned int getNumberOfNodes(const PartitionMap &partitionMap);

/**
 * @brief Calculate the Total Number of Edges in the Partitioned Graph
 *
 * @param partitionMap the Partition Map
 *
 * @return The number of Edges
 */
static unsigned int getNumberOfReplicatedEdges(
    const PartitionMap &partitionMap);

/**
 * @brief Calculate the Total Number of Nodes in the Partitioned Graph
 *
 * @param partitionMap the Partition Map
 *
 * @return The number of Nodes
 */
static unsigned int getNumberOfReplicatedNodes(
    const PartitionMap &partitionMap);

Partition::Partition() : Graph {
  partitionId = 0;
}

Partition::Partition(const CXXGraph::id_t partitionId) : Graph() {
  this->partitionId = partitionId;
}

Partition::Partition(const T_EdgeSet &edgeSet) : Graph(edgeSet) {
  partitionId = 0;
}

Partition::Partition(const CXXGraph::id_t partitionId,
                        const T_EdgeSet &edgeSet)
    : Graph(edgeSet) {
  this->partitionId = partitionId;
}

CXXGraph::id_t Partition::getPartitionId() const {
  return partitionId;
}

void Partition::setPartitionId(const CXXGraph::id_t partitionId) {
  this->partitionId = partitionId;
}

PartitioningStats getPartitionStats(const PartitionMap &partitionMap) {
  PartitioningStats result;
  result.numberOfPartitions = partitionMap.size();
  result.numberOfNodes = getNumberOfNodes(partitionMap);
  result.numberOfEdges = getNumberOfEdges(partitionMap);
  result.replicatedNodesCount = getNumberOfReplicatedNodes(partitionMap);
  result.replicatedEdgesCount = getNumberOfReplicatedEdges(partitionMap);
  result.maxEdgesLoad = getMaxEdgesLoad(partitionMap);
  result.minEdgesLoad = getMinEdgesLoad(partitionMap);
  result.maxNodesLoad = getMaxNodesLoad(partitionMap);
  result.minNodesLoad = getMinNodesLoad(partitionMap);
  result.edgesReplicationFactor =
      static_cast<double>(result.replicatedEdgesCount) / result.numberOfEdges;
  result.nodesReplicationFactor =
      static_cast<double>(result.replicatedNodesCount) / result.numberOfNodes;
  result.balanceEdgesFactor =
      static_cast<double>((result.maxEdgesLoad - result.minEdgesLoad)) /
      (result.maxEdgesLoad);
  result.balanceNodesFactor =
      static_cast<double>((result.maxNodesLoad - result.minNodesLoad)) /
      (result.maxNodesLoad);
  return result;
}

unsigned int getMaxEdgesLoad(const PartitionMap &partitionMap) {
  unsigned int maxLoad = 0;
  for (const auto &it : partitionMap) {
    if (it.second->getEdgeSet().size() > maxLoad) {
      maxLoad = it.second->getEdgeSet().size();
    }
  }
  return maxLoad;
}

unsigned int getMinEdgesLoad(const PartitionMap &partitionMap) {
  unsigned int minLoad = std::numeric_limits<unsigned int>::max();
  for (const auto &it : partitionMap) {
    if (it.second->getEdgeSet().size() < minLoad) {
      minLoad = it.second->getEdgeSet().size();
    }
  }
  return minLoad;
}

unsigned int getMaxNodesLoad(const PartitionMap &partitionMap) {
  unsigned int maxLoad = 0;
  for (const auto &it : partitionMap) {
    if (it.second->getNodeSet().size() > maxLoad) {
      maxLoad = it.second->getNodeSet().size();
    }
  }
  return maxLoad;
}

unsigned int getMinNodesLoad(const PartitionMap &partitionMap) {
  unsigned int minLoad = std::numeric_limits<unsigned int>::max();
  for (const auto &it : partitionMap) {
    if (it.second->getNodeSet().size() < minLoad) {
      minLoad = it.second->getNodeSet().size();
    }
  }
  return minLoad;
}

unsigned int getNumberOfEdges(const PartitionMap &partitionMap) {
  unsigned int numberOfEdges = 0;
  T_EdgeSet edgeSet;

  for (const auto &it : partitionMap) {
    const T_EdgeSet partitionEdgeSet = it.second->getEdgeSet();
    for (const auto &it2 : partitionEdgeSet) {
      edgeSet.insert(it2);
    }
  }

  return edgeSet.size();
}

unsigned int getNumberOfNodes(const PartitionMap &partitionMap) {
  unsigned int numberOfNodes = 0;
  std::unordered_set<shared<const Node>, nodeHash> nodeSet;

  for (const auto &it : partitionMap) {
    const std::unordered_set<shared<const Node>, nodeHash> partitionNodeSet = it.second->getNodeSet();
    for (const auto &it2 : partitionNodeSet) {
      // if (std::find_if(nodeSet.begin(), nodeSet.end(), [it2](const Node
      // *node)
      //                  { return (*it2 == *node); }) == nodeSet.end())
      // {
      nodeSet.insert(it2);
      // }
    }
  }

  return nodeSet.size();
}

unsigned int getNumberOfReplicatedEdges(const PartitionMap &partitionMap) {
  unsigned int numberOfEdges = 0;
  for (const auto &it : partitionMap) {
    numberOfEdges += it.second->getEdgeSet().size();
  }
  return numberOfEdges;
}

unsigned int getNumberOfReplicatedNodes(const PartitionMap &partitionMap) {
  unsigned int numberOfNodes = 0;
  for (const auto &it : partitionMap) {
    numberOfNodes += it.second->getNodeSet().size();
  }
  return numberOfNodes;
}

std::ostream &operator<<(std::ostream &os, const Partition &partition) {
  os << "Partition " << partition.getPartitionId() << ":\n";
  auto edgeList = partition.getEdgeSet();
  for (const auto &it : edgeList) {
    if (!(*it)->isDirected().has_value() && !(*it)->isWeighted().has_value()) {
      // Edge Case
      os << **it << "\n";
    } else if (((*it)->isDirected().has_value() &&
                (*it)->isDirected().value()) &&
               ((*it)->isWeighted().has_value() &&
                (*it)->isWeighted().value())) {
      os << std::static_pointer_cast<const DirectedWeightedEdge>(*it) << "\n";
    } else if ((it->isDirected().has_value() && it->isDirected().value()) &&
               !(it->isWeighted().has_value() && it->isWeighted().value())) {
      os << std::static_pointer_cast<const DirectedEdge>(*it) << "\n";
    } else if (!(it->isDirected().has_value() && it->isDirected().value()) &&
               (it->isWeighted().has_value() && it->isWeighted().value())) {
      os << std::static_pointer_cast<const UndirectedWeightedEdge>(*it) << "\n";
    } else if (!(it->isDirected().has_value() && it->isDirected().value()) &&
               !(it->isWeighted().has_value() && it->isWeighted().value())) {
      os << std::static_pointer_cast<const UndirectedEdge>(*it) << "\n";
    } else {
      // Should never happens
      os << "Wrong Edge Class"
         << "\n";
    }
  }
  return os;
}
}  // namespace Partitioning

}  // namespace CXXGraph

#endif  // __CXXGRAPH_PARTITION_H__
