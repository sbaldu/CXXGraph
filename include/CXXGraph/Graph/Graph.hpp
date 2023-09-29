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

/// Class that implement the Graph. ( This class is not Thread Safe )
template <typename T>
class Graph {
 private:
  T_EdgeSet<T> edgeSet = {};
  T_NodeSet<T> isolatedNodesSet = {};

  shared<AdjacencyMatrix<T>> cachedAdjMatrix;

  // Private non-const getter for the set of nodes
  std::unordered_set<shared<Node<T>>, nodeHash<T>> nodeSet();

  std::optional<std::pair<std::string, char>> getExtenstionAndSeparator(
      InputOutputFormat format) const;
  void writeGraphToStream(std::ostream &oGraph, std::ostream &oNodeFeat,
                          std::ostream &oEdgeWeight, const char &sep,
                          bool writeNodeFeat, bool writeEdgeWeight) const;
  void readGraphFromStream(std::istream &iGraph, std::istream &iNodeFeat,
                           std::istream &iEdgeWeight, bool readNodeFeat,
                           bool readEdgeWeight);
  int writeToDot(const std::string &workingDir, const std::string &OFileName,
                 const std::string &graphName) const;
  int readFromDot(const std::string &workingDir, const std::string &fileName);
  void recreateGraph(
      std::unordered_map<CXXGraph::id_t, std::pair<std::string, std::string>>
          &edgeMap,
      std::unordered_map<CXXGraph::id_t, bool> &edgeDirectedMap,
      std::unordered_map<std::string, T> &nodeFeatMap,
      std::unordered_map<CXXGraph::id_t, double> &edgeWeightMap);

#ifdef WITH_COMPRESSION
  int compressFile(const std::string &inputFile,
                   const std::string &outputFile) const;
  int decompressFile(const std::string &inputFile,
                     const std::string &outputFile) const;
#endif

 public:
  Graph();
  Graph(const T_EdgeSet<T> &edgeSet);
  virtual ~Graph() = default;
  /**
   * \brief
   * Function that return the Edge set of the Graph
   * Note: No Thread Safe
   *
   * @returns a list of Edges of the graph
   *
   */
  virtual const T_EdgeSet<T> &getEdgeSet() const;
  /**
   * \brief
   * Function set the Edge Set of the Graph
   * Note: No Thread Safe
   *
   * @param edgeSet The Edge Set
   *
   */
  virtual void setEdgeSet(const T_EdgeSet<T> &edgeSet);
  /**
   * \brief
   * Function add an Edge to the Graph Edge Set
   * First check if a pointer to a node with the same userId has
   * already been added, and if not add it
   * Note: No Thread Safe
   *
   * @param edge The Edge to insert
   *
   */
  virtual void addEdge(const Edge<T> *edge);
  /**
   * \brief
   * Function add an Edge to the Graph Edge Set
   * First check if a pointer to a node with the same userId has
   * already been added, and if not add it
   * Note: No Thread Safe
   *
   * @param edge The Edge to insert
   *
   */
  virtual void addEdge(shared<const Edge<T>> edge);
  /**
   * \brief
   * Function to add a Node to the Graph Node Set
   * Note: No Thread Safe
   *
   * @param pointer to the node
   *
   */
  virtual void addNode(const Node<T> *node);
  /**
   * \brief
   * Function to add a Node to the Graph Node Set
   * Note: No Thread Safe
   *
   * @param shared pointer to the node
   *
   */
  virtual void addNode(shared<const Node<T>> node);
  /**
   * \brief
   * Function remove an Edge from the Graph Edge Set
   * Note: No Thread Safe
   *
   * @param edgeId The Edge Id to remove
   *
   */
  virtual void removeEdge(const CXXGraph::id_t edgeId);
  /**
   * \brief
   * Function to remove a Node from the Graph Node Set
   * Note: No Thread Safe
   *
   * @param edgeId The Edge Id to remove
   *
   */
  virtual void removeNode(const std::string &nodeUserId);
  /**
   * \brief
   * Finds the given edge defined by v1 and v2 within the graph.
   *
   * @param v1 The first vertex.
   * @param v2 The second vertex.
   * @param id The edge id if the edge is found. Otherwise set to 0.
   * @return True if the edge exists in the graph.
   */
  virtual bool findEdge(const Node<T> *v1, const Node<T> *v2,
                        CXXGraph::id_t &id) const;
  /**
   * \brief
   * Overload of findEdge which takes shared pointers as parameters
   *
   * @param v1 The first vertex.
   * @param v2 The second vertex.
   * @param id The edge id if the edge is found. Otherwise set to 0.
   * @return True if the edge exists in the graph.
   */
  virtual bool findEdge(shared<const Node<T>> v1, shared<const Node<T>> v2,
                        CXXGraph::id_t &id) const;
  /**
   * \brief
   * Function that return the Node Set of the Graph
   * Note: No Thread Safe
   *
   * @returns a list of Nodes of the graph
   *
   */
  virtual const T_NodeSet<T> getNodeSet() const;
  /**
   * \brief
   * Function that return the Set of isolated nodes
   * in the Graph
   * Note: No Thread Safe
   *
   * @returns a list of Nodes of the graph
   *
   */
  virtual const T_NodeSet<T> getIsolatedNodeSet() const;
  /**
   * \brief
   * Function that sets the data contained in a node
   *
   * @param nodeUserId The userId string of the node whose data is to be changes
   * @param data The new value for the node data
   *
   */
  virtual void setNodeData(const std::string &nodeUserId, T data);
  /**
   * \brief
   * Function that sets the data contained in every node of the graph
   *
   * @param dataMap Map of the userId of every node with its new data value
   *
   */
  virtual void setNodeData(std::map<std::string, T> &dataMap);
  /**
   * \brief
   * Function that return an Edge with specific ID if Exist in the Graph
   * Note: No Thread Safe
   *
   * @param edgeId The Edge Id to return
   * @returns the Edge if exist
   *
   */
  virtual const std::optional<shared<const Edge<T>>> getEdge(
      const CXXGraph::id_t edgeId) const;
  /**
   * \brief
   * Function that return a Node with specific ID if Exist in the Graph
   * Note: No Thread Safe
   *
   * @param nodeId The Node Id to return
   * @returns the Node if exist
   *
   */
  virtual const std::optional<shared<const Node<T>>> getNode(
      const std::string &nodeUserId) const;
  /**
   * @brief This function generate a list of adjacency matrix with every element
   * of the matrix contain the node where is directed the link and the Edge
   * corrispondent to the link
   * Note: No Thread Safe
   */
  virtual shared<AdjacencyMatrix<T>> getAdjMatrix() const;

  virtual void cacheAdjMatrix();
  /**
   * \brief This function generates a set of nodes linked to the provided node
   * in a directed graph
   *
   * @param Pointer to the node
   *
   */
  virtual const std::unordered_set<shared<const Node<T>>, nodeHash<T>>
  outNeighbors(const Node<T> *node) const;
  /**
   * \brief This function generates a set of nodes linked to the provided node
   * in a directed graph
   *
   * @param Pointer to the node
   *
   */
  virtual const std::unordered_set<shared<const Node<T>>, nodeHash<T>>
  outNeighbors(shared<const Node<T>> node) const;
  /**
   * \brief This function generates a set of nodes linked to the provided node
   * in any graph
   *
   * @param Pointer to the node
   *
   */
  virtual const std::unordered_set<shared<const Node<T>>, nodeHash<T>>
  inOutNeighbors(const Node<T> *node) const;
  /**
   * \brief
   * \brief This function generates a set of nodes linked to the provided node
   * in any graph
   *
   * @param Pointer to the node
   *
   */
  virtual const std::unordered_set<shared<const Node<T>>, nodeHash<T>>
  inOutNeighbors(shared<const Node<T>> node) const;
  /**
   * \brief
   * \brief This function generates a set of Edges going out of a node
   * in any graph
   *
   * @param Pointer to the node
   *
   */
  virtual const std::unordered_set<shared<const Edge<T>>, edgeHash<T>> outEdges(
      const Node<T> *node) const;
  /**
   * \brief
   * \brief This function generates a set of Edges going out of a node
   * in any graph
   *
   * @param Shared pointer to the node
   *
   */
  virtual const std::unordered_set<shared<const Edge<T>>, edgeHash<T>> outEdges(
      shared<const Node<T>> node) const;
  /**
   * \brief
   * \brief This function generates a set of Edges coming in or going out of
   * a node in any graph
   *
   * @param Pointer to the node
   *
   */
  virtual const std::unordered_set<shared<const Edge<T>>, edgeHash<T>>
  inOutEdges(const Node<T> *node) const;
  /**
   * \brief
   * \brief This function generates a set of Edges coming in or going out of
   * a node in any graph
   *
   * @param Shared pointer to the node
   *
   */
  virtual const std::unordered_set<shared<const Edge<T>>, edgeHash<T>>
  inOutEdges(shared<const Node<T>> node) const;
  /**
   * @brief This function finds the subset of given a nodeId
   * Subset is stored in a map where keys are the hash-id of the node & values
   * is the subset.
   * @param subset query subset, we want to find target in this subset
   * @param elem elem that we wish to find in the subset
   *
   * @return parent node of elem
   * Note: No Thread Safe
   */
  virtual CXXGraph::id_t setFind(std::unordered_map<CXXGraph::id_t, Subset> *,
                                 const CXXGraph::id_t elem) const;
  /**
   * @brief This function finds the subset of given a nodeId
   * Subset is stored in a map where keys are the hash-id of the node & values
   * is the subset.
   * @param shared pointer to subset query subset, we want to find target in
   * this subset
   * @param elem elem that we wish to find in the subset
   *
   * @return parent node of elem
   * Note: No Thread Safe
   */
  virtual CXXGraph::id_t setFind(
      shared<std::unordered_map<CXXGraph::id_t, Subset>>,
      const CXXGraph::id_t elem) const;
  /**
   * @brief This function modifies the original subset array
   * such that it the union of two sets a and b
   * @param subset original subset is modified to obtain union of a & b
   * @param a parent id of set1
   * @param b parent id of set2
   * NOTE: Original subset is no longer available after union.
   * Note: No Thread Safe
   */
  virtual void setUnion(std::unordered_map<CXXGraph::id_t, Subset> *,
                        const CXXGraph::id_t set1,
                        const CXXGraph::id_t elem2) const;
  /**
   * @brief This function modifies the original subset array
   * such that it the union of two sets a and b
   * @param subset original subset is modified to obtain union of a & b
   * @param a parent id of set1
   * @param b parent id of set2
   * NOTE: Original subset is no longer available after union.
   * Note: No Thread Safe
   */
  virtual void setUnion(shared<std::unordered_map<CXXGraph::id_t, Subset>>,
                        const CXXGraph::id_t set1,
                        const CXXGraph::id_t elem2) const;
  /**
   * @brief This function finds the eulerian path of a directed graph using
   * hierholzers algorithm
   *
   * @return a vector containing nodes in eulerian path
   * Note: No Thread Safe
   */
  virtual std::shared_ptr<std::vector<Node<T>>> eulerianPath() const;
  /**
   * @brief Function runs the dijkstra algorithm for some source node and
   * target node in the graph and returns the shortest distance of target
   * from the source.
   * Note: No Thread Safe
   *
   * @param source source vertex
   * @param target target vertex
   *
   * @return shortest distance if target is reachable from source else ERROR in
   * case if target is not reachable from source or there is error in the
   * computation.
   */
  virtual const DijkstraResult dijkstra(const Node<T> &source,
                                        const Node<T> &target) const;
  /**
   * @brief This function runs the tarjan algorithm and returns different types
   * of results depending on the input parameter typeMask.
   *
   * @param typeMask each bit of typeMask within valid range represents a kind
   * of results should be returned.
   *
   * Note: No Thread Safe
   *
   * @return The types of return include strongly connected components
   * (only for directed graphs) and cut vertices、 bridges、edge
   * biconnected components and vertice biconnected components
   * (only for undirected graphs).
   */
  virtual const TarjanResult<T> tarjan(const unsigned int typeMask) const;
  /**
   * @brief Function runs the bellman-ford algorithm for some source node and
   * target node in the graph and returns the shortest distance of target
   * from the source. It can also detect if a negative cycle exists in the
   * graph. Note: No Thread Safe
   *
   * @param source source vertex
   * @param target target vertex
   *
   * @return shortest distance if target is reachable from source else ERROR in
   * case if target is not reachable from source. If there is no error then also
   * returns if the graph contains a negative cycle.
   */
  virtual const BellmanFordResult bellmanford(const Node<T> &source,
                                              const Node<T> &target) const;
  /**
   * @brief This function computes the transitive reduction of the graph,
   * returning a graph with the property of transitive closure satisfied. It
   * removes the "short-circuit" paths from a graph, leaving only the longest
   * paths. Commonly used to remove duplicate edges among nodes that do not pass
   * through the entire graph.
   * @return A copy of the current graph with the transitive closure property
   * satisfied.
   *
   */
  virtual const Graph<T> transitiveReduction() const;
  /**
   * @brief Function runs the floyd-warshall algorithm and returns the shortest
   * distance of all pair of nodes. It can also detect if a negative cycle
   * exists in the graph. Note: No Thread Safe
   * @return a map whose keys are node ids and values are the shortest distance.
   * If there is no error then also returns if the graph contains a negative
   * cycle.
   */
  virtual const FWResult floydWarshall() const;
  /**
   * @brief Function runs the prim algorithm and returns the minimum spanning
   * tree if the graph is undirected. Note: No Thread Safe
   * @return a vector containing id of nodes in minimum spanning tree & cost of
   * MST
   */
  virtual const MstResult prim() const;
  /**
   * @brief Function runs the boruvka algorithm and returns the minimum spanning
   * tree & cost if the graph is undirected. Note: No Thread Safe
   * @return struct of type MstResult with following fields
   * success: true if algorithm completed successfully ELSE false
   * mst: vector containing id of nodes in minimum spanning tree & cost of MST
   * mstCost: Cost of MST
   * errorMessage: "" if no error ELSE report the encountered error
   */
  virtual const MstResult boruvka() const;
  /**
   * @brief Function runs the kruskal algorithm and returns the minimum spanning
   * tree if the graph is undirected. Note: No Thread Safe
   * @return struct of type MstResult with following fields
   * success: true if algorithm completed successfully ELSE false
   * mst: vector containing id of nodes in minimum spanning tree & cost of MST
   * mstCost: Cost of MST
   * errorMessage: "" if no error ELSE report the encountered error
   */
  virtual const MstResult kruskal() const;
  /**
   * \brief
   * Function runs the best first search algorithm over the graph
   * using an evaluation function to decide which adjacent node is
   * most promising to explore
   * Note: No Thread Safe
   *
   * @param source source node
   * @param target target node
   * @returns a struct with a vector of Nodes if target is reachable else ERROR
   * in case if target is not reachable or there is an error in the computation.
   *
   */
  virtual BestFirstSearchResult<T> best_first_search(
      const Node<T> &source, const Node<T> &target) const;
  /**
   * \brief
   * Function performs the breadth first search algorithm over the graph
   * Note: No Thread Safe
   *
   * @param start Node from where traversing starts
   * @returns a vector of Node indicating which Node were visited during the
   * search.
   *
   */
  virtual const std::vector<Node<T>> breadth_first_search(
      const Node<T> &start) const;
  /**
   * \brief
   * The multithreaded version of breadth_first_search
   * It turns out to be two indepentent functions because of implemntation
   * differences
   *
   * @param start Node from where traversing starts
   * @param num_threads number of threads
   * @returns a vector of Node indicating which Node were visited during the
   * search.
   *
   */
  virtual const std::vector<Node<T>> concurrency_breadth_first_search(
      const Node<T> &start, size_t num_threads) const;
  /**
   * \brief
   * Function performs the depth first search algorithm over the graph
   * Note: No Thread Safe
   *
   * @param start Node from where traversing starts
   * @returns a vector of Node indicating which Node were visited during the
   * search.
   *
   */
  virtual const std::vector<Node<T>> depth_first_search(
      const Node<T> &start) const;

  /**
   * \brief
   * This function uses DFS to check for cycle in the graph.
   * Pay Attention, this function work only with directed Graph
   * Note: No Thread Safe
   *
   * @return true if a cycle is detected, else false. ( false is returned also
   * if the graph in indirected)
   */
  virtual bool isCyclicDirectedGraphDFS() const;

  /**
   * \brief
   * This function uses BFS to check for cycle in the graph.
   * Pay Attention, this function work only with directed Graph
   * Note: No Thread Safe
   *
   * @return true if a cycle is detected, else false. ( false is returned also
   * if the graph in indirected)
   */
  virtual bool isCyclicDirectedGraphBFS() const;

  /**
   * @brief
   * This function checks if the given set of edges
   * forms a cycle or not using union-find method.
   *
   * @return true if a cycle is detected, else false
   */
  virtual bool containsCycle(const T_EdgeSet<T> *) const;
  /**
   * @brief
   * This function checks if the given set of edges
   * forms a cycle or not using union-find method.
   *
   * @return true if a cycle is detected, else false
   */
  virtual bool containsCycle(shared<const T_EdgeSet<T>>) const;
  /**
   * @brief
   * This function checks if the given Subset
   * forms a cycle or not using union-find method.
   *
   * @return true if a cycle is detected, else false
   */
  virtual bool containsCycle(
      shared<const T_EdgeSet<T>> edgeSet,
      shared<std::unordered_map<CXXGraph::id_t, Subset>>) const;

  /**
   * \brief
   * This function checks if a graph is directed
   * Note: No Thread Safe
   *
   * @return true if the graph is directed, else false.
   */
  virtual bool isDirectedGraph() const;

  /**
   * \brief
   * This function checks if a graph is undirected
   * Note: No Thread Safe
   *
   * @return true if the graph is undirected, else false.
   */
  virtual bool isUndirectedGraph() const;

  /**
   * @brief This function reverse the direction of the edges in a directed graph
   */
  virtual void reverseDirectedGraph();

  /**
   * @brief This function checks if the graph is connected or not
   * 	Applicable for Undirected Graph, for Directed Graph use the
   * isStronglyConnectedGraph() function
   *
   * @return true if the graph is connected
   * @return false otherwise
   */
  virtual bool isConnectedGraph() const;

  /**
   * @brief This function checks if the graph is strongly connected or not
   * 	Applicable for Directed Graph, for Undirected Graph use the
   * isConnectedGraph() function
   *
   * @return true if the graph is connected
   * @return false otherwise
   */
  virtual bool isStronglyConnectedGraph() const;

  /**
   * @brief This function sort nodes in topological order.
   * 	Applicable for Directed Acyclic Graph
   *
   * @return a struct with a vector of Nodes ordered topologically else ERROR in
   * case of undirected or cyclic graph
   */
  virtual TopoSortResult<T> topologicalSort() const;

  /**
   * @brief This function sort nodes in topological order using kahn's algorithm
   * 	Applicable for Directed Acyclic Graph
   *
   * @return a struct with a vector of Nodes ordered topologically else ERROR in
   * case of undirected or cyclic graph
   */
  virtual TopoSortResult<T> kahn() const;

  /**
  * \brief
  * This function performs performs the kosaraju algorthm on the graph to find
  the strongly connected components.
  *
  * Mathematical definition of the problem:
  * A strongly connected component (SCC) of a directed graph is a maximal
  strongly connected subgraph.

  * Note: No Thread Safe
  * @return a vector of vector of strongly connected components.
  */
  virtual SCCResult<T> kosaraju() const;

  /**
  * \brief
  * This function performs Graph Slicing based on connectivity
  *
  * Mathematical definition of the problem:
  *
  * Let G be the set of nodes in a graph and n be a given node in that set.
  * Let C be the non-strict subset of G containing both n and all nodes
  reachable
  * from n, and let C' be its complement. There's a third set M, which is the
  * non-strict subset of C containing all nodes that are reachable from any node
  in C'.
  * The problem consists of finding all nodes that belong to C but not to M.

  * Note: No Thread Safe
  * @param start Node from where traversing starts
  * @return a vector of nodes that belong to C but not to M.
  */
  virtual const std::vector<Node<T>> graph_slicing(const Node<T> &start) const;

  /**
   * @brief Function runs the Dial algorithm  (Optimized Dijkstra for small
   * range weights) for some source node and target node in the graph and
   * returns the shortest distance of target from the source. Note: No Thread
   * Safe
   *
   * @param source source vertex
   * @param maxWeight maximum weight of the edge
   *
   * @return shortest distance for all nodes reachable from source else ERROR in
   * case there is error in the computation.
   */
  virtual const DialResult dial(const Node<T> &source, int maxWeight) const;

  /**
   * @brief Function runs the Ford-Fulkerson algorithm for some source node and
   * target node in the graph and returns the maximum flow of the graph
   *
   * @param source source vertex
   * @param target  target vertex
   * @return double Max-Flow value or -1 in case of error
   */
  virtual double fordFulkersonMaxFlow(const Node<T> &source,
                                      const Node<T> &target) const;

  /**
   * \brief
   * This function writes the graph to an output file
   * Note: Not threadsafe
   *
   * @param format The output format of the file
   * @param workingDir The parent directory of the output file
   * @param OFileName The output filename
   * @param compress Enables compression (requires zlib)
   * @param writeNodeFeat Indicates if export also Node Features
   * @param writeEdgeWeight Indicates if export also Edge Weights
   * @return 0 if OK, else return a negative value
   */
  virtual int writeToFile(
      InputOutputFormat format = InputOutputFormat::STANDARD_CSV,
      const std::string &workingDir = ".",
      const std::string &OFileName = "graph", bool compress = false,
      bool writeNodeFeat = false, bool writeEdgeWeight = false) const;

  virtual int writeToDotFile(const std::string &workingDir,
                             const std::string &OFileName,
                             const std::string &graphName) const;

  virtual int writeToMTXFile(const std::string &workingDir,
                             const std::string &OFileName, char delimier) const;

  /**
   * \brief
   * This function reads the graph from an input file
   * Note: Not threadsafe
   *
   * @param format The input format of the file
   * @param workingDir The parent directory of the input
   * file
   * @param OFileName The input filename
   * @param compress Enables compression (requires zlib)
   * @param readNodeFeat Indicates if import also Node Features
   * @param readEdgeWeight Indicates if import also Edge Weights
   * @return 0 if OK, else return a negative value
   */
  virtual int readFromFile(
      InputOutputFormat format = InputOutputFormat::STANDARD_CSV,
      const std::string &workingDir = ".",
      const std::string &OFileName = "graph", bool compress = false,
      bool readNodeFeat = false, bool readEdgeWeight = false);

  virtual int readFromDotFile(const std::string &workingDir,
                              const std::string &fileName);

  virtual int readFromMTXFile(const std::string &workingDir,
                              const std::string &fileName);

  /**
   * \brief
   * This function partition a graph in a set of partitions
   * Note: No Thread Safe
   *
   * @param algorithm The partition algorithm
   * @param numberOfPartition The number of partitions
   * @return The partiton Map of the partitioned graph
   */
  virtual PartitionMap<T> partitionGraph(
      const Partitioning::PartitionAlgorithm algorithm,
      const unsigned int numberOfPartitions, const double param1 = 0.0,
      const double param2 = 0.0, const double param3 = 0.0,
      const unsigned int numberOfthreads =
          std::thread::hardware_concurrency()) const;

  friend std::ostream &operator<< <>(std::ostream &os, const Graph<T> &graph);
  friend std::ostream &operator<< <>(std::ostream &os,
                                     const AdjacencyMatrix<T> &adj);
};

template <typename T>
Graph<T>::Graph() {
  /* Caching the adjacency matrix */
  cacheAdjMatrix();
}

template <typename T>
Graph<T>::Graph(const T_EdgeSet<T> &edgeSet) {
  for (auto edgeIt : edgeSet) {
    this->edgeSet.insert(edgeIt);
  }
  /* Caching the adjacency matrix */
  cacheAdjMatrix();
}

template <typename T>
const T_EdgeSet<T> &Graph<T>::getEdgeSet() const {
  return edgeSet;
}

template <typename T>
void Graph<T>::setEdgeSet(const T_EdgeSet<T> &edgeSet) {
  this->edgeSet.clear();
  for (auto edgeIt : edgeSet) {
    this->edgeSet.insert(edgeIt);
  }
  /* Caching the adjacency matrix */
  cacheAdjMatrix();
}

template <typename T>
void Graph<T>::addEdge(const Edge<T> *edge) {
  if (edge->isDirected().has_value() && edge->isDirected().value()) {
    if (edge->isWeighted().has_value() && edge->isWeighted().value()) {
      auto edge_shared = make_shared<DirectedWeightedEdge<T>>(*edge);
      this->edgeSet.insert(edge_shared);

      std::pair<shared<const Node<T>>, shared<const Edge<T>>> elem = {
          edge_shared->getNodePair().second, edge_shared};
      (*cachedAdjMatrix)[edge_shared->getNodePair().first].push_back(
          std::move(elem));
    } else {
      auto edge_shared = make_shared<DirectedEdge<T>>(*edge);
      this->edgeSet.insert(edge_shared);

      std::pair<shared<const Node<T>>, shared<const Edge<T>>> elem = {
          edge_shared->getNodePair().second, edge_shared};
      (*cachedAdjMatrix)[edge_shared->getNodePair().first].push_back(
          std::move(elem));
    }
  } else {
    if (edge->isWeighted().has_value() && edge->isWeighted().value()) {
      auto edge_shared = make_shared<UndirectedWeightedEdge<T>>(*edge);
      this->edgeSet.insert(edge_shared);

      std::pair<shared<const Node<T>>, shared<const Edge<T>>> elem = {
          edge_shared->getNodePair().second, edge_shared};
      (*cachedAdjMatrix)[edge_shared->getNodePair().first].push_back(
          std::move(elem));
      std::pair<shared<const Node<T>>, shared<const Edge<T>>> elem1 = {
          edge_shared->getNodePair().first, edge_shared};
      (*cachedAdjMatrix)[edge_shared->getNodePair().second].push_back(
          std::move(elem1));
    } else {
      auto edge_shared = make_shared<UndirectedEdge<T>>(*edge);
      this->edgeSet.insert(edge_shared);

      std::pair<shared<const Node<T>>, shared<const Edge<T>>> elem = {
          edge_shared->getNodePair().second, edge_shared};
      (*cachedAdjMatrix)[edge_shared->getNodePair().first].push_back(
          std::move(elem));
      std::pair<shared<const Node<T>>, shared<const Edge<T>>> elem1 = {
          edge_shared->getNodePair().first, edge_shared};
      (*cachedAdjMatrix)[edge_shared->getNodePair().second].push_back(
          std::move(elem1));
    }
  }
}

template <typename T>
void Graph<T>::addEdge(shared<const Edge<T>> edge) {
  this->edgeSet.insert(edge);

  /* Adding new edge in cached adjacency matrix */
  if (edge.get()->isDirected().has_value() &&
      edge.get()->isDirected().value()) {
    std::pair<shared<const Node<T>>, shared<const Edge<T>>> elem = {
        edge.get()->getNodePair().second, edge};
    (*cachedAdjMatrix)[edge.get()->getNodePair().first].push_back(
        std::move(elem));
  } else {
    std::pair<shared<const Node<T>>, shared<const Edge<T>>> elem = {
        edge.get()->getNodePair().second, edge};
    (*cachedAdjMatrix)[edge.get()->getNodePair().first].push_back(
        std::move(elem));
    std::pair<shared<const Node<T>>, shared<const Edge<T>>> elem1 = {
        edge.get()->getNodePair().first, edge};
    (*cachedAdjMatrix)[edge.get()->getNodePair().second].push_back(
        std::move(elem1));
  }
}

template <typename T>
void Graph<T>::addNode(const Node<T> *node) {
  auto node_shared = make_shared<const Node<T>>(*node);
  this->isolatedNodesSet.insert(node_shared);
}

template <typename T>
void Graph<T>::addNode(shared<const Node<T>> node) {
  this->isolatedNodesSet.insert(node);
}

template <typename T>
void Graph<T>::removeEdge(const CXXGraph::id_t edgeId) {
  auto edgeOpt = Graph<T>::getEdge(edgeId);
  if (edgeOpt.has_value()) {
    /*
    edgeSet.erase(std::find_if(this->edgeSet.begin(), this->edgeSet.end(),
    [edgeOpt](const Edge<T> *edge) { return (*(edgeOpt.value()) == *edge); }));
    */
    edgeSet.erase(edgeSet.find(edgeOpt.value()));
    int delIndex = -1;
    int i = 0;
    /* Removing the edge from the cached adjacency matrix */
    for (auto elem :
         (*cachedAdjMatrix)[edgeOpt.value().get()->getNodePair().first]) {
      if (elem.second.get()->getId() == edgeId) {
        delIndex = i;
        break;
      }
      i++;
    }
    if (delIndex != -1) {
      (*cachedAdjMatrix)[edgeOpt.value().get()->getNodePair().first].erase(
          (*cachedAdjMatrix)[edgeOpt.value().get()->getNodePair().first]
              .begin() +
          delIndex);
    }

    delIndex = -1;
    i = 0;
    for (auto elem :
         (*cachedAdjMatrix)[edgeOpt.value().get()->getNodePair().second]) {
      if (elem.second.get()->getId() == edgeId) {
        delIndex = i;
        break;
      }
      i++;
    }
    if (delIndex != -1) {
      (*cachedAdjMatrix)[edgeOpt.value().get()->getNodePair().second].erase(
          (*cachedAdjMatrix)[edgeOpt.value().get()->getNodePair().second]
              .begin() +
          delIndex);
    }
  }
}

template <typename T>
void Graph<T>::removeNode(const std::string &nodeUserId) {
  auto nodeOpt = getNode(nodeUserId);
  auto isolatedNodeIt = isolatedNodesSet.find(nodeOpt.value());

  if (nodeOpt.has_value() && isolatedNodeIt != isolatedNodesSet.end()) {
    // The node is isolated
    isolatedNodesSet.erase(isolatedNodeIt);
  } else if (nodeOpt.has_value()) {
    // The node is not isolated
    // Remove the edges containing the node
    auto edgeset = edgeSet;
    for (auto edgeIt : edgeset) {
      if (edgeIt->getNodePair().first->getUserId() == nodeUserId ||
          edgeIt->getNodePair().second->getUserId() == nodeUserId) {
        this->removeEdge(edgeIt->getId());
      }
    }
  }
}

template <typename T>
bool Graph<T>::findEdge(const Node<T> *v1, const Node<T> *v2,
                        CXXGraph::id_t &id) const {
  auto v1_shared = make_shared<const Node<T>>(*v1);
  auto v2_shared = make_shared<const Node<T>>(*v2);

  return findEdge(v1_shared, v2_shared, id);
}

template <typename T>
bool Graph<T>::findEdge(shared<const Node<T>> v1, shared<const Node<T>> v2,
                        CXXGraph::id_t &id) const {
  // This could be made faster by looking for the edge hash, assuming we hash
  // based on node data, instead of a unique integer
  if (cachedAdjMatrix.get() != NULL && cachedAdjMatrix->size() != 0) {
    /* Searching for the edge using cached adjacency matrix */

    for (auto elem : (*cachedAdjMatrix)[v1]) {
      if (elem.first == v2) {
        id = elem.second.get()->getId();
        return true;
      }
    }
  } else {
    /* Searching for the edge using edgeset */

    for (auto e : this->edgeSet) {
      if ((e->getNodePair().first == v1) && (e->getNodePair().second == v2)) {
        id = e->getId();
        return true;
      }
      if (!e->isDirected() &&
          ((e->getNodePair().second == v1) && (e->getNodePair().first == v2))) {
        id = e->getId();
        return true;
      }
    }
  }
  /*for (auto e : this->edgeSet) {
    if ((e->getNodePair().first == v1) && (e->getNodePair().second == v2)) {
      id = e->getId();
      return true;
    }
    if (!e->isDirected() &&
        ((e->getNodePair().second == v1) && (e->getNodePair().first == v2))) {
      id = e->getId();
      return true;
    }
  }*/

  id = 0;
  return false;
}

template <typename T>
const T_NodeSet<T> Graph<T>::getNodeSet() const {
  T_NodeSet<T> nodeSet;

  for (const auto &edgeSetIt : edgeSet) {
    nodeSet.insert(edgeSetIt->getNodePair().first);
    nodeSet.insert(edgeSetIt->getNodePair().second);
  }
  // Merge with the isolated nodes
  nodeSet.insert(this->isolatedNodesSet.begin(), this->isolatedNodesSet.end());

  return nodeSet;
}

template <typename T>
const T_NodeSet<T> Graph<T>::getIsolatedNodeSet() const {
  return isolatedNodesSet;
}

template <typename T>
void Graph<T>::setNodeData(const std::string &nodeUserId, T data) {
  auto nodeSet = this->nodeSet();
  auto nodeIt = std::find_if(
      nodeSet.begin(), nodeSet.end(),
      [&nodeUserId](auto node) { return node->getUserId() == nodeUserId; });
  std::const_pointer_cast<Node<T>>(*nodeIt)->setData(std::move(data));
}

template <typename T>
void Graph<T>::setNodeData(std::map<std::string, T> &dataMap) {
  // Construct the set of all the nodes in the graph
  for (auto &nodeSetIt : this->nodeSet()) {
    nodeSetIt->setData(std::move(dataMap[nodeSetIt->getUserId()]));
  }
}

template <typename T>
const std::optional<shared<const Edge<T>>> Graph<T>::getEdge(
    const CXXGraph::id_t edgeId) const {
  for (const auto &it : edgeSet) {
    if (it->getId() == edgeId) {
      return it;
    }
  }

  return std::nullopt;
}

template <typename T>
const std::optional<shared<const Node<T>>> Graph<T>::getNode(
    const std::string &nodeUserId) const {
  for (const auto &it : getNodeSet()) {
    if (it->getUserId() == nodeUserId) {
      return it;
    }
  }

  return std::nullopt;
}

template <typename T>
std::unordered_set<shared<Node<T>>, nodeHash<T>> Graph<T>::nodeSet() {
  std::unordered_set<shared<Node<T>>, nodeHash<T>> nodeSet;
  for (auto &edgeSetIt : edgeSet) {
    nodeSet.insert(
        std::const_pointer_cast<Node<T>>(edgeSetIt->getNodePair().first));
    nodeSet.insert(
        std::const_pointer_cast<Node<T>>(edgeSetIt->getNodePair().second));
  }
  for (auto &isNodeIt : isolatedNodesSet) {
    nodeSet.insert(std::const_pointer_cast<Node<T>>(isNodeIt));
  }

  return nodeSet;
}

template <typename T>
std::optional<std::pair<std::string, char>> Graph<T>::getExtenstionAndSeparator(
    InputOutputFormat format) const {
  if (format == InputOutputFormat::STANDARD_CSV) {
    return std::pair<std::string, char>(".csv", ',');
  } else if (format == InputOutputFormat::STANDARD_TSV) {
    return std::pair<std::string, char>(".tsv", '\t');
  } else {
    return std::nullopt;
  }
}

template <typename T>
int Graph<T>::writeToDot(const std::string &workingDir,
                         const std::string &OFileName,
                         const std::string &graphName) const {
  const std::string linkSymbol = "--";
  const std::string directedLinkSymbol = "->";

  const std::string completePathToFileGraph =
      workingDir + '/' + OFileName + ".dot";
  std::ofstream ofileGraph;
  ofileGraph.open(completePathToFileGraph);
  if (!ofileGraph.is_open()) {
    // ERROR File Not Open
    return -1;
  }

  // Write the header of the DOT file
  std::string headerLine;
  if (this->isDirectedGraph()) {
    headerLine = "digraph " + graphName + " {\n";
  } else {
    headerLine = "graph " + graphName + " {\n";
  }
  ofileGraph << headerLine;

  for (auto const &edgePtr : edgeSet) {
    std::string edgeLine = "";
    if (edgePtr->isDirected().has_value() && edgePtr->isDirected().value()) {
      auto directedPtr =
          std::static_pointer_cast<const DirectedEdge<T>>(edgePtr);
      edgeLine += '\t' + directedPtr->getFrom().getUserId() + ' ';
      edgeLine += directedLinkSymbol + ' ';
      edgeLine += directedPtr->getTo().getUserId();
    } else {
      edgeLine += '\t' + edgePtr->getNodePair().first->getUserId() + ' ';
      edgeLine += linkSymbol + ' ';
      edgeLine += edgePtr->getNodePair().second->getUserId();
    }
    if (edgePtr->isWeighted().has_value() && edgePtr->isWeighted().value()) {
      // Weights in dot files must be integers
      edgeLine += " [weight=" +
                  std::to_string(static_cast<int>(
                      std::dynamic_pointer_cast<const Weighted>(edgePtr)
                          ->getWeight())) +
                  ']';
    }
    edgeLine += ";\n";
    ofileGraph << edgeLine;
  }
  ofileGraph << '}';
  ofileGraph.close();

  return 0;
}

// This ctype facet classifies ',' and '\t' as whitespace
struct csv_whitespace : std::ctype<char> {
  static const mask *make_table() {
    // make a copy of the "C" locale table
    static std::vector<mask> v(classic_table(), classic_table() + table_size);
    v[','] |= space;  // comma will be classified as whitespace
    v['\t'] |= space;
    v[' '] &= ~space;  // space will not be classified as whitespace
    return &v[0];
  }
  csv_whitespace(std::size_t refs = 0) : ctype(make_table(), false, refs) {}
};

template <typename T>
void Graph<T>::writeGraphToStream(std::ostream &oGraph, std::ostream &oNodeFeat,
                                  std::ostream &oEdgeWeight, const char &sep,
                                  bool writeNodeFeat,
                                  bool writeEdgeWeight) const {
  for (const auto &edge : edgeSet) {
    oGraph << edge->getId() << sep << edge->getNodePair().first->getUserId()
           << sep << edge->getNodePair().second->getUserId() << sep
           << ((edge->isDirected().has_value() && edge->isDirected().value())
                   ? 1
                   : 0)
           << std::endl;
  }

  if (writeNodeFeat) {
    auto nodeSet = getNodeSet();
    for (const auto &node : nodeSet) {
      oNodeFeat << node->getUserId() << sep << node->getData() << std::endl;
    }
  }

  if (writeEdgeWeight) {
    for (const auto &edge : edgeSet) {
      oEdgeWeight
          << edge->getId() << sep
          << (edge->isWeighted().has_value() && edge->isWeighted().value()
                  ? (std::dynamic_pointer_cast<const Weighted>(edge))
                        ->getWeight()
                  : 0.0)
          << sep
          << (edge->isWeighted().has_value() && edge->isWeighted().value() ? 1
                                                                           : 0)
          << std::endl;
    }
  }
}

template <typename T>
void Graph<T>::readGraphFromStream(std::istream &iGraph,
                                   std::istream &iNodeFeat,
                                   std::istream &iEdgeWeight, bool readNodeFeat,
                                   bool readEdgeWeight) {
  std::unordered_map<CXXGraph::id_t, std::pair<std::string, std::string>>
      edgeMap;
  std::unordered_map<CXXGraph::id_t, bool> edgeDirectedMap;
  std::unordered_map<std::string, T> nodeFeatMap;
  std::unordered_map<CXXGraph::id_t, double> edgeWeightMap;

  CXXGraph::id_t edgeId;
  std::string nodeId1;
  std::string nodeId2;
  bool directed;
  while (iGraph >> edgeId >> nodeId1 >> nodeId2 >>
         directed) { /* loop continually */
    edgeMap[edgeId] = std::pair<std::string, std::string>(nodeId1, nodeId2);
    edgeDirectedMap[edgeId] = directed;
  }

  if (readNodeFeat) {
    std::string nodeId;
    T nodeFeat;
    while (iNodeFeat >> nodeId >> nodeFeat) {
      nodeFeatMap[nodeId] = nodeFeat;
    }
  }

  if (readEdgeWeight) {
    CXXGraph::id_t edgeId;
    double weight;
    bool weighted;
    while (iEdgeWeight >> edgeId >> weight >> weighted) { /* loop continually */
      if (weighted) {
        edgeWeightMap[edgeId] = weight;
      }
    }
  }

  recreateGraph(edgeMap, edgeDirectedMap, nodeFeatMap, edgeWeightMap);
}

template <typename T>
int Graph<T>::readFromDot(const std::string &workingDir,
                          const std::string &fileName) {
  // Define the edge maps
  std::unordered_map<CXXGraph::id_t, std::pair<std::string, std::string>>
      edgeMap;
  std::unordered_map<std::string, T> nodeFeatMap;
  std::unordered_map<CXXGraph::id_t, bool> edgeDirectedMap;
  std::unordered_map<CXXGraph::id_t, double> edgeWeightMap;

  // Define the node strings and the "garbage collector" temp string
  std::string node1;
  std::string node2;
  std::string temp;

  // Get full path to the file and open it
  const std::string completePathToFileGraph =
      workingDir + '/' + fileName + ".dot";

  // Check if the graph is directed
  bool directed = false;
  std::ifstream fileContentStream(completePathToFileGraph);
  std::string fileContent(std::istreambuf_iterator<char>{fileContentStream},
                          {});
  if (fileContent.find("->") != std::string::npos) {
    directed = true;
  }
  // Check if the graph is weighted
  bool weighted = false;
  if (fileContent.find("weight") != std::string::npos) {
    weighted = true;
  }
  fileContentStream.close();

  std::ifstream iFile(completePathToFileGraph);
  // Write the header of the DOT file in the temp string
  getline(iFile, temp);

  CXXGraph::id_t edgeId = 0;
  std::string fileRow;
  while (getline(iFile, fileRow)) {
    // If you've reached the end of the file, stop
    if (fileRow == "}") {
      break;
    }

    // Remove the whitespaces before the definition of the edge
    while (*fileRow.begin() == ' ' || *fileRow.begin() == '\t') {
      fileRow.erase(fileRow.begin());
    }

    std::stringstream row_stream(fileRow);
    getline(row_stream, node1, ' ');
    // Store the symbol representing the edge inside temp
    getline(row_stream, temp, ' ');
    if (weighted) {
      getline(row_stream, node2, '[');
      // Remove any whitespaces or tabs from the node string
      node2.erase(std::remove(node2.begin(), node2.end(), ' '), node2.end());
      node2.erase(std::remove(node2.begin(), node2.end(), '\t'), node2.end());

      getline(row_stream, temp, '=');
      std::string weight;
      getline(row_stream, weight, ']');
      // Erase any whitespaces
      weight.erase(std::remove(weight.begin(), weight.end(), ' '),
                   weight.end());
      edgeWeightMap[edgeId] = std::stod(weight);
    } else {
      getline(row_stream, node2, ';');
    }

    // Save the edge and increment the edge counter
    edgeMap[edgeId] = std::pair<std::string, std::string>(node1, node2);
    edgeDirectedMap[edgeId] = directed;
    ++edgeId;
  }
  iFile.close();

  recreateGraph(edgeMap, edgeDirectedMap, nodeFeatMap, edgeWeightMap);
  return 0;
}

template <typename T>
void Graph<T>::recreateGraph(
    std::unordered_map<CXXGraph::id_t, std::pair<std::string, std::string>>
        &edgeMap,
    std::unordered_map<CXXGraph::id_t, bool> &edgeDirectedMap,
    std::unordered_map<std::string, T> &nodeFeatMap,
    std::unordered_map<CXXGraph::id_t, double> &edgeWeightMap) {
  std::unordered_map<std::string, shared<Node<T>>> nodeMap;
  for (const auto &edgeIt : edgeMap) {
    shared<Node<T>> node1(nullptr);
    shared<Node<T>> node2(nullptr);
    if (nodeMap.find(edgeIt.second.first) == nodeMap.end()) {
      // Create new Node
      T feat;
      if (nodeFeatMap.find(edgeIt.second.first) != nodeFeatMap.end()) {
        feat = std::move(nodeFeatMap.at(edgeIt.second.first));
      }
      node1 = make_shared<Node<T>>(edgeIt.second.first, feat);
      nodeMap[edgeIt.second.first] = node1;
    } else {
      node1 = nodeMap.at(edgeIt.second.first);
    }
    if (nodeMap.find(edgeIt.second.second) == nodeMap.end()) {
      // Create new Node
      T feat;
      if (nodeFeatMap.find(edgeIt.second.second) != nodeFeatMap.end()) {
        feat = std::move(nodeFeatMap.at(edgeIt.second.second));
      }
      node2 = make_shared<Node<T>>(edgeIt.second.second, feat);
      nodeMap[edgeIt.second.second] = node2;
    } else {
      node2 = nodeMap.at(edgeIt.second.second);
    }

    if (edgeWeightMap.find(edgeIt.first) != edgeWeightMap.end()) {
      if (edgeDirectedMap.find(edgeIt.first) != edgeDirectedMap.end() &&
          edgeDirectedMap.at(edgeIt.first)) {
        auto edge = make_shared<DirectedWeightedEdge<T>>(
            edgeIt.first, node1, node2, edgeWeightMap.at(edgeIt.first));
        addEdge(edge);
      } else {
        auto edge = make_shared<UndirectedWeightedEdge<T>>(
            edgeIt.first, node1, node2, edgeWeightMap.at(edgeIt.first));
        addEdge(edge);
      }
    } else {
      if (edgeDirectedMap.find(edgeIt.first) != edgeDirectedMap.end() &&
          edgeDirectedMap.at(edgeIt.first)) {
        auto edge = make_shared<DirectedEdge<T>>(edgeIt.first, node1, node2);
        addEdge(edge);
      } else {
        auto edge = make_shared<UndirectedEdge<T>>(edgeIt.first, node1, node2);
        addEdge(edge);
      }
    }
  }
}

#ifdef WITH_COMPRESSION
template <typename T>
int Graph<T>::compressFile(const std::string &inputFile,
                           const std::string &outputFile) const {
  std::ifstream ifs;
  ifs.open(inputFile);
  if (!ifs.is_open()) {
    // ERROR File Not Open
    return -1;
  }
  std::string content((std::istreambuf_iterator<char>(ifs)),
                      (std::istreambuf_iterator<char>()));

  const char *content_ptr = content.c_str();
  gzFile outFileZ = gzopen(outputFile.c_str(), "wb");
  if (outFileZ == NULL) {
    // printf("Error: Failed to gzopen %s\n", outputFile.c_str());
    return -1;
  }

  unsigned int zippedBytes;
  zippedBytes =
      gzwrite(outFileZ, content_ptr, (unsigned int)strlen(content_ptr));

  ifs.close();
  gzclose(outFileZ);
  return 0;
}

template <typename T>
int Graph<T>::decompressFile(const std::string &inputFile,
                             const std::string &outputFile) const {
  gzFile inFileZ = gzopen(inputFile.c_str(), "rb");
  if (inFileZ == NULL) {
    // printf("Error: Failed to gzopen %s\n", inputFile.c_str());
    return -1;
  }
  unsigned char unzipBuffer[8192];
  unsigned int unzippedBytes;
  std::vector<unsigned char> unzippedData;
  std::ofstream ofs;
  ofs.open(outputFile);
  if (!ofs.is_open()) {
    // ERROR File Not Open
    return -1;
  }
  while (true) {
    unzippedBytes = gzread(inFileZ, unzipBuffer, 8192);
    if (unzippedBytes > 0) {
      unzippedData.insert(unzippedData.end(), unzipBuffer,
                          unzipBuffer + unzippedBytes);
    } else {
      break;
    }
  }
  for (const auto &c : unzippedData) {
    ofs << c;
  }
  ofs << std::endl;
  ofs.close();
  gzclose(inFileZ);
  return 0;
}
#endif

template <typename T>
CXXGraph::id_t Graph<T>::setFind(
    std::unordered_map<CXXGraph::id_t, Subset> *subsets,
    const CXXGraph::id_t nodeId) const {
  auto subsets_ptr =
      make_shared<std::unordered_map<CXXGraph::id_t, Subset>>(*subsets);
  // find root and make root as parent of i
  // (path compression)
  if ((*subsets)[nodeId].parent != nodeId) {
    (*subsets)[nodeId].parent =
        Graph<T>::setFind(subsets_ptr, (*subsets)[nodeId].parent);
  }

  return (*subsets)[nodeId].parent;
}

template <typename T>
CXXGraph::id_t Graph<T>::setFind(
    shared<std::unordered_map<CXXGraph::id_t, Subset>> subsets,
    const CXXGraph::id_t nodeId) const {
  // find root and make root as parent of i
  // (path compression)
  if ((*subsets)[nodeId].parent != nodeId) {
    (*subsets)[nodeId].parent =
        Graph<T>::setFind(subsets, (*subsets)[nodeId].parent);
  }

  return (*subsets)[nodeId].parent;
}

template <typename T>
void Graph<T>::setUnion(std::unordered_map<CXXGraph::id_t, Subset> *subsets,
                        const CXXGraph::id_t elem1,
                        const CXXGraph::id_t elem2) const {
  /* auto subsets_ptr = make_shared<std::unordered_map<CXXGraph::id_t,
   * Subset>>(*subsets); */
  // if both sets have same parent
  // then there's nothing to be done
  /* if ((*subsets_ptr)[elem1].parent == (*subsets_ptr)[elem2].parent) return;
   */
  /* auto elem1Parent = Graph<T>::setFind(subsets_ptr, elem1); */
  /* auto elem2Parent = Graph<T>::setFind(subsets_ptr, elem2); */
  /* if ((*subsets_ptr)[elem1Parent].rank < (*subsets_ptr)[elem2Parent].rank) */
  /*   (*subsets_ptr)[elem1].parent = elem2Parent; */
  /* else if ((*subsets_ptr)[elem1Parent].rank >
   * (*subsets_ptr)[elem2Parent].rank) */
  /*   (*subsets_ptr)[elem2].parent = elem1Parent; */
  /* else { */
  /*   (*subsets_ptr)[elem2].parent = elem1Parent; */
  /*   (*subsets_ptr)[elem1Parent].rank++; */
  /* } */
  if ((*subsets)[elem1].parent == (*subsets)[elem2].parent) return;
  auto elem1Parent = Graph<T>::setFind(subsets, elem1);
  auto elem2Parent = Graph<T>::setFind(subsets, elem2);
  if ((*subsets)[elem1Parent].rank < (*subsets)[elem2Parent].rank)
    (*subsets)[elem1].parent = elem2Parent;
  else if ((*subsets)[elem1Parent].rank > (*subsets)[elem2Parent].rank)
    (*subsets)[elem2].parent = elem1Parent;
  else {
    (*subsets)[elem2].parent = elem1Parent;
    (*subsets)[elem1Parent].rank++;
  }
}

template <typename T>
void Graph<T>::setUnion(
    shared<std::unordered_map<CXXGraph::id_t, Subset>> subsets,
    const CXXGraph::id_t elem1, const CXXGraph::id_t elem2) const {
  // if both sets have same parent
  // then there's nothing to be done
  if ((*subsets)[elem1].parent == (*subsets)[elem2].parent) return;
  auto elem1Parent = Graph<T>::setFind(subsets, elem1);
  auto elem2Parent = Graph<T>::setFind(subsets, elem2);
  if ((*subsets)[elem1Parent].rank < (*subsets)[elem2Parent].rank)
    (*subsets)[elem1].parent = elem2Parent;
  else if ((*subsets)[elem1Parent].rank > (*subsets)[elem2Parent].rank)
    (*subsets)[elem2].parent = elem1Parent;
  else {
    (*subsets)[elem2].parent = elem1Parent;
    (*subsets)[elem1Parent].rank++;
  }
}

template <typename T>
std::shared_ptr<std::vector<Node<T>>> Graph<T>::eulerianPath() const {
  const auto nodeSet = Graph<T>::getNodeSet();
  const auto adj = Graph<T>::getAdjMatrix();

  std::shared_ptr<std::vector<Node<T>>> eulerPath =
      std::make_shared<std::vector<Node<T>>>();

  bool undirected = this->isUndirectedGraph();

  std::vector<shared<const Node<T>>> currentPath;
  // The starting node is the only node which has more outgoing than ingoing
  // links
  auto firstNodeIt =
      std::max_element(nodeSet.begin(), nodeSet.end(), [adj](auto n1, auto n2) {
        return adj->at(n1).size() < adj->at(n2).size();
      });
  auto currentNode = *(firstNodeIt);
  currentPath.push_back(currentNode);

  while (currentPath.size() > 0) {
    auto &edges = adj->at(currentNode);
    // we keep removing the edges that
    // have been traversed from the adjacency list
    if (edges.size()) {
      auto firstEdge = edges.back().second;

      shared<const Node<T>> nextNodeId;
      nextNodeId = firstEdge->getOtherNode(currentNode);

      currentPath.push_back(nextNodeId);
      currentNode = nextNodeId;
      edges.pop_back();
    } else {
      eulerPath->push_back(*currentNode);
      currentNode = currentPath.back();
      currentPath.pop_back();
    }
  }
  return eulerPath;
}

template <typename T>
shared<AdjacencyMatrix<T>> Graph<T>::getAdjMatrix() const {
  auto adj = std::make_shared<AdjacencyMatrix<T>>();
  auto addElementToAdjMatrix = [&adj](shared<const Node<T>> nodeFrom,
                                      shared<const Node<T>> nodeTo,
                                      shared<const Edge<T>> edge) {
    std::pair<shared<const Node<T>>, shared<const Edge<T>>> elem = {nodeTo,
                                                                    edge};
    (*adj)[nodeFrom].push_back(std::move(elem));
  };
  for (const auto &edgeSetIt : edgeSet) {
    if (edgeSetIt->isDirected().has_value() &&
        edgeSetIt->isDirected().value()) {
      shared<const DirectedEdge<T>> d_edge =
          std::static_pointer_cast<const DirectedEdge<T>>(edgeSetIt);
      addElementToAdjMatrix(d_edge->getNodePair().first,
                            d_edge->getNodePair().second, d_edge);
    } else if (edgeSetIt->isDirected().has_value() &&
               !edgeSetIt->isDirected().value()) {
      shared<const UndirectedEdge<T>> ud_edge =
          std::static_pointer_cast<const UndirectedEdge<T>>(edgeSetIt);
      ;
      addElementToAdjMatrix(ud_edge->getNodePair().first,
                            ud_edge->getNodePair().second, ud_edge);
      addElementToAdjMatrix(ud_edge->getNodePair().second,
                            ud_edge->getNodePair().first, ud_edge);
    } else {  // is a simple edge we cannot create adj matrix
      return adj;
    }
  }
  return adj;
}

template <typename T>
void Graph<T>::cacheAdjMatrix() {
  const auto adj = Graph<T>::getAdjMatrix();
  this->cachedAdjMatrix = adj;
}

template <typename T>
const std::unordered_set<shared<const Node<T>>, nodeHash<T>>
Graph<T>::outNeighbors(const Node<T> *node) const {
  auto node_shared = make_shared<const Node<T>>(*node);

  return outNeighbors(node_shared);
}

template <typename T>
const std::unordered_set<shared<const Node<T>>, nodeHash<T>>
Graph<T>::outNeighbors(shared<const Node<T>> node) const {
  auto adj = getAdjMatrix();
  if (adj->find(node) == adj->end()) {
    return std::unordered_set<shared<const Node<T>>, nodeHash<T>>();
  }
  auto nodeEdgePairs = adj->at(node);

  std::unordered_set<shared<const Node<T>>, nodeHash<T>> outNeighbors;
  for (auto pair : nodeEdgePairs) {
    if (pair.second->isDirected().has_value() &&
        pair.second->isDirected().value()) {
      outNeighbors.insert(pair.first);
    }
  }

  return outNeighbors;
}

template <typename T>
const std::unordered_set<shared<const Node<T>>, nodeHash<T>>
Graph<T>::inOutNeighbors(const Node<T> *node) const {
  auto node_shared = make_shared<const Node<T>>(*node);

  return inOutNeighbors(node_shared);
}

template <typename T>
const std::unordered_set<shared<const Node<T>>, nodeHash<T>>
Graph<T>::inOutNeighbors(shared<const Node<T>> node) const {
  auto adj = Graph<T>::getAdjMatrix();
  if (adj->find(node) == adj->end()) {
    return std::unordered_set<shared<const Node<T>>, nodeHash<T>>();
  }
  auto nodeEdgePairs = adj->at(node);

  std::unordered_set<shared<const Node<T>>, nodeHash<T>> inOutNeighbors;
  for (auto pair : nodeEdgePairs) {
    inOutNeighbors.insert(pair.first);
  }

  return inOutNeighbors;
}

template <typename T>
const std::unordered_set<shared<const Edge<T>>, edgeHash<T>> Graph<T>::outEdges(
    const Node<T> *node) const {
  auto node_shared = make_shared<const Node<T>>(*node);

  return outEdges(node_shared);
}

template <typename T>
const std::unordered_set<shared<const Edge<T>>, edgeHash<T>> Graph<T>::outEdges(
    shared<const Node<T>> node) const {
  if (cachedAdjMatrix->find(node) == cachedAdjMatrix->end()) {
    return std::unordered_set<shared<const Edge<T>>, edgeHash<T>>();
  }
  auto nodeEdgePairs = cachedAdjMatrix->at(node);

  std::unordered_set<shared<const Edge<T>>, edgeHash<T>> outEdges;
  for (auto pair : nodeEdgePairs) {
    if (pair.second->isDirected().has_value() &&
        pair.second->isDirected().value()) {
      outEdges.insert(pair.second);
    }
  }

  return outEdges;
}

template <typename T>
const std::unordered_set<shared<const Edge<T>>, edgeHash<T>>
Graph<T>::inOutEdges(const Node<T> *node) const {
  auto node_shared = make_shared<const Node<T>>(*node);

  return outEdges(node_shared);
}

template <typename T>
const std::unordered_set<shared<const Edge<T>>, edgeHash<T>>
Graph<T>::inOutEdges(shared<const Node<T>> node) const {
  if (cachedAdjMatrix->find(node) == cachedAdjMatrix->end()) {
    return std::unordered_set<shared<const Edge<T>>, edgeHash<T>>();
  }
  auto nodeEdgePairs = cachedAdjMatrix->at(node);

  std::unordered_set<shared<const Edge<T>>, edgeHash<T>> outEdges;
  for (auto pair : nodeEdgePairs) {
    outEdges.insert(pair.second);
  }

  return outEdges;
}

template <typename T>
int Graph<T>::writeToFile(InputOutputFormat format,
                          const std::string &workingDir,
                          const std::string &OFileName, bool compress,
                          bool writeNodeFeat, bool writeEdgeWeight) const {
  int result = 0;

  // Open streams and write
  auto extSep = getExtenstionAndSeparator(format);
  if (!extSep) {
    std::cerr << "Unknown format\n";
    return -1;
  }
  auto &[extension, separator] = *extSep;

  std::ofstream ofileGraph;
  std::ofstream ofileNodeFeat;
  std::ofstream ofileEdgeWeight;

  std::string completePathToFileGraph =
      workingDir + "/" + OFileName + extension;
  ofileGraph.open(completePathToFileGraph);
  if (!ofileGraph.is_open()) {
    // ERROR File Not Open
    return -1;
  }

  if (writeNodeFeat) {
    std::string completePathToFileNodeFeat =
        workingDir + "/" + OFileName + "_NodeFeat" + extension;
    ofileNodeFeat.open(completePathToFileNodeFeat);
    if (!ofileNodeFeat.is_open()) {
      // ERROR File Not Open
      return -1;
    }
  }

  if (writeEdgeWeight) {
    std::string completePathToFileEdgeWeight =
        workingDir + "/" + OFileName + "_EdgeWeight" + extension;
    ofileEdgeWeight.open(completePathToFileEdgeWeight);
    if (!ofileEdgeWeight.is_open()) {
      // ERROR File Not Open
      std::cout << "ERROR File Not Open" << std::endl;
      return -1;
    }
  }

  writeGraphToStream(ofileGraph, ofileNodeFeat, ofileEdgeWeight, separator,
                     writeNodeFeat, writeEdgeWeight);

  // Cleanup from writing
  ofileGraph.close();
  if (writeNodeFeat) ofileNodeFeat.close();
  if (writeEdgeWeight) ofileEdgeWeight.close();

#ifdef WITH_COMPRESSION
  if (result == 0 && compress) {
    auto _compress = [this, &workingDir, &OFileName, &writeNodeFeat,
                      &writeEdgeWeight](const std::string &extension) {
      std::string completePathToFileGraph =
          workingDir + "/" + OFileName + extension;
      std::string completePathToFileGraphCompressed =
          workingDir + "/" + OFileName + extension + ".gz";
      int _result = compressFile(completePathToFileGraph,
                                 completePathToFileGraphCompressed);
      if (_result == 0) {
        _result = remove(completePathToFileGraph.c_str());
      }
      if (_result == 0) {
        if (writeNodeFeat) {
          std::string completePathToFileNodeFeat =
              workingDir + "/" + OFileName + "_NodeFeat" + extension;
          std::string completePathToFileNodeFeatCompressed =
              workingDir + "/" + OFileName + "_NodeFeat" + extension + ".gz";
          _result = compressFile(completePathToFileNodeFeat,
                                 completePathToFileNodeFeatCompressed);
          if (_result == 0) {
            _result = remove(completePathToFileNodeFeat.c_str());
          }
        }
      }
      if (_result == 0) {
        if (writeEdgeWeight) {
          std::string completePathToFileEdgeWeight =
              workingDir + "/" + OFileName + "_EdgeWeight" + extension;
          std::string completePathToFileEdgeWeightCompressed =
              workingDir + "/" + OFileName + "_EdgeWeight" + extension + ".gz";
          _result = compressFile(completePathToFileEdgeWeight,
                                 completePathToFileEdgeWeightCompressed);
          if (_result == 0) {
            _result = remove(completePathToFileEdgeWeight.c_str());
          }
        }
      }
      return _result;
    };
    if (format == InputOutputFormat::STANDARD_CSV) {
      auto result = _compress(".csv");
    } else if (format == InputOutputFormat::STANDARD_TSV) {
      auto result = _compress(".tsv");
    } else {
      // OUTPUT FORMAT NOT RECOGNIZED
      result = -1;
    }
  }
#endif

  return result;
}

template <typename T>
int Graph<T>::readFromFile(InputOutputFormat format,
                           const std::string &workingDir,
                           const std::string &OFileName, bool compress,
                           bool readNodeFeat, bool readEdgeWeight) {
  int result = 0;

#ifdef WITH_COMPRESSION
  if (compress) {
    auto decompress = [this, &workingDir, &OFileName, &readNodeFeat,
                       &readEdgeWeight](const std::string &extension) {
      std::string completePathToFileGraph =
          workingDir + "/" + OFileName + extension;
      std::string completePathToFileGraphCompressed =
          workingDir + "/" + OFileName + extension + ".gz";
      int _result = decompressFile(completePathToFileGraphCompressed,
                                   completePathToFileGraph);
      if (_result == 0) {
        if (readNodeFeat) {
          std::string completePathToFileNodeFeat =
              workingDir + "/" + OFileName + "_NodeFeat" + extension;
          std::string completePathToFileNodeFeatCompressed =
              workingDir + "/" + OFileName + "_NodeFeat" + extension + ".gz";
          _result = decompressFile(completePathToFileNodeFeatCompressed,
                                   completePathToFileNodeFeat);
        }
      }
      if (_result == 0) {
        if (readEdgeWeight) {
          std::string completePathToFileEdgeWeight =
              workingDir + "/" + OFileName + "_EdgeWeight" + extension;
          std::string completePathToFileEdgeWeightCompressed =
              workingDir + "/" + OFileName + "_EdgeWeight" + extension + ".gz";
          _result = decompressFile(completePathToFileEdgeWeightCompressed,
                                   completePathToFileEdgeWeight);
        }
      }
      return _result;
    };
    if (format == InputOutputFormat::STANDARD_CSV) {
      result = decompress(".csv");
    } else if (format == InputOutputFormat::STANDARD_TSV) {
      result = decompress(".tsv");
    } else {
      // INPUT FORMAT NOT RECOGNIZED
      result = -1;
    }

    if (result != 0) {
      return result;
    }
  }
#endif
  // Open streams and read
  auto extSep = getExtenstionAndSeparator(format);
  if (!extSep) {
    std::cerr << "Unknown format\n";
    return -1;
  }
  auto &[extension, separator] = *extSep;

  std::string completePathToFileGraph =
      workingDir + "/" + OFileName + extension;
  std::string completePathToFileNodeFeat;
  std::string completePathToFileEdgeWeight;

  std::ifstream ifileGraph;
  std::ifstream ifileNodeFeat;
  std::ifstream ifileEdgeWeight;

  ifileGraph.open(completePathToFileGraph);
  if (!ifileGraph.is_open()) {
    // ERROR File Not Open
    // std::cout << "ERROR File Not Open : " << completePathToFileGraph <<
    // std::endl;
    return -1;
  }
  ifileGraph.imbue(std::locale(ifileGraph.getloc(), new csv_whitespace));

  if (readNodeFeat) {
    completePathToFileNodeFeat =
        workingDir + "/" + OFileName + "_NodeFeat" + extension;
    ifileNodeFeat.open(completePathToFileNodeFeat);
    if (!ifileNodeFeat.is_open()) {
      // ERROR File Not Open
      // std::cout << "ERROR File Not Open" << std::endl;
      return -1;
    }
    ifileNodeFeat.imbue(std::locale(ifileGraph.getloc(), new csv_whitespace));
  }

  if (readEdgeWeight) {
    completePathToFileEdgeWeight =
        workingDir + "/" + OFileName + "_EdgeWeight" + extension;
    ifileEdgeWeight.open(completePathToFileEdgeWeight);
    if (!ifileEdgeWeight.is_open()) {
      // ERROR File Not Open
      // std::cout << "ERROR File Not Open" << std::endl;
      return -1;
    }
    ifileEdgeWeight.imbue(std::locale(ifileGraph.getloc(), new csv_whitespace));
  }

  readGraphFromStream(ifileGraph, ifileNodeFeat, ifileEdgeWeight, readNodeFeat,
                      readEdgeWeight);

  // Cleanup
  ifileGraph.close();
#ifdef WITH_COMPRESSION
  if (compress) remove(completePathToFileGraph.c_str());
#endif

  if (readNodeFeat) {
    ifileNodeFeat.close();
#ifdef WITH_COMPRESSION
    if (compress) remove(completePathToFileNodeFeat.c_str());
#endif
  }

  if (readEdgeWeight) {
    ifileEdgeWeight.close();
#ifdef WITH_COMPRESSION
    if (compress) remove(completePathToFileEdgeWeight.c_str());
#endif
  }

  return result;
}

template <typename T>
int Graph<T>::writeToDotFile(const std::string &workingDir,
                             const std::string &OFileName,
                             const std::string &graphName) const {
  return writeToDot(workingDir, OFileName, graphName);
}

template <typename T>
int Graph<T>::writeToMTXFile(const std::string &workingDir,
                             const std::string &OFileName,
                             char delimitier) const {
  // Get the full path and open the file
  const std::string completePathToFileGraph =
      workingDir + '/' + OFileName + ".mtx";
  std::ofstream iFile(completePathToFileGraph);

  // Write the header of the file
  std::string header = "%%MatrixMarket graph";
  // Check if the adjacency matrix is symmetric, i.e., if all the edges are
  // undirected
  bool symmetric = !std::any_of(edgeSet.begin(), edgeSet.end(), [](auto edge) {
    return (edge->isDirected().has_value() && edge->isDirected().value());
  });
  // Write in the header whether the adj matrix is symmetric or not
  if (symmetric) {
    header += " symmetric\n";
  } else {
    header += '\n';
  }
  iFile << header;

  // Write the line containing the number of nodes and edges
  const std::string firstLine =
      std::to_string(getNodeSet().size()) + delimitier +
      std::to_string(getNodeSet().size()) + delimitier +
      std::to_string(getEdgeSet().size()) + '\n';
  iFile << firstLine;

  // Construct the edges
  for (const auto &edgeIt : edgeSet) {
    std::string line;
    line += edgeIt->getNodePair().first->getUserId() + delimitier;
    line += edgeIt->getNodePair().second->getUserId() + delimitier;
    if (edgeIt->isWeighted().has_value() && edgeIt->isWeighted().value()) {
      line += std::to_string(edgeIt->isWeighted().value()) + '\n';
    } else {
      line += std::to_string(1.) + '\n';
    }
    iFile << line;
  }

  iFile.close();
  return 0;
}

template <typename T>
int Graph<T>::readFromDotFile(const std::string &workingDir,
                              const std::string &fileName) {
  return readFromDot(workingDir, fileName);
}

template <typename T>
int Graph<T>::readFromMTXFile(const std::string &workingDir,
                              const std::string &fileName) {
  // Define the edge maps
  std::unordered_map<CXXGraph::id_t, std::pair<std::string, std::string>>
      edgeMap;
  std::unordered_map<std::string, T> nodeFeatMap;
  std::unordered_map<CXXGraph::id_t, bool> edgeDirectedMap;
  std::unordered_map<CXXGraph::id_t, double> edgeWeightMap;

  // Get full path to the file and open it
  const std::string completePathToFileGraph =
      workingDir + '/' + fileName + ".mtx";
  std::ifstream iFile(completePathToFileGraph);
  // Check that the file is open
  if (!iFile.is_open()) {
    return -1;
  }

  // Define the number of columns and rows in the matrix
  int n_cols, n_rows;
  int n_edges;
  bool undirected = false;

  // From the first line of the file read the number of rows, columns and edges
  std::string row_content;
  getline(iFile, row_content);
  if (row_content.find("symmetric") != std::string::npos) {
    undirected = true;
  }

  // Get rid of any commented lines between the header and the size line
  while (row_content.find('%') != std::string::npos) {
    getline(iFile, row_content);
  }

  // From the size line of the file read the number of rows, columns and edges
  std::stringstream header_stream(row_content);
  std::string value;
  getline(header_stream, value, ' ');
  n_rows = std::stoi(value);
  getline(header_stream, value, ' ');
  n_cols = std::stoi(value);
  getline(header_stream, value, ' ');
  n_edges = std::stoi(value);

  // Since the matrix represents the adjacency matrix, it must be square
  if (n_rows != n_cols) {
    return -1;
  }

  // Read the content of each line
  std::string node1;
  std::string node2;
  std::string edge_weight;
  CXXGraph::id_t edge_id = 0;
  while (getline(iFile, row_content)) {
    std::stringstream row_stream(row_content);

    // Read the content of the node ids and the weight into strings
    getline(row_stream, node1, ' ');
    getline(row_stream, node2, ' ');
    getline(row_stream, edge_weight);

    edgeMap[edge_id] = std::pair<std::string, std::string>(node1, node2);
    edgeWeightMap[edge_id] = std::stod(edge_weight);
    edgeDirectedMap[edge_id] = !undirected;

    // If the edge is a self-link, it must be undirected
    if (node1 == node2) {
      edgeDirectedMap[edge_id] = false;
    }

    // Increase the edge id
    ++edge_id;
  }

  if (n_edges != edgeMap.size()) {
    std::cout << "Error: The number of edges does not match the value provided "
                 "in the size line.\n";
    return -1;
  }

  iFile.close();
  recreateGraph(edgeMap, edgeDirectedMap, nodeFeatMap, edgeWeightMap);
  return 0;
}

template <typename T>
PartitionMap<T> Graph<T>::partitionGraph(
    const Partitioning::PartitionAlgorithm algorithm,
    const unsigned int numberOfPartitions, const double param1,
    const double param2, const double param3,
    const unsigned int numberOfThreads) const {
  PartitionMap<T> partitionMap;
  Partitioning::Globals globals(numberOfPartitions, algorithm, param1, param2,
                                param3, numberOfThreads);
  auto edgeSet_ptr = make_shared<const T_EdgeSet<T>>(getEdgeSet());
  globals.edgeCardinality = edgeSet_ptr->size();
  globals.vertexCardinality = this->getNodeSet().size();
  Partitioning::Partitioner<T> partitioner(edgeSet_ptr, globals);
  Partitioning::CoordinatedPartitionState<T> partitionState =
      partitioner.performCoordinatedPartition();
  partitionMap = partitionState.getPartitionMap();
  return partitionMap;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const Graph<T> &graph) {
  os << "Graph:\n";
  auto edgeList = graph.getEdgeSet();
  auto it = edgeList.begin();
  for (it; it != edgeList.end(); ++it) {
    if (!(*it)->isDirected().has_value() && !(*it)->isWeighted().has_value()) {
      // Edge Case
      os << **it << "\n";
    } else if (((*it)->isDirected().has_value() &&
                (*it)->isDirected().value()) &&
               ((*it)->isWeighted().has_value() &&
                (*it)->isWeighted().value())) {
      os << std::static_pointer_cast<const DirectedWeightedEdge<T>>(*it)
         << "\n";
    } else if (((*it)->isDirected().has_value() &&
                (*it)->isDirected().value()) &&
               !((*it)->isWeighted().has_value() &&
                 (*it)->isWeighted().value())) {
      os << std::static_pointer_cast<const DirectedEdge<T>>(*it) << "\n";
    } else if (!(it->isDirected().has_value() && it->isDirected().value()) &&
               (it->isWeighted().has_value() && it->isWeighted().value())) {
      os << std::static_pointer_cast<const UndirectedWeightedEdge<T>>(*it)
         << "\n";
    } else if (!(it->isDirected().has_value() && it->isDirected().value()) &&
               !(it->isWeighted().has_value() && it->isWeighted().value())) {
      os << std::static_pointer_cast<const UndirectedEdge<T>>(*it) << "\n";
    } else {
      os << *it << "\n";
    }
  }
  return os;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const AdjacencyMatrix<T> &adj) {
  os << "Adjacency Matrix:\n";
  unsigned long max_column = 0;
  for (const auto &it : adj) {
    if (it.second.size() > max_column) {
      max_column = (unsigned long)it.second.size();
    }
  }
  if (max_column == 0) {
    os << "ERROR in Print\n";
    return os;
  } else {
    os << "|--|";
    for (unsigned long i = 0; i < max_column; ++i) {
      os << "-----|";
    }
    os << "\n";
    for (const auto &it : adj) {
      os << "|N" << it.first->getId() << "|";
      for (const auto &it2 : it.second) {
        os << "N" << it2.first->getId() << ",E" << it2.second->getId() << "|";
      }
      os << "\n|--|";
      for (unsigned long i = 0; i < max_column; ++i) {
        os << "-----|";
      }
      os << "\n";
    }
  }
  return os;
}

}  // namespace CXXGraph
#endif  // __CXXGRAPH_GRAPH_H__
