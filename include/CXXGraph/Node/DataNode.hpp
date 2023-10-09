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

#ifndef __CXXGRAPH_DATANODE_H__
#define __CXXGRAPH_DATANODE_H__

#pragma once

#include "CXXGraph/Utility/id_t.hpp"

#include <iomanip>
#include <iostream>

#include "Node.hpp"

namespace CXXGraph {
class Node;
std::ostream &operator<<(std::ostream &os, const Node &node);

template <typename T>
class DataNode : public Node {
  private:
	T data;
  public:
	DataNode(const std::string &, const T& data);
	// Move constructor
	DataNode(const std::string &, T&& data) noexcept;
	
	void setData(T&& new_data);
	const T &getData() const;

	// operator
	bool operator==(const DataNode<T> &b) const;
	bool operator<(const DataNode<T> &b) const;
};

template <typename T>
DataNode<T>::DataNode(const std::string& id, const T& data) : Node(id) {
  this->data = data;
}

template <typename T>
DataNode<T>::DataNode(const std::string& id, T&& data) noexcept : Node(id) {
  this->data = std::move(data);
}

template <typename T>
const T &DataNode<T>::getData() const {
  return data;
}

template <typename T>
void DataNode<T>::setData(T&& new_data) {
  this->data = std::move(new_data);
}

// The data type T must have an overload of the equality operator
template <typename T>
bool DataNode<T>::operator==(const DataNode<T> &b) const {
  return (this->id == b.id && this->data == b.data);
}

template <typename T>
bool DataNode<T>::operator<(const DataNode<T> &b) const {
  return (this->id < b.id);
}

// ostream overload
// The data type T must have an overload of the ostream operator
template <typename T>
std::ostream &operator<<(std::ostream &os, const DataNode<T> &node) {
  os << "Node: {\n"
     << "  Id:\t" << node.userId << "\n  Data:\t" << node.data << "\n}";
  return os;
}
}  // namespace CXXGraph

#endif  // __CXXGRAPH_NODE_H__
