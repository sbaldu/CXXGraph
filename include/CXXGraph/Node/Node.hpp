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

#ifndef __CXXGRAPH_NODE_H__
#define __CXXGRAPH_NODE_H__

#pragma once

#include "CXXGraph/Utility/id_t.hpp"

#include <iomanip>
#include <iostream>

namespace CXXGraph {
class Node;
std::ostream &operator<<(std::ostream &os, const Node &node);

class Node {
 private:
  CXXGraph::id_t id = 0;
  std::string userId = "";
  void setId(const std::string &);

 public:
  Node(const std::string &);
  ~Node() = default;
  const CXXGraph::id_t &getId() const;
  const std::string &getUserId() const;
  // operator
  bool operator==(const Node &b) const;
  bool operator<(const Node &b) const;
  friend std::ostream &operator<<(std::ostream &os, const Node &node);
};

Node::Node(const std::string& id) {
  this->userId = id;
  // the userid is set as sha512 hash of the user provided id
  setId(id);
}

void Node::setId(const std::string &inpId) {
  // const unsigned char* userId = reinterpret_cast<const unsigned char
  // *>((*inpId).c_str() ); unsigned char obuf[64]; unsigned long long obuf[8];
  // SHA512(userId, (*inpId).length(), reinterpret_cast<unsigned char*>(obuf));
  /**
  // Transform byte-array to string
  std::stringstream shastr;
  shastr << std::hex << std::setfill('0');
  int i = 0;
  //unsigned long can only store 8 bytes so we truncate the hash to 8 bytes
  for (const auto &byte: obuf)
  {
          shastr << std::setw(2) << static_cast<int>(byte);
          i++;
          if (i==8) break;
  }
  auto idStr =  shastr.str();
  // convert hex string to unsigned long long
  std::istringstream iss(idStr);
  iss >> std::hex >> this->id;

  **/
  this->id = std::hash<std::string>{}(inpId);
}

const CXXGraph::id_t &Node::getId() const {
  return id;
}

const std::string &Node::getUserId() const {
  return userId;
}

// The data type T must have an overload of the equality operator
bool Node::operator==(const Node &b) const {
  return this->id == b.id;
}

bool Node::operator<(const Node &b) const {
  return this->id < b.id;
}

// ostream overload
// The data type T must have an overload of the ostream operator
std::ostream &operator<<(std::ostream &os, const Node &node) {
  os << "Node: {\n"
     << "  Id:\t" << node.userId << "\n}";
  return os;
}
}  // namespace CXXGraph

#endif  // __CXXGRAPH_NODE_H__
