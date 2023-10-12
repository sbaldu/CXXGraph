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

#ifndef __CXXGRAPH_PARTITIONING_COORDINATEDRECORD_H__
#define __CXXGRAPH_PARTITIONING_COORDINATEDRECORD_H__

#pragma once

#include <set>

#include "Record.hpp"
#include "CXXGraph/Utility/Typedef.hpp"

namespace CXXGraph {
namespace Partitioning {

class CoordinatedRecord : public Record {
 private:
  std::set<int> partitions = {};
  std::mutex *lock = nullptr;
  int degree = 0;

 public:
  CoordinatedRecord();
  ~CoordinatedRecord();

  const std::set<int> &getPartitions() const override;
  void addPartition(const int m) override;
  bool hasReplicaInPartition(const int m) const override;
  bool getLock() override;
  bool releaseLock() override;
  int getReplicas() const override;
  int getDegree() const override;
  void incrementDegree() override;

  void addAll(const std::set<int> &set);
  std::set<int> partition_intersection(
      const std::shared_ptr<CoordinatedRecord> other) const;
  std::set<int> partition_union(
      const std::shared_ptr<CoordinatedRecord> other) const;
  std::set<int> partition_difference(
      const std::shared_ptr<CoordinatedRecord> other) const;
};

std::set<int> CoordinatedRecord::partition_intersection(
    std::shared_ptr<CoordinatedRecord> other) const {
  std::set<int> result;
  set_intersection(this->partitions.begin(), this->partitions.end(),
                   other->partitions.begin(), other->partitions.end(),
                   std::inserter(result, result.begin()));
  return result;
}

std::set<int> CoordinatedRecord::partition_union(
    std::shared_ptr<CoordinatedRecord> other) const {
  std::set<int> result;
  set_union(this->partitions.begin(), this->partitions.end(),
            other->partitions.begin(), other->partitions.end(),
            std::inserter(result, result.begin()));
  return result;
}

std::set<int> CoordinatedRecord::partition_difference(
    std::shared_ptr<CoordinatedRecord> other) const {
  std::set<int> result;
  set_difference(this->partitions.begin(), this->partitions.end(),
                 other->partitions.begin(), other->partitions.end(),
                 std::inserter(result, result.begin()));
  return result;
}

CoordinatedRecord::CoordinatedRecord() : partitions() {
  lock = new std::mutex();
  degree = 0;
}

CoordinatedRecord::~CoordinatedRecord() {
  // std::cout << "CoordinatedRecord::~CoordinatedRecord()" << std::endl;
  // TODOOOOOOOO
  if (lock != nullptr) {
    delete lock;
  }
}

const std::set<int> &CoordinatedRecord::getPartitions() const {
  return partitions;
}

void CoordinatedRecord::addPartition(int m) {
  if (m == -1) {
    std::cout << "ERROR! record.addPartition(-1)" << std::endl;
    exit(-1);
  }
  partitions.insert(m);
}

bool CoordinatedRecord::hasReplicaInPartition(const int m) const {
  return partitions.find(m) != partitions.end();
}

bool CoordinatedRecord::getLock() {
  return lock->try_lock();
}

bool CoordinatedRecord::releaseLock() {
  lock->unlock();
  return true;
}

int CoordinatedRecord::getReplicas() const {
  return (int)partitions.size();
}

int CoordinatedRecord::getDegree() const {
  return degree;
}

void CoordinatedRecord::incrementDegree() {
  degree++;
}

void CoordinatedRecord::addAll(const std::set<int> &set) {
  partitions.insert(set.begin(), set.end());
}
}  // namespace Partitioning
}  // namespace CXXGraph

#endif  // __CXXGRAPH_PARTITIONING_COORDINATEDRECORD_H__
