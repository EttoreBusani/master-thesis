#ifndef PARTIAL_SOLUTION_CLASS_HPP
#define PARTIAL_SOLUTION_CLASS_HPP

#include"ResourceDistribution.hpp"
#include<iostream>
#include<vector>

template<typename T>
class PartialSolution : public ResourceDistribution<T> {
// Class template for Solution of relaxed optimization problem.
// In addition to resource allocation contains information about every user's
// deployment decision and the extra cost charged by the platform, as well as
// the total platform profit.
// Template parameter T is the same as of ResourceDistribution

public:
  // Map-keys are users ID and the value indicates if the corresponging user chooses
  // 1° or 2° deployment
  std::map<unsigned, bool> userSecondDepl;
  std::map<unsigned, bool> userFirstDepl;
  // Order in which the platform assigns applications in the edge
  std::vector<std::string> order={};

  // Extra cost charged by the platform for running 2° deployment
  double r=0;
  // Platform profit
  double profit=0;

  // Constructors
  PartialSolution(double r): r(r) {};
  PartialSolution()=default;
  virtual ~PartialSolution()=default;

  // Update Resources
  // Sets ResourceDistribution members of current object equal to the ones in input.
  // INPUTS: res -> ResourceDistribution object with new allocation of edge and cloud resources
  // OUTPUTS: -
  void updateResources(const ResourceDistribution<T>& res){
    this->edgeServers = res.edgeServers;
    this->edgeLoads = res.edgeLoads;
    this->cloudServers = res.cloudServers;
    this->cloudLoads = res.cloudLoads;
  }

  // Utility methods
  void print() const {
    std::cout << "Applications order: ";
    for(auto name : order) std::cout << " --> " << name;
    std::cout << std::endl;
    std::cout << "Resources Used: " << std::endl;
    for(auto name : order){
      std::cout << name << " uses :" << std::endl;
      std::cout << "  " << this->edgeServers.find(name)->second << " edge servers" << std::endl;
      std::cout << "  " << this->cloudServers.find(name)->second << " cloud servers" << std::endl;
      std::cout << "  Edge Load: " << this->edgeLoads.find(name)->second << std::endl;
      std::cout << "  Cloud Load: " << this->cloudLoads.find(name)->second << std::endl;
    }
    std::cout << "r: " << r << std::endl;
    std::cout << "Profit: " << profit << std::endl;
    std::cout << "\t" << std::endl;
  }
  void printLess() const {
    int count1=0;
    int count2=0;
    for(auto [key,val] : userFirstDepl){
      if(val) ++count1;
    }
    for(auto [key,val] : userSecondDepl){
      if(val) ++count2;
    }
    std::cout << "For r: " << r << std::endl;
    std::cout << " n. of users 1° depl: " << count1 << std::endl;
    std::cout << " n. of users 2° depl: " << count2 << std::endl;
    std::cout << "\t" << std::endl;
  };
};

// Operator <
// Defines order relation on Partial Solution based on profit
template<typename T>
inline bool operator<(const PartialSolution<T> & lhs,
const PartialSolution<T> & rhs){
  return lhs.profit < rhs.profit;
};




#endif
