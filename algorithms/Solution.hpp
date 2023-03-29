#ifndef SOLUTION_CLASS_HPP
#define SOLUTION_CLASS_HPP

#include"PartialSolution.hpp"
#include"Application.hpp"
#include<memory>

class Solution : public PartialSolution<int> {
// Class for the Solution of the platform Stackelberg Game (complete optimization problem).
// In addition to members of PartialSolution<int> class members, contains information
// about single user assignment to edge or cloud and load per time execution of each
// application.

public:

  // Map-keys are users ID and the value indicates if the corresponding user is
  // assigned to edge or cloud.
  std::map<unsigned, bool> userEdgeAllocation={};
  std::map<unsigned, bool> userCloudAllocation={};
  // Load per time execution-map of edge and cloud
  std::map<std::string,double> L_e={};
  std::map<std::string,double> L_c={};

  // Constructors
  Solution()=default;

  template<typename T>
  Solution(const PartialSolution<T>& partSol){
    r = partSol.r;
    order = partSol.order;
    userFirstDepl = partSol.userFirstDepl;
    userSecondDepl = partSol.userSecondDepl;
    edgeLoads = partSol.edgeLoads;
    cloudLoads = partSol.cloudLoads;
  };

  // Initialize Application Values
  // Sets the loads, servers and VMs corresponding to the application in input
  // to null values.
  // INPUTS: name -> string containing name of application whose values are going to be initialized
  // OUTPUTS: -
  void initializeAppValues(const std::shared_ptr<Application>& app);
  void print() const;
};

// Operator <
// Defines order relation based on solution profit.
bool operator< (const Solution& lhs, const Solution& rhs);

// Comparison Operator
// Defines order relation based on profit and other parameters. Contrary to
// operator <, two solutions having the same profit are not considered equal,
// unless all other members are also the same. It is used in Tabu Search
// to distinguish solutions with the same profit but different resource and users
// allocation.
bool comparisonOperator(const Solution&, const Solution&);


class TabuAttribute {
public:
  std::map<std::string,int> edgeServers={};
  std::map<std::string,int> cloudServers={};
  int edgeUsers=0;
  int cloudUsers=0;

  TabuAttribute(const Solution& sol){
    edgeServers = sol.edgeServers;
    cloudServers = sol.cloudServers;
    for(auto [id,val] : sol.userEdgeAllocation){
      if(val) ++edgeUsers;
    }
    for(auto [id,val] : sol.userCloudAllocation){
      if(val) ++cloudUsers;
    }
  };
  TabuAttribute() = default;

  bool operator==(const TabuAttribute& rhs) const{
    return (this->edgeServers==rhs.edgeServers && this->cloudServers==rhs.cloudServers
    && this->edgeUsers==rhs.edgeUsers && this->cloudUsers==rhs.cloudUsers);
  };
  bool compOp(const TabuAttribute& rhs) const{
    return (this->edgeServers==rhs.edgeServers && this->cloudServers==rhs.cloudServers);
  }
  void print() const{
    for(auto [name,n] : edgeServers) std::cout << name << ": n_e=" << n
    << ", n_c=" << cloudServers.find(name)->second << "\n";
    std::cout << "Users on edge: " << edgeUsers << ", users on cloud: " << cloudUsers << std::endl;
  };
};



#endif
