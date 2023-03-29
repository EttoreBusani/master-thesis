#include"Solution.hpp"

void Solution::print() const{
  std::cout << " \n Applications order: ";
  for(auto name : order) std::cout << " --> " << name;
  std::cout << std::endl;
  std::cout << "Resources Used: " << std::endl;
  for(auto name : order){
    std::cout << name << " uses :" << std::endl;
    std::cout << "  " << this->edgeServers.find(name)->second << " edge servers" << std::endl;
    std::cout << "  " << this->cloudServers.find(name)->second << " cloud servers" << std::endl;
  }
  int sumEdge = 0;
  int sumCloud = 0;
  for(auto [id, val] : userEdgeAllocation){
    if(val)
      sumEdge += 1;
  }
  for(auto [id, val] : userCloudAllocation){
    if(val)
      sumCloud += 1;
  }
  std::cout << "Users served on Edge: " << sumEdge << std::endl;
  std::cout << "Users served on Cloud: " << sumCloud << std::endl;
  std::cout << "r: " << r << std::endl;
  std::cout << "Profit: " << profit << std::endl;
  std::cout << "----------------------------------------------------------------" << std::endl;
}
// Initialize Application Values
void Solution::initializeAppValues(const std::shared_ptr<Application>& app){
  this->edgeServers[app->name]=0;
  this->cloudServers[app->name]=0;
  L_e[app->name]=app->D_e*edgeLoads.find(app->name)->second;
  L_c[app->name]=app->D_c*cloudLoads.find(app->name)->second;
}
// Operator <
bool operator<(const Solution& lhs, const Solution& rhs){
  if(lhs.profit == rhs.profit) return lhs.r > rhs.r;
  return lhs.profit < rhs.profit;
}
// Comparison Operator
bool comparisonOperator(const Solution& lhs, const Solution& rhs){
  for(auto it=lhs.userEdgeAllocation.cbegin();
  it!=lhs.userEdgeAllocation.cend(); ++it){
    if(it->second && !rhs.userEdgeAllocation.find(it->first)->second) return false;
    if(!it->second && rhs.userEdgeAllocation.find(it->first)->second) return true;
  }
  for(auto it=lhs.userCloudAllocation.cbegin();
  it!=lhs.userCloudAllocation.cend(); ++it){
    if(it->second && !rhs.userCloudAllocation.find(it->first)->second) return true;
    if(!it->second && rhs.userCloudAllocation.find(it->first)->second) return false;
  }
  return lhs.profit<rhs.profit;
}
