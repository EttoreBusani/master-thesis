#ifndef RESOURCE_DISTRIBUTION_CLASS_HPP
#define RESOURCE_DISTRIBUTION_CLASS_HPP

#include<map>
#include<string>


template<typename T>
class ResourceDistribution {
// Class Template for Allocation of Resources (edge and cloud servers) among
// applications.
// Template parameter T is the type of the number of edge and cloud servers allocated
// for each application. In the first stage of our heuristic approach to the platform
// optimization problem we neglect the integrality constraint of the number of edge
// and cloud servers used by the platform, hence allowing T to be float or double.

public:
  // Container that maps each application to its load assigned to edge by the platform
  std::map<std::string,double> edgeLoads={};
  // Container that maps each application to its load assigned to cloud by the platform
  std::map<std::string,double> cloudLoads={};
  // Container that maps each application to the servers allocated by the platform to serve
  // its load in the edge
  std::map<std::string,T> edgeServers={};
  // Container that maps each application to the cloud VMs allocated by the platform to serve
  // its load in the cloud
  std::map<std::string,T> cloudServers={};
};

#endif
