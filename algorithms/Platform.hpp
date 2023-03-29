#ifndef PLATFORM_CLASS_HPP
#define PLATFORM_CLASS_HPP

#include"User.hpp"
#include"vector"
#include<tuple>
#include"set"
#include"array"
#include"functional"
#include<map>
#include"PartialSolution.hpp"
#include<cmath>
#include"Solution.hpp"
#include"Permutations.hpp"
#include"Methods.hpp"
#include"list"
#ifdef GNUPLOT
#include "gnuplot-iostream.hpp"
#endif

class Platform
{
public:

  ///----------------------------------------------------------------------------
  /// PLATFORM DATA STRUCTURES AND PARAMETERS

  std::map<unsigned, User> users={};
  std::map<std::string,std::vector<unsigned>> userMap={};
  std::map<std::string,std::shared_ptr<Application>> apps={};

  // Global ID counter used by the platform to assign IDs when new users are added
  unsigned id_counter = 0;

  // Time Orizon for Resource Allocation (seconds)
  int Time;
  // Max number of edge servers available to the platform
  unsigned Ne;
  // Cost of single edge server
  double ce;
  // Cost of single cloud VM
  double c;
  // Interval boundaries for the extra cost r.
  double rMin;
  double rMax;


  // Constructor
  Platform(unsigned Ne, double ce, double c, double rMin, double rMax, int T);

  //----------------------------------------------------------------------------
  // UPDATE OF DATA STRUCTURES

  // Add single user
  void addUser(const User& user);
  // Add vector of users
  void addUsers(const std::vector<User>& user_vec);
  void print() const;
  void printLong() const;


  //----------------------------------------------------------------------------
  /// UTILITIES
  std::set<double> getChangePoints(double eps) const;

  double getNewR(const std::string& name,
  const std::map<unsigned,bool>& userSecondDepl) const;

  template<typename T>
  double getProfit(const PartialSolution<T>& partSol) const;

  std::shared_ptr<User> getUserByLD(const Solution& sol,
  std::string app, std::string type,long unsigned int idx) const;

  bool checkSolution(const Solution&) const;

  //----------------------------------------------------------------------------
  // THEOREMS - KKT
  ResourceDistribution<double> onlyEdge(const std::map<unsigned,bool>& userSecondDepl,
  const std::map<std::string,double>& appLoad) const;

  ResourceDistribution<double> onlyCloud(const std::map<unsigned,bool>& userSecondDepl,
  const std::map<std::string,double>& appLoad) const;

  ResourceDistribution<double> oneApp(const std::map<unsigned,bool>& userSecondDepl,
  const std::string& name, double load, double maxN) const;

  //----------------------------------------------------------------------------

  // Optimal Price Algorithm
  std::set<PartialSolution<double>> algorithm2(int n,const
  std::vector<std::vector<std::string>>& orders, bool verbose) const;

  void ComputeOptimalVMs(PartialSolution<double>& partSol,
  const std::map<std::string,double>& appLoad) const;

  //----------------------------------------------------------------------------

  // Users' Assignment Algorithm
  Solution algorithm3(const std::set<PartialSolution<double>> & eliteSol,
  bool) const;

  void CountEdgeVMs(Solution & sol, const PartialSolution<double> & partSol,
  const std::string & name, bool) const;

  void MoveUser(Solution & sol, unsigned id, double minLD) const;

  void MoveApplication(Solution & sol, const std::string & name) const;

  void CountCloudVMs(Solution& sol, const std::string& name) const;

  //----------------------------------------------------------------------------
  Solution algorithm4(const Solution&,int,std::default_random_engine&,
  bool,Method,int,bool,int,int) const;

  std::vector<Solution> MoveCloudToEdge(const Solution&,
  const std::string&, bool) const;

  std::vector<Solution> MoveEdgeToCloud(const Solution& sol,
  const std::string& name, int diff_ne, bool) const;

  //----------------------------------------------------------------------------
  Solution heuristicApproach(int n,const std::vector<std::vector<std::string>>&,
  int,bool,std::default_random_engine&,Method,int,bool,int,int) const;
};

#endif
