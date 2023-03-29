#include"Platform.hpp"

Platform::Platform(unsigned Ne, double ce, double c, double rMin, double rMax, int T):
Ne(Ne), ce(ce), c(c), rMin(rMin), rMax(rMax) , Time(T) {}

void Platform::addUser(const User& user){
  users[user.id] = user;
  auto it = userMap.find(user.app->name);
  if(it == userMap.cend()) userMap[user.app->name] = {user.id};
  else it->second.push_back(user.id);
}

void Platform::addUsers(const std::vector<User>& user_vec){
  for(auto u : user_vec)
    addUser(u);
}

//------------------------------------------------------------------------------
// UTILITIES

// Platform Print
void Platform::print() const{
  std::cout << "Platform Info: " << std::endl;
  std::cout << " N° of applications: " << apps.size() << std::endl;
  std::cout << " N° of users: " << users.size() << ", distributed as follows: " << std::endl;
  for(auto it=userMap.cbegin(); it!=userMap.cend(); ++it){
    std::cout << "  " << it->second.size() << " users of " <<
    it->first << std::endl;
  }
}

void Platform::printLong() const{
  using std::cout;
  using std::endl;
  cout << "Cost edge: " << ce << endl;
  cout << "Cost cloud: " << c << endl;
  cout << "Max nodes: " << Ne << endl;
  cout << "rMax: " << rMax << endl;
  cout << "rMin: " << rMin << endl;
  for(auto [name,app] : apps){
    cout << "Users of " << name << ": ";
    for(auto id : userMap.find(name)->second) cout << id << ", ";
    cout << endl;
    cout << "\tD_e: " << app->D_e << endl;
    cout << "\tD_c: " << app->D_c << endl;
    cout << "\tlambda: " << app->lambda << endl;
    cout << "\tR: " << app->R << endl;
    cout << "\tdelta: " << app->delta << endl;
    cout << "\tgamma: " << app->gamma << endl;
    cout << "\tr_1: " << app->r_1 << endl;
    cout << endl;
  }
  cout << "----------------------------------------------------------------------" << endl;
  for(auto [id,u] : users){
    bool s_1 = u.app->lambda*(u.T*u.T)*u.p_1<=u.E && u.app->m_1<=u.M && u.D_1<u.app->R;
    bool s_2 = u.app->lambda*(u.T*u.T)*u.p_2<=u.E && u.app->m_2<=u.M && u.D_2+u.app->delta/u.B+u.app->D_e < u.app->R;

    cout << "User: " << id << endl;
    cout << "\tD_1: " << u.D_1 << endl;
    cout << "\tD_2: " << u.D_2 << endl;
    cout << "\tp_1: " << u.p_1 << endl;
    cout << "\tp_2: " << u.p_2 << endl;
    cout << "\tbeta: " << u.beta << endl;
    cout << "\tB: " << u.B << endl;
    cout << "\talpha: " << u.alpha << endl;
    cout << "\tE: " << u.E << endl;
    cout << "\tM: " << u.M << endl;
    cout << "\tzeta: " << u.zeta << endl;
    cout << "\tT: " << u.T << endl;
    cout << "\tU: " << u.U << endl;
    cout << "\ts_1: " << s_1 << endl;
    cout << "\ts_2: " << s_2 << endl;
    cout << endl;

  }
}

// Select Random Element of vector in input.
int selectRandIdx(const std::vector<int>& idxVec,std::default_random_engine& engine){
  std::uniform_int_distribution<> dist(1,idxVec.size());
  int n = dist(engine);
  return idxVec[n-1];
}

// Selects one neighbor among vector of neighbors, according to method selected
int selectNeighbor(const std::vector<Solution>& neighbors, int method,
const std::vector<int>& idxVec,std::default_random_engine& engine){
  // random
  if(method==RANDOM){
    int randIdx = selectRandIdx(idxVec,engine);
    return randIdx;
  }
  else{
    double bestProfit = 0;
    int bestIdx = 0;
    for(auto idx:idxVec){
      if(neighbors[idx].profit>bestProfit) bestIdx = idx;
    }
    return bestIdx;
  }
}


// Check Solution
// Checks that members of solution in input are consistent (all load constraints
// are satisfied) and that edge servers allocated do not exceed platform avaiability.
// INPUTS: solution to be checked.
// OUTPUTS: false if some error is found, true otherwise.
bool Platform::checkSolution(const Solution& sol) const{
  bool check = true;
  // Check Idle Nodes are >= 0
  unsigned usedVMs = 0;
  for(auto it=sol.edgeServers.cbegin();
  it!=sol.edgeServers.cend();++it) usedVMs += it->second;
  if(usedVMs > Ne){
    std::cerr << " \n !!! ATTENTION: Solution uses more edge servers than what are available. "
    << std::endl;
    return false;
  }
  // Check Users' assignment and edge/cloud vms are coherent for each app:
  for(auto it=userMap.cbegin(); it!=userMap.cend(); ++it){
    // Count how many Users on Edge and on Cloud
    auto app = apps.find(it->first)->second;
    auto Lc = sol.L_c.find(it->first)->second;
    auto Le = sol.L_e.find(it->first)->second;
    auto LambdaC = sol.cloudLoads.find(it->first)->second;
    auto LambdaE = sol.edgeLoads.find(it->first)->second;
    auto n_e = sol.edgeServers.find(it->first)->second;
    auto n_c = sol.cloudServers.find(it->first)->second;
    int nE = 0;
    int nC = 0;
    for(auto id : it->second){
      if(sol.userEdgeAllocation.find(id)->second) ++nE;
      if(sol.userCloudAllocation.find(id)->second) ++nC;
    }
    bool c1 = std::abs(LambdaE-nE*app->lambda)>0.0000001;
    bool c2 = std::abs(Le-LambdaE*app->D_e)>0.0000001;
    bool c3 = std::abs(LambdaC-nC*app->lambda)>0.0000001;
    bool c4 = std::abs(Lc-LambdaC*app->D_c)>0.0000001;
    // std::cout << Le << " " << LambdaE*app->D_e << std::endl;
    // std::cout << Lc << " " << LambdaC*app->D_c << std::endl;
    if(c1 or c2 or c3 or c4){
      std::cout << " \n !!! ATTENTION: Resources allocation is not coherent with Users' assignment"
      << std::endl;
      std::cout << it->first << std::endl;
      std::cout << "LambdaE/lambda: " << LambdaE/app->lambda << ", Le/De: " << Le/app->D_e << std::endl;
      std::cout << "LambdaC/lambda: " << LambdaC/app->lambda << ", Le/De: " << Lc/app->D_c << std::endl;
      std::cout << "Le,LambdaE,n_e,nE: " << Le << " " << LambdaE << " " << n_e << " " << nE << std::endl;
      std::cout << "Lc,LambdaC,n_c,nC: " << Lc << " " << LambdaC << " " << n_c << " " << nC << std::endl;
      std::cout << c1 << " " << c2 << " " << c3 << " " << c4 << std::endl;
      check = false;
    }
  }
  return check;
}

// Get Decision Change Points of Users
// Returns a set of relevant values for the extra cost r, such that one or
// multiple users change their deployment decision.
// INPUTS: eps -> infinitesimal value that is added and subtracted to critical
//                points of user deployment decisions
// OUTPUTS: set of candidates for optimal extra-cost r.
std::set<double> Platform::getChangePoints(double eps=0.000000001) const{
  std::set<double> points={};
  for(auto it=users.cbegin(); it!=users.cend(); ++it){
    //points.insert(it->second.changeDeplPoint()+eps);
    points.insert(it->second.changeDeplPoint()-eps);
    //points.insert(it->second.dropDeplPoint()+eps);
    points.insert(it->second.dropDeplPoint()-eps);
  }
  return points;
}

// Get New Time-Response Constraint for application
// Returns new time-response constraint used in relazed optimization problem for
// selected application.
// INPUTS: name -> name of considered application
//         userSecondDepl -> Map container that maps users' ID to bool value
//                           indicating whether he/she chose 2° depl.
// OUTPUTS: new time response limit.
double Platform::getNewR(const std::string& name,
const std::map<unsigned,bool>& userSecondDepl) const{
  double sum=0;
  unsigned count=0;
  for(auto id : userMap.find(name)->second){
    if(userSecondDepl.find(id)->second){
      sum += users.find(id)->second.getV();
      count += 1;
    }
  }
  auto app = apps.find(name)->second;
  if(count==0)
    return app->R;
  else
    return app->R - sum/count;
}

// Get Profit
// Returns profit of Partial Solution in Input
// INPUTS: partSol -> partial solution containing indormation of resources used by platform
// OUTPUTS: profit
template<typename T>
double Platform::getProfit(const PartialSolution<T>& partSol) const{
  double revenue = 0;
  double cost = 0;
  for(auto it=userMap.cbegin(); it!=userMap.cend(); ++it){
    auto app = apps.find(it->first)->second;
    for(auto id : it->second){
      User u = users.find(id)->second;
      revenue += (app->r_1*partSol.userFirstDepl.find(id)->second +
      (app->r_1+app->gamma*partSol.r)*partSol.userSecondDepl.find(id)->second)*u.T/3600;
    }
    cost += (ce*partSol.edgeServers.find(app->name)->second +
    c*partSol.cloudServers.find(app->name)->second)*this->Time/3600;
  }
  return revenue-cost;
}

// Get User By Local Delay
// Ranks users of selected application assigned to selected type by local delay,
// and returns a pointer to the user in position idx.
// INPUTS: sol -> solution containing information about users' assignment
//         app -> name of application selected
//         type -> must be "edge" or "cloud"
//         idx -> user ranking by local delay
// OUTPUTS: pointer to user in position idx
std::shared_ptr<User> Platform::getUserByLD(const Solution& sol,
std::string app, std::string type, long unsigned int idx) const{
  std::multimap<double,unsigned> usersByLD={};
  for(auto id : userMap.find(app)->second){
    if((sol.userCloudAllocation.find(id)->second && type == "cloud") or
    (sol.userEdgeAllocation.find(id)->second && type == "edge"))
      usersByLD.insert({users.find(id)->second.getLD(),id});
  }
  if(usersByLD.size()-1<idx){
    std::cout << "  !!! There are no users in " << type
    << " for index " << idx << std::endl;
    std::exit(1);
  }
  std::vector<unsigned> vecID = {};
  // Copy IDs in a vector ordered by LD
  for(auto it=usersByLD.cbegin(); it!=usersByLD.cend(); ++it)
    vecID.push_back(it->second);
  auto ID = vecID[idx];
  auto res = std::make_shared<User>(users.find(ID)->second);
  // std::cout << "User ID: " << res->id << " , LD: " << res->getLD() << std::endl;
  return res;
}

//-----------------------------------------------------------------------------
// THEOREMS - KKT
// Only Edge Scenario
// Returns optimal resource allocation under assumption of only edge scenario
// INPUTS: usedSecondDepl -> map container where the value is true if the user
//                           corresponding to the ID of the key selected 2° deployment
//         appLoad -> map container that maps applications to the total load of users
//                    who chose the 2° depl.
// OUTPUTS: optimal resource allocation
ResourceDistribution<double> Platform::onlyEdge(const std::map<unsigned,bool>&
userSecondDepl, const std::map<std::string,double>& appLoad) const {
  ResourceDistribution<double> res;
  for(auto it=appLoad.cbegin();it!=appLoad.cend();++it){
    auto app = apps.find(it->first)->second;
    double barR = this->getNewR(it->first, userSecondDepl);
    std::string appName = it->first;
    double load = it->second;
    res.edgeLoads[appName] = load;
    res.cloudLoads[appName] = 0;
    res.edgeServers[appName] = barR*app->D_e*load/(barR-app->D_e);
    res.cloudServers[appName] = 0;
  }
  return res;
}

// Only Cloud Scenario
// Returns optimal resource allocation under assumption of only cloud scenario
// INPUTS: usedSecondDepl -> map container where the value is true if the user
//                           corresponding to the ID of the key selected 2° deployment
//         appLoad -> map container that maps applications to the total load of users
//                    who chose the 2° depl.
// OUTPUTS: optimal resource allocation
ResourceDistribution<double> Platform::onlyCloud(const std::map<unsigned,bool>&
userSecondDepl, const std::map<std::string,double>& appLoad) const {
  ResourceDistribution<double> res;
  for(auto it=appLoad.cbegin();it!=appLoad.cend();++it){
    auto app = apps.find(it->first)->second;
    double barR = this->getNewR(it->first, userSecondDepl);
    std::string name = it->first;
    double load = it->second;
    res.edgeLoads[name] = 0;
    res.cloudLoads[name] = load;
    res.edgeServers[name] = 0;
    res.cloudServers[name] = barR*app->D_c*load/(barR-app->D_c);
  }
  return res;
}

// One Application - Mixed Scenario
// Returns optimal resource allocation under assumption of one application.
// INPUTS: usedSecondDepl -> map container where the value is true if the user
//                           corresponding to the ID of the key selected 2° deployment
//         name -> name of application to be considered in this scenario
//         load -> total load of users of input-application who chose 2° depl.
//         maxN -> max number of edge servers available to platform.
// OUTPUTS: optimal resource allocation
ResourceDistribution<double> Platform::oneApp(const std::map<unsigned,bool>&
userSecondDepl, const std::string& name, double load, double maxN) const{
  ResourceDistribution<double> res;
  auto app = apps.find(name)->second;
  double barR = this->getNewR(name, userSecondDepl);
  if(load <= maxN*(barR-app->D_e)/(barR*app->D_e)){
    res.edgeLoads[name] = load;
    res.cloudLoads[name] = 0;
    res.edgeServers[name] = barR*app->D_e*load/(barR-app->D_e);
    res.cloudServers[name] = 0;
  }
  else{
    double lambda_e = maxN*load*(barR-std::sqrt(app->D_c*app->D_e))/
    (maxN*app->D_e + barR*load*app->D_e - maxN*std::sqrt(app->D_e*app->D_c));
    double n_c = (app->D_c*load*(barR*app->D_e*load-maxN*(barR-app->D_e)))/
    (maxN*std::pow(std::sqrt(app->D_e)-std::sqrt(app->D_c),2) + app->D_e*load*
    (barR-app->D_c));
    res.edgeLoads[name] = lambda_e;
    res.cloudLoads[name] = load-lambda_e;
    res.edgeServers[name] = maxN;
    res.cloudServers[name] = n_c;
  }
  return res;
}

//------------------------------------------------------------------------------
// ALGORITHM 2

// Compute Optimal VMs
// Computes optimal number of cloud and edge servers according to input-app load
// and updates partSol members.
// INPUTS: - partial Solution to be updated
//         - Map container mapping each application by their name to their total load
// OUTPUTS: -
void Platform::ComputeOptimalVMs(PartialSolution<double>& partSol,
const std::map<std::string,double>& appLoad) const{
  partSol.updateResources(onlyEdge(partSol.userSecondDepl, appLoad));
  double usedVMs = 0;
  for(auto it=partSol.edgeServers.cbegin();
  it!=partSol.edgeServers.cend(); ++it) usedVMs += it->second;
  //std::cout << " Used Edge Nodes: " << usedVMs << std::endl;
  if(usedVMs > Ne){
    ResourceDistribution<double> resCloud = onlyCloud(partSol.userSecondDepl, appLoad);
    std::size_t idx = partSol.order.size();
    while(usedVMs>Ne){
      idx -= 1;
      auto app = apps.find(partSol.order[idx])->second;
      //std::cout << "Moving app: " << app->name << std::endl;
      usedVMs -= partSol.edgeServers[app->name];
      //std::cout << " Updated used VMS: " << usedVMs << std::endl;
      partSol.edgeServers[app->name] = 0;
      partSol.edgeLoads[app->name] = 0;
      partSol.cloudServers[app->name] = resCloud.cloudServers.find(app->name)->second;
      partSol.cloudLoads[app->name] = resCloud.cloudLoads.find(app->name)->second;
    }
    double leftVMs = Ne-usedVMs;
    auto appName = partSol.order[idx];
    ResourceDistribution<double> solMix = oneApp(partSol.userSecondDepl, appName,
    appLoad.find(appName)->second,leftVMs);
    partSol.edgeServers[appName] = solMix.edgeServers[appName];
    partSol.edgeLoads[appName] = solMix.edgeLoads[appName];
    partSol.cloudServers[appName] = solMix.cloudServers[appName];
    partSol.cloudLoads[appName] = solMix.cloudLoads[appName];
  }
}

// ALGORITHM 2 - Optimal Prices Algorithm
// Returns set of Elite Solutions: partial solutions with largest profit of relaxed
// platform optimization problem.
// INPUTS: - number of Elite Solutions to be returned
//         - method to be used, must be "tabu" or "permutations"
//         - verbosity option
// OUTPUTS set of Partial Solutions with largest profit.
std::set<PartialSolution<double>> Platform::algorithm2(int n,
const std::vector<std::vector<std::string>>& orders, bool verbose) const{
  if(verbose){
    std::cout << "________________________________________________________________" << std::endl;
    std::cout << "OPTIMAL PRICES ALGORITHM" << std::endl;
  }
  std::set<PartialSolution<double>> solutions;
  std::set<double> changePoints = getChangePoints();
  if(verbose) std::cout << "Number of r-values to inspect: " << changePoints.size() << "\n" << std::endl;
  for(double r : changePoints){
    if(r<=rMax && r>=rMin){
      for(auto order:orders){
        PartialSolution<double> partSol(r);
        partSol.order = order;
        std::map<std::string,double> appLoad;
        for(auto it=userMap.cbegin(); it!=userMap.cend(); ++it){ // Compute appLoad
          auto app = apps.find(it->first)->second;
          appLoad[app->name] = 0;
          //partSol.order.push_back(app->name);
          for(auto id : it->second){
            std::array<bool,2> userDecision = users.find(id)->second.algorithm1(r);
            partSol.userFirstDepl[id] = userDecision[0];
            partSol.userSecondDepl[id] = userDecision[1];
            appLoad[app->name] += userDecision[1]*app->lambda;
          }
        }
        ComputeOptimalVMs(partSol, appLoad);
        partSol.profit = getProfit(partSol);
        solutions.insert(partSol);
      }
    }
  }
  int count = 0;
  std::set<PartialSolution<double>> res={};
  if(verbose) std::cout << "Take only the " << n << " partial solutions with best profit: \n" << std::endl;
  // Filter out useless Solutions (same users' decisions, but lower r)
  for(auto it=solutions.crbegin(); it!=solutions.crend() && count<n; ++it){
    bool useless_sol = false;
    for(auto it1=res.cbegin(); it1!=res.cend(); ++it1){
      if(it1->userFirstDepl==it->userFirstDepl or
        it1->userSecondDepl==it->userSecondDepl) useless_sol = true;
    }
    if(not useless_sol){
      if(verbose) it->printLess();
      res.insert(*it);
      ++count;
    }
  }
  return res;
}

//------------------------------------------------------------------------------
// ALGORITHM 3

// ALGORITHM 3 - Users' Assignment Algorithm
// Computes optimal users' assignment of partial solutions in input and returns
// best solution found.
// INPUTS: - set of elite partial solutions
//         - verbosity option
// OUTPUTS: solution with largest profit found.
Solution Platform::algorithm3(const std::set<PartialSolution<double>>& eliteSol,
bool verbose) const{

  if(verbose){
    std::cout << "________________________________________________________________" << std::endl;
    std::cout << "USERS' RESOURCE ASSIGNMENT " << std::endl;
  }
  std::set<Solution> solutions={};
  for(auto it=eliteSol.cbegin(); it!=eliteSol.cend(); ++it){
    if(verbose){
      std::cout << "\n## Current Partial Solution ##" << std::endl;
      it->print();
      std::cout << "Now computing optimal users' resource assigment:" << std::endl;
    }
    unsigned int barNe = 0;
    Solution sol(*it); // Solution associated to partial solution pointed by current it.
    for(auto it1=apps.cbegin(); it1!=apps.cend(); ++it1){
      auto app=it1->second; // app is pointer to application
      sol.initializeAppValues(app); // Set n. of servers to 0
      CountEdgeVMs(sol, *it, app->name, verbose);
      barNe += sol.edgeServers.find(app->name)->second;
    }
    if(barNe > Ne){
      unsigned int usedVMs = 0;
      std::size_t idx = 0;
      while(usedVMs <= Ne){
        barNe = usedVMs;
        std::string appName = it->order[idx];
        usedVMs += sol.edgeServers.find(appName)->second;
        ++idx;
      }
      if(verbose){
        std::cout << "Exceeding number of Edge Servers available:" << std::endl;
        if(barNe<Ne) std::cout << "  --> moving users of "<< it->order[idx-1] << " to cloud.";
      }
      int leftVMs = Ne-barNe;
      if(leftVMs==0) MoveApplication(sol,it->order[idx-1]);
      else{
        auto app = apps.find(it->order[idx-1])->second; // Pointer to application whose user need to me moved in c
        int userMoved = 0;
        while(sol.edgeServers.find(app->name)->second > leftVMs){
          std::multimap<double,unsigned> userByViolation={};
          for(auto id : userMap.find(app->name)->second){
            if(sol.userEdgeAllocation.find(id)->second){
              double V = users.find(id)->second.getV();
              userByViolation.insert({V,id});
            }
          }
          unsigned idLowestV = userByViolation.cbegin()->second;
          double minLD = users.find(userByViolation.crbegin()->second)->second.getLD();
          MoveUser(sol, idLowestV, minLD);
          ++userMoved;
        }
        if(verbose){
          std::cout << " Users moved: " << userMoved << std::endl;
        }
      }
      for(auto i=idx; i<apps.size(); ++i){
        if(it->edgeLoads.find(it->order[i])->second > 0 or
        it->cloudLoads.find(it->order[i])->second > 0){
          if(verbose) std::cout << "  --> Assigning " << it->order[i] << " entirely on cloud." << std::endl;
          MoveApplication(sol, it->order[i]);
        }
      }
    }
    for(auto name : it->order)
      CountCloudVMs(sol, name);
    sol.profit = getProfit(sol);
    bool check = checkSolution(sol);
    if(!check) std::exit(1);
    solutions.insert(sol);
    if(verbose){
      std::cout << "\n## Final Resource Assignment for Current Solution ##" << std::endl;
      sol.print();
      std::cout << "\t" << std::endl;
    }
  }
  return *solutions.rbegin();
}

// Count Edge VMs
// Updates edge servers assigned to input-application of input-solution according to resource allocation of
// input-partial solution, such that all load and time-response constraints are satisfied.
// INPUTS: - solution to be updated
//         - partial solution with approximated resource allocation
//         - name of application to be considered
//         - verbosity option
// OUTPUTS: -
void Platform::CountEdgeVMs(Solution& sol, const PartialSolution<double>& partSol,
const std::string & name, bool verbose) const{
  // Initialization
  auto app = apps.find(name)->second; // app is a pointer to application with name
  std::multimap<double, unsigned> violations={}; // Order users by Violations
  unsigned countEdge=0;
  unsigned countCloud=0; // n. of users that need to be served on edge or cloud
  if(partSol.cloudServers.find(name)->second==0 &&
  partSol.edgeServers.find(name)->second > 0){ // No Cloud Servers Case
    for(auto id : userMap.find(name)->second){
      if(partSol.userSecondDepl.find(id)->second){
        sol.userEdgeAllocation[id] = 1;
        sol.userCloudAllocation[id] = 0;
        countEdge += 1;
      }
    }
    sol.L_e[name] = app->D_e*app->lambda*countEdge;
    int n_e = std::ceil(partSol.edgeServers.find(name)->second);
    for(auto id : userMap.find(name)->second){
      if(partSol.userSecondDepl.find(id)->second){
        User u = users.find(id)->second;
        double V = u.D_2+app->delta/u.B + app->D_e/
        (1-sol.L_e.find(name)->second/n_e)-app->R;
        // violations[V] = *it; operator [] not defined for multimaps
        violations.insert({V,id});
      }
    }
    double maxV = violations.rbegin()->first;
    if(maxV <= 0){
      sol.edgeServers[name] = n_e;
    }
    else{
      User uMaxV = users.find(violations.rbegin()->second)->second;
      double LD = uMaxV.getLD();
      sol.edgeServers[name] =
      std::ceil(sol.L_e.find(name)->second*LD/(LD-app->D_e));
    }
    sol.edgeLoads[name] = partSol.edgeLoads.find(name)->second;
  }
  if(partSol.cloudServers.find(name)->second > 0 &&
  partSol.edgeServers.find(name)->second > 0){ // Case with cloud Servers
    for(auto id : userMap.find(name)->second){
      if(partSol.userSecondDepl.find(id)->second){
        User u = users.find(id)->second;
        double V = u.getV();
        //violations[V] = *it;
        violations.insert({V,id});
      }
    }
    std::vector<unsigned> v={}; // vector of users' id ordered by V
    for(auto it=violations.crbegin(); it!=violations.crend(); ++it){
      v.push_back(it->second);
    }
    std::size_t idx = partSol.cloudLoads.find(name)->second/app->lambda;
    for(std::size_t i=0; i<idx; ++i){
      countCloud += 1;
      sol.userCloudAllocation[v[i]] = 1;
      sol.userEdgeAllocation[v[i]] = 0;
    }
    for(std::size_t i=idx; i<v.size(); ++i){
      countEdge += 1;
      sol.userCloudAllocation[v[i]] = 0;
      sol.userEdgeAllocation[v[i]] = 1;
    }
    User uMaxV = users.find(v[idx])->second;
    sol.L_c[name] = app->D_c*app->lambda*countCloud;
    sol.L_e[name] = app->D_e*app->lambda*countEdge;
    sol.edgeLoads[name] = app->lambda*countEdge;
    sol.cloudLoads[name] = app->lambda*countCloud;
    double LD = uMaxV.getLD();
    sol.edgeServers[name] =
    std::ceil(sol.L_e.find(name)->second*LD/(LD-app->D_e));
  }
  else if(partSol.cloudServers.find(name)->second > 0 &&
  partSol.edgeServers.find(name)->second==0){
    for(auto id : userMap.find(name)->second){
      if(partSol.userSecondDepl.find(id)->second){
        sol.userEdgeAllocation[id] = 0;
        sol.userCloudAllocation[id] = 1;
        countCloud += 1;
      }
    }
    sol.L_c[name] = app->D_c*app->lambda*countCloud;
  }
  // Else: do Nothing.
  if(verbose){
    std::cout << "  " << name << " requires " <<
    sol.edgeServers.find(name)->second << " edge servers." << std::endl;
  }
}

// Move User
// Moves input-user from edge to cloud, updating all the relevant values.
// INPUTS: - solution to be updated
//         - ID of users to be moved
//         - lowest local delay of existig users in the edge.
// OUTPUTS: -
void Platform::MoveUser(Solution& sol, unsigned id, double minLD) const{
  std::string appName = users.find(id)->second.app->name;
  auto app = apps.find(appName)->second; // ptr to Application
  sol.userEdgeAllocation[id] = 0;
  sol.userCloudAllocation[id] = 1;
  sol.edgeLoads[appName] -= app->lambda;
  sol.cloudLoads[appName] += app->lambda;
  sol.L_e[appName] = app->D_e*sol.edgeLoads.find(appName)->second;
  sol.L_c[appName] = app->D_c*sol.cloudLoads.find(appName)->second;
  sol.edgeServers[appName] =
  std::ceil(app->D_e*sol.edgeLoads.find(appName)->second*minLD/(minLD-app->D_e));
}

// Move Application
// Moves all users of input-application in the cloud and updates relevant values.
// INPUTS: - solution to be updated
//         - name of application to be moved
// OUTPUTS: -
void Platform::MoveApplication(Solution& sol, const std::string& name) const{
  auto app = apps.find(name)->second;
  for(auto id : userMap.find(name)->second){
    if(sol.userSecondDepl.find(id)->second){
      sol.userEdgeAllocation[id]=0;
      sol.userCloudAllocation[id]=1;
    }
  }
  sol.L_c[name] += app->D_c*sol.edgeLoads.find(name)->second;
  sol.cloudLoads[name] += sol.edgeLoads.find(name)->second;
  sol.edgeLoads[name]=0;
  sol.edgeServers[name]=0;
  sol.L_e[name]=0;
}

// Cout Cloud VMS
// Updates cloud VMs assigned to input-application of input-solution so that
// all load and time-response constraints are satisfied.
// INPUTS: - solution to be updated
//         - name of application to be considered
// OUTPUTS: -
void Platform::CountCloudVMs(Solution& sol, const std::string& name) const{
  auto app = apps.find(name)->second;
  if(sol.cloudLoads.find(name)->second<0.0000001)
    sol.cloudServers[name]=0;
  else{
    double LD = getUserByLD(sol,name,"cloud",0)->getLD();
    double Lc = sol.L_c.find(name)->second;
    sol.cloudServers[name]=std::ceil(Lc*LD/(LD-app->D_c));
  }
}

//------------------------------------------------------------------------------
// ALGORITHM 4

// Move Cloud To Edge
// Generates vector of neighbors of initial solution by forcing move of users
// from cloud to edge.
// INPUTS: - initial solution
//         - name of application whose users are moved
//         - verbosity option
// OUTPUTS: vector of neighbor solutions
std::vector<Solution> Platform::MoveCloudToEdge(const Solution& sol,
const std::string& name, bool verbose) const{
  if(verbose){
    std::cout << "\n \t -- Move Cloud To Edge -- " << std::endl;
    std::cout << " We reduce by 1 the cloud resources used by " << name
    << std::endl;
    std::cout << " Current n_e: " << sol.edgeServers.find(name)->second
    << ", n_c: " << sol.cloudServers.find(name)->second << std::endl;
  }
  auto app = apps.find(name)->second;
  // Count Users currently served on Cloud
  int nUserCloud = 0;
  for(auto id : userMap.find(name)->second){
    if(sol.userCloudAllocation.find(id)->second) ++nUserCloud;
  }
  int n_c_prime = sol.cloudServers.find(name)->second-1;
  auto userPtr = getUserByLD(sol,name,"cloud",0);
  auto minLDc = userPtr->getLD();
  int maxUsersCloud = std::floor((minLDc-app->D_c)*n_c_prime)
  /(app->lambda*minLDc*app->D_c);
  int usersToMove = nUserCloud-maxUsersCloud;
  if(verbose) std::cout << " --> Need to move " << usersToMove << " users in the edge" << std::endl;
  // Move UsersToMove users with highest LD (= faster devices):
  Solution neighbor = sol;
  // if(neighbor.edgeLoads.find(name)->second>0){
  //   std::cout << " Initial minLD in edge: " <<
  //   getUserByLD(neighbor,name,"edge",0)->getLD() << std::endl;
  // }
  for(int i=0; i<usersToMove; ++i){
    int idx = nUserCloud-1-i;
    auto ID = getUserByLD(sol,name,"cloud",idx)->id;
    neighbor.userCloudAllocation[ID] = false;
    neighbor.userEdgeAllocation[ID] = true;
    neighbor.L_c[name] -= app->lambda*app->D_c;
    neighbor.L_e[name] += app->lambda*app->D_e;
    neighbor.cloudLoads[name] -= app->lambda;
    neighbor.edgeLoads[name] += app->lambda;
  }
  auto minLDe = getUserByLD(neighbor,name,"edge",0)->getLD();
  // std::cout << " Post minLD in edge "  << minLDe << std::endl;
  int n_e_prime =
  std::ceil(neighbor.L_e.find(name)->second*minLDe/(minLDe-app->D_e));
  int diff_ne = n_e_prime - sol.edgeServers.find(name)->second;
  // Compute IDLE nodes in edge
  int usedVMs = 0;
  for(auto it = sol.edgeServers.cbegin(); it!= sol.edgeServers.cend(); ++it)
    usedVMs += it->second;
  int idle = Ne-usedVMs;
  diff_ne -= idle;
  if(usedVMs-sol.edgeServers.find(name)->second<diff_ne) return {};
  if(diff_ne<=0){
    if(verbose) std::cout << " Idle edge nodes are enough to cover extra resources needed." << std::endl;
    neighbor.edgeServers[name] = n_e_prime;
    CountCloudVMs(neighbor,name);
    //neighbor.cloudServers[name] = n_c_prime;
    neighbor.profit = getProfit(neighbor);
    return {neighbor};
  }
  else{
    if(verbose){
      std::cout << " To serve these users, " << diff_ne+idle <<
      " additional edge servers are needed." << std::endl;
      std::cout << " --> Edge servers must be freed from other applications."
      << std::endl;
    }
    CountCloudVMs(neighbor,name);
    //neighbor.cloudServers[name] = n_c_prime*(maxUsersCloud>0);
    neighbor.edgeServers[name] = n_e_prime;
    return MoveEdgeToCloud(neighbor, name, diff_ne, verbose);
  }
}

// Move Edge to Cloud
// Generate vector of neighbors of initial solution by forcing move of users
// from edge to cloud.
// INPUTS: - initial solution
//         - name of application whose users are moved
//         - number of edge servers that need to be freed from other applications
//         - verbosity option
// OUTPUTS: vector of neighbor solutions
std::vector<Solution> Platform::MoveEdgeToCloud(const Solution& sol,
const std::string& name, int diff_ne, bool verbose) const{
  if(verbose) std::cout << "\n \t -- Move Edge To Cloud -- " << std::endl;
  std::vector<Solution> res = {};
  // Loop on applications with Edge Load
  for(auto it=apps.cbegin(); it!=apps.cend(); ++it){
    auto cName = it->first;
    auto cApp = it->second;
    auto curr_n_e = sol.edgeServers.find(cName)->second;
    if(curr_n_e > 0 && cName != name){
      Solution solCopy = sol;
      if(verbose){
        std::cout << " # Now Considering " << cName << "\n  Resources used: "
        << "n_e = "<< curr_n_e << ", n_c = " << sol.cloudServers.find(cName)->second << std::endl;
      }
      int nUserEdge = 0; // Count Users currently served on edge of current app
      for(auto id : userMap.find(cName)->second){
        if(sol.userEdgeAllocation.find(id)->second) ++nUserEdge;
      }
      // std::cout << "  Users currently using edge: " << nUserEdge << std::endl;
      if(curr_n_e >= diff_ne){
        double minLDe = getUserByLD(solCopy,cName,"edge",0)->getLD();
        int barNe = std::floor((minLDe-cApp->D_e)*(curr_n_e-diff_ne)
        /(cApp->lambda*minLDe*cApp->D_e));
        int usersToMove=nUserEdge-barNe; // Users to move on Cloud
        if(verbose){
          std::cout << "  Current Users served on Edge: " << nUserEdge <<
          ", Users to be Moved: " << usersToMove << std::endl;
        }
        if(usersToMove <= 0) solCopy.edgeServers[cName] = curr_n_e-diff_ne;
        else if(barNe==0){
          MoveApplication(solCopy,cName);
          CountCloudVMs(solCopy,cName);
        }
        else{
          for(int i=0; i<usersToMove; ++i){
            int idx = nUserEdge-1-i;
            auto ID = getUserByLD(solCopy,cName,"edge",idx)->id;
            MoveUser(solCopy,ID,minLDe);
            CountCloudVMs(solCopy,cName);
          }
        }
        if(verbose){
          std::cout << "  New n_e: " << solCopy.edgeServers.find(cName)->second;
          std::cout << "  New n_c: " << solCopy.cloudServers.find(cName)->second << std::endl;
        }
        solCopy.profit = getProfit(solCopy);
        res.push_back(solCopy);
      }
      else{
        // Move all users running cName to cloud:
        if(verbose) std::cout << "  Moving entire " << cName << " to cloud." << std::endl;
        MoveApplication(solCopy,cName);
        CountCloudVMs(solCopy,cName);
        if(verbose){
          std::cout << "  New n_e: " << solCopy.edgeServers.find(cName)->second;
          std::cout << "  New n_c: " << solCopy.cloudServers.find(cName)->second << std::endl;
        }
        // Select randomly other application
        int diff_ne1 = diff_ne - curr_n_e;
        for(auto it2=apps.cbegin(); it2!=apps.cend(); ++it2){
          auto cName2 = it2->first;
          auto cApp2 = it2->second;
          auto curr_n_e2 = sol.edgeServers.find(cName2)->second;
          if(curr_n_e2 >= diff_ne1 && cName2!=cName && cName2!=name){
            if(verbose){
              std::cout << "  # Considering additional application: " << cName2 << std::endl;
              std::cout <<  "  Resources used: " << "n_e = "<< curr_n_e2 <<
              ", n_c = " << sol.cloudServers.find(cName2)->second << std::endl;
            }
            int nUserEdge2 = 0; // Count Users currently served on edge of current app
            for(auto id : userMap.find(cName2)->second){
              if(sol.userEdgeAllocation.find(id)->second) ++nUserEdge2;
            }
            double minLDe2 = getUserByLD(solCopy,cName2,"edge",0)->getLD();
            int barNe2 = std::floor((minLDe2-cApp2->D_e)*(curr_n_e2-diff_ne1)
            /(cApp2->lambda*minLDe2*cApp2->D_e));
            int usersToMove2=nUserEdge2-barNe2; // Users to move on Cloud
            if(usersToMove2<=0) solCopy.edgeServers[cName2] = curr_n_e2-diff_ne1;
            else if(barNe2==0){
              MoveApplication(solCopy,cName2);
              CountCloudVMs(solCopy,cName2);
            }
            else{
              //double minLDe2 = getUserByLD(solCopy,cName2,"edge",0)->getLD();
              for(int i=0; i<usersToMove2; ++i){
                int idx = nUserEdge2-1-i;
                auto ID = getUserByLD(solCopy,cName2,"edge",idx)->id;
                MoveUser(solCopy,ID,minLDe2);
                CountCloudVMs(solCopy,cName2);
              }
            }
            if(verbose){
              std::cout << "  New n_e: " << solCopy.edgeServers.find(cName2)->second;
              std::cout << " New n_c: " << solCopy.cloudServers.find(cName2)->second << std::endl;
            }
            solCopy.profit = getProfit(solCopy);
            res.push_back(solCopy);
          }
        }
      }
    }
  }
  return res;
}


// Plots the Profit Prfiles of Current and Best Solution
void profitPlot(const std::vector<int>& iterVec, const std::vector<double>&
currP, const std::vector<double>& bestP, const std::string& title){
  #ifdef GNUPLOT
  Gnuplot gp;
  gp << "TITLE = " << title << std::endl;
  gp << "set title TITLE" << std::endl;
  gp << "set key right bottom" << std::endl;
  gp << "set xlabel 'n. of iterations'" << std::endl;
  gp << "set ylabel 'Platform profit'" << std::endl;
  gp << "plot" << gp.file1d(std::tie(iterVec,currP)) << "w lp ps 0.5 lc rgb 'blue' lw 0.5 title 'Current Solution Profit profile',"
  << gp.file1d(std::tie(iterVec,bestP)) << "w lp ps 0.5 lc rgb 'red' lw 0.5 title 'Best Solution Profit profile'" << std::endl;
  #endif
}



// ALGORITHM 4 - Tabu Search
// Implements Tabu Search attempting to improve the profit of initial solution
// INPUTS: - sol: Initial solution
//         - maxIter: max number of iterations
//         - engine: engine for random selection
//         - verbose: verbosity option
//         - method: method for random selection of neighbors
//         - rank: rank (for plotting profit profile of tabu search).
// OUTPUTS: best solution found (largest profit)
Solution Platform::algorithm4(const Solution& sol, int maxIter,
std::default_random_engine& engine, bool verbose,Method method,
int rank,bool plot, int tabuAttribute,int tabuSize) const{

  if(verbose){
  std::cout << "________________________________________________________________" << std::endl;
  std::cout << "TABU SEARCH \n" << std::endl;
  }
  // Title for plot.
  std::string title;
  switch(method){
    case RANDOM: title = "'Tabu Search - RANDOM neighbor selection"; break;
    case BEST_PROFIT: title = "'Tabu Search - BEST PROFIT neighbor selection'"; break;
  }

  int count = 0; // Counter for how many times Algorithm improves Our Solution.
  Solution bestSol = sol;
  Solution currSol = sol;
  // Initialize vectors for current profit and best-current profit
  std::vector<double> currP={};
  std::vector<double> bestP={};
  std::vector<int> iterVec={};
  // Defining Tabu List with Custom Comparison Operator:
  // auto cmp = [](const Solution& lhs, const Solution& rhs){return comparisonOperator(lhs,rhs);};
  // std::set<Solution, decltype(cmp)> tabuList = {};
  std::list<TabuAttribute> tabuList={};
  for(int i=0; i<maxIter; ++i){
    if(verbose) std::cout << " \n Current Iter : " << i << std::endl;
    // updates
    currP.push_back(currSol.profit);
    bestP.push_back(bestSol.profit);
    iterVec.push_back(i);
    // Find App with min cloud usage:
    std::map<double,std::string> VMuse = {};
    for(auto [name,val] : currSol.L_c){
      if(val > 0.00000001) {
        double ratio = val/currSol.cloudServers.find(name)->second;
        VMuse[ratio] = name;
      }
    }
    if(VMuse.empty()){
      if(verbose){
        std::cout << " **** There are no applications using cloud resources - Tabu Search not able to improve current Solution. **** EXIT 1 \n       Returning Best Solution found "
        << std::endl;
      }
      if(plot) profitPlot(iterVec,currP,bestP,title);
      return bestSol;
    }

    // std::string name = VMuse.cbegin()->second;
    // std::vector<Solution> neighbors = MoveCloudToEdge(currSol,name,verbose);

    std::vector<Solution> neighbors = {};
    for(auto [ratio,name] : VMuse){
      auto curr = MoveCloudToEdge(currSol,name,verbose);
      for(auto neigh : curr) neighbors.push_back(neigh);
    }

    if(verbose) std::cout << "\n  " << neighbors.size() << " neighbors found" << std::endl;
    while(neighbors.empty()){ // Check for neighbors
      auto iter = ++VMuse.cbegin();
      if(iter == VMuse.cend()){
        if(verbose){
          std::cout << " **** No neighbors found - Tabu Search not able to improve current Solution. **** EXIT 2 \n       Returning Best Solution found"
          << std::endl;
          std::cout << " Tabu Search stopped after " << i << " iterations." << std::endl;
          std::cout << " Overall Tabu Serch improved Initial Solution "
          << count << " times. " << std::endl;
        }
        if(plot) profitPlot(iterVec,currP,bestP,title);
        return bestSol;
      }
      // If no neighbors are possible for app with lowest cloud usage try with second
      neighbors = MoveCloudToEdge(currSol,iter->second,verbose);
    }

    std::vector<int> idxVec(neighbors.size(),0);
    for(int i=0;i<neighbors.size();++i) idxVec[i]=i;
    auto idxNeigh = selectNeighbor(neighbors,method,idxVec,engine);
    auto neighbor = neighbors[idxNeigh];
    bool check = checkSolution(neighbor);
    if(!check) std::exit(1);
    TabuAttribute neighborAttribute(neighbor);

    while(true){
      // if(tabuList.find(neighbor)!=tabuList.end()){
      bool tabuListed = false;
      for(const auto & el : tabuList){
        if(el==neighborAttribute or (el.compOp(neighborAttribute) && tabuAttribute==0)){
          tabuListed = true;
          break;
        }
      }

      if(tabuListed){
        auto newEnd = std::remove(idxVec.begin(), idxVec.end(),idxNeigh);
        idxVec.erase(newEnd,idxVec.end());
        if(idxVec.empty()) {
          if(verbose){
            std::cout << "\n **** Tabu Search is stuck with no more neighbors to visit. **** EXIT 3 \n      Returning Best Solution found. "
            << std::endl;
            std::cout << " Overall Tabu Serch improved Initial Solution "
            << count << " times. " << std::endl;
          }
          if(plot) profitPlot(iterVec,currP,bestP,title);
          return bestSol;
        }
        else{
          idxNeigh = selectNeighbor(neighbors,method,idxVec,engine);
          neighbor = neighbors[idxNeigh];
          neighborAttribute = TabuAttribute(neighbor);
        }
      }
      else{
        // tabuList.insert(neighbor);
        tabuList.push_back(neighborAttribute);
        if(tabuList.size() > tabuSize){
          tabuList.pop_front();
          if(verbose) std::cout << "  Reached max size of Tabu List, removing first element of the list" << std::endl;
        }
        currSol = neighbor;
        if(currSol.profit > bestSol.profit) {
          ++count;
          if(verbose) std::cout << " ** New Best Solution found ** " << std::endl;
          bestSol=currSol;
        }
        break;
      }
    }
  }
  if(verbose){
    std::cout << "\n **** Tabu Search has reached max iteration **** EXIT 4 \n      Returning Best Solution found. "
    << std::endl;
    std::cout << " Overall Tabu Serch improved Initial Solution "
    << count << " times. " << std::endl;
  }
  if(plot) profitPlot(iterVec,currP,bestP,title);
  return bestSol;
}


// HEURISTIC APPROACH
// Solves the Complete Platform Optimization problem by running
// all previous algorithms
// INPUTS:  - n: size of elite solutions
//          - orders: vector of permutations to be considered
//          - maxIter: max n. of iterations for Tabu Search
//          - verbose: verbosity option
//          - engine: engine for Tabu Search
//          - rank: number of rank for plot
//          - plot: plot option
// OUTPUTS: Final Solution found by heuristic Approach
Solution Platform::heuristicApproach(int n,const std::vector<std::vector<std::string>>& orders,
int maxIter,bool verbose,std::default_random_engine& engine,Method method,
int rank,bool plot,int tabuAttribute,int tabuSize) const{
  if(orders.size()==0) return Solution();
  auto eliteSol = algorithm2(n,orders,verbose); // Optimal Price Algorithm
  auto sol = algorithm3(eliteSol,verbose); // Users' Resource Assignment Algorithm
  return algorithm4(sol,maxIter,engine,verbose,method,rank,plot,tabuAttribute,tabuSize); // Tabu Search
}
