#ifndef CONFIG_PARAMS_HPP
#define CONFIG_PARAMS_HPP

#include<random>
#include"Application.hpp"
#include<memory>
#include"User.hpp"
#include"GetPot"
#include"Platform.hpp"
#include"json.hpp"
#include<fstream>

// Creates Applications
inline
std::vector<std::shared_ptr<Application>> create_apps(int n,double lR,double uR,double lL,
double uL,double lr1,double ur1,std::default_random_engine& engine){
  // std::cout << "Creating Apps" << std::endl;
  // Initialization
  using std::cout;
  using std::endl;
  std::vector<std::shared_ptr<Application>> res={};
  res.resize(n);
  // Random Distributions:
  std::uniform_real_distribution dist(0.,1.);
  std::uniform_real_distribution noise(-1.,1.);
  for(int i=0;i<n;++i){
    std::string name = "App" + std::to_string(i+1);
    // Generate Random Params
    double R = lR + (uR-lR)*dist(engine);
    double lambda = lL + (uL-lL)*dist(engine);
    double r1 = lr1 + (ur1-lr1)*dist(engine);
    // Compute Other Params
    double D_c = R/15 *(1+0.1*noise(engine));
    double D_e = 6./5*D_c;
    double delta = 0.66016*(1+0.1*noise(engine));
    double m_2 = 1 + 0.1*noise(engine);
    double m_1 = m_2+0.5;
    double gamma = 1.3*(1+0.1*noise(engine));
    //cout << name << " " << "R: " << R << " ,D_c: " << D_c << ", lambda: " << lambda << endl;
    // Create Pointer
    auto appPtr = std::make_shared<Application>(name,D_e,D_c,lambda,delta,
    m_1,m_2,R,r1,gamma);
    res[i] = appPtr;
  }
  return res;
};

// Creates Users
inline
std::vector<User> create_users(int n,const std::shared_ptr<Application>& appPtr,
double lB,double uB,double lM,double uM,int lT,int uT,double lU,double uU,
std::default_random_engine& engine,unsigned& counterID,double var){
  // std::cout << "Creating Users" << std::endl;
  std::vector<User> res={};
  res.resize(n);
  // Random Distributions
  std::uniform_real_distribution dist(0.,1.);
  std::uniform_real_distribution noise(-1.,1.);
  std::uniform_int_distribution<int> distT(lT,uT);
  // Cut users in half for different time usage
  int midN = std::floor(n/2);
  for(int i=0;i<n;++i){
    using std::cout;
    using std::endl;
    unsigned ID = counterID++;
    double D_1 = std::numeric_limits<double>::infinity();
    while(D_1 + appPtr->D_e > appPtr->R) D_1 = appPtr->R*3/4*(1+var*noise(engine));
    double D_2 = D_1 * 0.522/3.211;
    double p_1 = 5*D_1;
    double p_2 = 10*D_2;
    double beta = 0.491 + (0.521-0.491)*dist(engine);
    double M = lM + (uM-lM)*dist(engine);
    double B = lB + (uB-lB)*dist(engine);
    double U = lU + (uU-lU)*dist(engine);
    int T = distT(engine) + (i>=midN)*distT(engine);
    double alpha = 0.65;
    double zeta = 5.3e-5 + (6.0e-5-5.3e-5)*dist(engine);
    double lE = 0.9*appPtr->lambda*T*T*p_2;
    double uE = 1.2*appPtr->lambda*T*T*p_1;
    double E = lE + (uE-lE)*dist(engine);
    User u(ID,D_1,D_2,p_1,p_2,beta,B,alpha,E,M,T,zeta,U);
    u.app = appPtr;
    res[i]=u;
  }
  return res;
};


// Creates Platform and adds Users and Applications
inline
Platform create_platform(const std::string& filename,
std::default_random_engine& engine){
  // Check that filename exists
  std::ifstream file(filename);
  if(not file) {
    std::cerr << "File " << filename << " cannot be found" << std::endl;
    std::exit(1);
  }
  // The json object type
  using json=nlohmann::json;
  // Read from a file
  std::ifstream ifile(filename);
  json data;
  ifile >>data;
  // Read values:
  auto n_users = data["n_users"].get<int>();
  auto n_apps = data["n_apps"].get<int>();
  // auto Ne = data["Ne"].get<int>();
  auto nodes_per_user_ratio = data["nodes_per_user_ratio"].get<double>();
  int Ne = n_users*nodes_per_user_ratio;
  auto c = data["c"].get<double>();
  auto ce = c*0.2;
  auto rMin = data["rMin"].get<double>();
  auto rMax = data["rMax"].get<double>();
  auto equally_dist = data["equally_dist"].get<bool>();
  auto Time = data["T"].get<int>();
  Platform platform(Ne,ce,c,rMin,rMax,Time);

  auto lR = data["lR"].get<double>();
  auto uR = data["uR"].get<double>();
  auto lL = data["lL"].get<double>();
  auto uL = data["uL"].get<double>();
  auto lr1 = data["lr1"].get<double>();
  auto ur1 = data["ur1"].get<double>();
  auto appsVec = create_apps(n_apps,lR,uR,lL,uL,lr1,ur1,engine);

  // Add Applications
  for(auto it=appsVec.cbegin();it!=appsVec.cend();++it){
    platform.apps[(*it)->name]=*it;
  }
  // Create Users for each application
  std::vector<int> n_users_vec(n_apps,0);
  if(equally_dist){
    for(int i=0;i<n_apps;++i) n_users_vec[i]= n_users/n_apps;
  }
  else{
    int usersLeft = n_users;
    for(int i=0;i<n_apps;++i){
      std::uniform_int_distribution d(1,usersLeft-5*(n_apps-i));
      int rand = d(engine);
      while(rand>usersLeft-n_apps) rand=d(engine);
      n_users_vec[i] = (i<n_apps-1) ? rand : usersLeft;
      usersLeft -= rand;
    }
  }
  auto lB = data["lB"].get<double>();
  auto uB = data["uB"].get<double>();
  auto lM = data["lM"].get<double>();
  auto uM = data["uM"].get<double>();
  auto lU = data["lU"].get<double>();
  auto uU = data["uU"].get<double>();
  auto var = data["var"].get<double>();
  int lT = data["lT"].get<int>();
  int uT = data["uT"].get<int>();
  for(int i=0;i<n_apps;++i){
    platform.addUsers(create_users(n_users_vec[i],appsVec[i],
    lB,uB,lM,uM,lT,uT,lU,uU,engine,platform.id_counter,var));
  }
  return platform;
};

inline
Platform create_platform(const std::string& filename,
std::default_random_engine& engine,double var){
  // Check that filename exists
  std::ifstream file(filename);
  if(not file) {
    std::cerr << "File " << filename << " cannot be found" << std::endl;
    std::exit(1);
  }
  // The json object type
  using json=nlohmann::json;
  // Read from a file
  std::ifstream ifile(filename);
  json data;
  ifile >>data;
  // Read values:
  auto n_users = data["n_users"].get<int>();
  auto n_apps = data["n_apps"].get<int>();
  // auto Ne = data["Ne"].get<int>();
  auto nodes_per_user_ratio = data["nodes_per_user_ratio"].get<double>();
  int Ne = n_users*nodes_per_user_ratio;
  auto c = data["c"].get<double>();
  auto ce = c*0.2;
  auto rMin = data["rMin"].get<double>();
  auto rMax = data["rMax"].get<double>();
  auto equally_dist = data["equally_dist"].get<bool>();
  auto Time = data["T"].get<int>();
  Platform platform(Ne,ce,c,rMin,rMax,Time);

  auto lR = data["lR"].get<double>();
  auto uR = data["uR"].get<double>();
  auto lL = data["lL"].get<double>();
  auto uL = data["uL"].get<double>();
  auto lr1 = data["lr1"].get<double>();
  auto ur1 = data["ur1"].get<double>();
  auto appsVec = create_apps(n_apps,lR,uR,lL,uL,lr1,ur1,engine);

  // Add Applications
  for(auto it=appsVec.cbegin();it!=appsVec.cend();++it){
    platform.apps[(*it)->name]=*it;
  }
  // Create Users for each application
  std::vector<int> n_users_vec(n_apps,0);
  if(equally_dist){
    for(int i=0;i<n_apps;++i) n_users_vec[i]= n_users/n_apps;
  }
  else{
    int usersLeft = n_users;
    for(int i=0;i<n_apps;++i){
      std::uniform_int_distribution d(1,usersLeft-5*(n_apps-i));
      int rand = d(engine);
      while(rand>usersLeft-n_apps) rand=d(engine);
      n_users_vec[i] = (i<n_apps-1) ? rand : usersLeft;
      usersLeft -= rand;
    }
  }

  auto lB = data["lB"].get<double>();
  auto uB = data["uB"].get<double>();
  auto lM = data["lM"].get<double>();
  auto uM = data["uM"].get<double>();
  auto lU = data["lU"].get<double>();
  auto uU = data["uU"].get<double>();
  //auto var = data["var"].get<double>();
  int lT = data["lT"].get<int>();
  int uT = data["uT"].get<int>();
  for(int i=0;i<n_apps;++i){
    platform.addUsers(create_users(n_users_vec[i],appsVec[i],
    lB,uB,lM,uM,lT,uT,lU,uU,engine,platform.id_counter,var));
  }
  return platform;
};


inline
Platform create_platform(const std::string& filename,
std::default_random_engine& engine,int n_users,int n_apps){
  // Check that filename exists
  std::ifstream file(filename);
  if(not file) {
    std::cerr << "File " << filename << " cannot be found" << std::endl;
    std::exit(1);
  }
  // The json object type
  using json=nlohmann::json;
  // Read from a file
  std::ifstream ifile(filename);
  json data;
  ifile >>data;
  // Read values:
  // auto n_users = data["n_users"].get<int>();
  // auto n_apps = data["n_apps"].get<int>();
  // auto Ne = data["Ne"].get<int>();
  auto nodes_per_user_ratio = data["nodes_per_user_ratio"].get<double>();
  int Ne = n_users*nodes_per_user_ratio;
  auto c = data["c"].get<double>();
  auto ce = c*0.2;
  auto rMin = data["rMin"].get<double>();
  auto rMax = data["rMax"].get<double>();
  auto equally_dist = data["equally_dist"].get<bool>();
  auto Time = data["T"].get<int>();
  Platform platform(Ne,ce,c,rMin,rMax,Time);

  auto lR = data["lR"].get<double>();
  auto uR = data["uR"].get<double>();
  auto lL = data["lL"].get<double>();
  auto uL = data["uL"].get<double>();
  auto lr1 = data["lr1"].get<double>();
  auto ur1 = data["ur1"].get<double>();
  auto appsVec = create_apps(n_apps,lR,uR,lL,uL,lr1,ur1,engine);

  // Add Applications
  for(auto it=appsVec.cbegin();it!=appsVec.cend();++it){
    platform.apps[(*it)->name]=*it;
  }
  // Create Users for each application
  std::vector<int> n_users_vec(n_apps,0);
  // if(equally_dist){
  //   for(int i=0;i<n_apps;++i) n_users_vec[i]= n_users/n_apps;
  // }
  if(equally_dist){
    for(int i=0;i<n_users;++i){
      auto index = i%n_apps;
      ++n_users_vec[index];
    }
  }
  else{
    int usersLeft = n_users;
    for(int i=0;i<n_apps;++i){
      std::uniform_int_distribution d(1,usersLeft-5*(n_apps-i));
      int rand = d(engine);
      while(rand>usersLeft-n_apps) rand=d(engine);
      n_users_vec[i] = (i<n_apps-1) ? rand : usersLeft;
      usersLeft -= rand;
    }
  }
  auto lB = data["lB"].get<double>();
  auto uB = data["uB"].get<double>();
  auto lM = data["lM"].get<double>();
  auto uM = data["uM"].get<double>();
  auto lU = data["lU"].get<double>();
  auto uU = data["uU"].get<double>();
  auto var = data["var"].get<double>();
  int lT = data["lT"].get<int>();
  int uT = data["uT"].get<int>();
  for(int i=0;i<n_apps;++i){
    platform.addUsers(create_users(n_users_vec[i],appsVec[i],
    lB,uB,lM,uM,lT,uT,lU,uU,engine,platform.id_counter,var));
  }
  return platform;
};

#endif
