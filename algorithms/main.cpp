#include <mpi.h>
#include"User.hpp"
#include"iostream"
#include"Platform.hpp"
#include<fstream>
#include"GetPot"
#include<chrono>
// #include "gnuplot-iostream.hpp"
#include"ConfigParams.hpp"

void printHelp(){
  using namespace std;
  cout << "_____________________________________________________________________________________________" << endl;
  cout << "Hello! \nThis programm was developed by Ettore Busani as a project for PACS course"
  << " with the supervision \nof Prof. Danilo Ardagna and Hamta Sedghani, academic year 2021/2022." << endl;
  cout << "\nThe code is the implementation of the algorithms we studied to solve the "
  << "optimal allocation \nproblem of a platform having to manage the "
  << "load of Users of different applications." << endl;
  cout << "---------------------------------------------------------------------------------------------" << endl;
  cout << "Please see below the instructions for running the code: \n" << endl;
  cout << "main.cpp creates an instance of the problem using values in 'parameters.json' and solves it\n";
  cout << "using the heuristic approach that we studied and implemented. \n" << endl;
  cout << "- You can change variables 'n_users' and 'n_apps' inside 'parameters.json' to create \n";
  cout << "  different instances of the problem; you can also turn 'equally_dist' to false\n  to have a non uniform distribution of users among apps" << endl;
  cout << "- './main' runs the code in serial" << endl;
  cout << "- Option -v turns on text outputs during programm execution." << endl;
  cout << "- 'mpirun -n' followed by the number of ranks and './main' executes \n  the programm in parallel for better performance."<< endl;
  cout << "  ! Beware that executing the programm with multiple ranks causes the text outputs to be missplaced." << endl;
  cout << "- Option -p can be used to show plots of the solution profit profile of Tabu Search." << endl;
  cout << "- Option -o followed by a number (from 0 up to the number of apps) runs an optimized \n";
  cout << "  version of the algorithm: the lower the number the faster the execution, although the profit might be lower as well." << endl;
  cout << "_____________________________________________________________________________________________\n" << endl;
  exit(1);
}

void writeParams(const Platform& platform,int instance){
  using namespace std;
  int nUsers = platform.users.size();
  int nApps = platform.apps.size();
  // Create Output File
  string nameFile = "configFiles/"+to_string(nUsers)
  +"users"+to_string(nApps)+"apps" + to_string(instance)+".yaml";
  ofstream file(nameFile);
  // Write Users and Apps
  file << "Users: [";
  for(auto [id,u] : platform.users) file << "u"+to_string(id)+",";
  file << "]" << endl;
  file << "Applications: [";
  for(auto [name,app] : platform.apps) file << name+",";
  file << "]" << endl;
  file << "Users_per_Application: {";
  for(auto [name,users] : platform.userMap){
    file << "'"+name+"':[";
    for(auto id : users) file << "u"+to_string(id)+",";
    file << "],";
  }
  file << "}" << endl;
  file << endl;

  // Write System Params
  file << "cloud_VM_cost: " << platform.c << endl;
  file << "edge_VM_cost: " << platform.ce << endl;
  file << "max_edge_VM_number: " << platform.Ne << endl;
  file << "M: " << 10000 << endl;
  file << "r_max: " << platform.rMax << endl;
  file << "r_min: " << platform.rMin << endl;
  file << "T: " << platform.Time << endl;
  file << "edge_demand_vector: {";
  for(auto [name,app] : platform.apps) file << "'"+name+"':" + to_string(app->D_e)+",";
  file << "}" << endl;
  file << "cloud_demand_vector: {";
  for(auto [name,app] : platform.apps) file << "'"+name+"':" + to_string(app->D_c)+",";
  file << "}" << endl;
  file << "Lambda: {";
  for(auto [name,app] : platform.apps) file << "'"+name+"':" + to_string(app->lambda)+",";
  file << "}" << endl;
  file << "R_constraints: {";
  for(auto [name,app] : platform.apps) file << "'"+name+"':" + to_string(app->R)+",";
  file << "}" << endl;
  file << "data_size: {";
  for(auto [name,app] : platform.apps) file << "'"+name+"':" + to_string(app->delta)+",";
  file << "}" << endl;
  file << "gamma: {";
  for(auto [name,app] : platform.apps) file << "'"+name+"':" + to_string(app->gamma)+",";
  file << "}" << endl;
  file << "r_1: {";
  for(auto [name,app] : platform.apps) file << "'"+name+"':" + to_string(app->r_1)+",";
  file << "}" << endl;
  file << endl;

  // Write Users' Params
  file << "user_demand1_vector: {";
  for(auto [id,u] : platform.users) file << "'u" + to_string(id) + "':" + to_string(u.D_1) + ",";
  file << "}" << endl;
  file << "user_demand2_vector: {";
  for(auto [id,u] : platform.users) file << "'u" + to_string(id) + "':" + to_string(u.D_2) + ",";
  file << "}" << endl;
  file << "p_1: {";
  for(auto [id,u] : platform.users) file << "'u" + to_string(id) + "':" + to_string(u.p_1) + ",";
  file << "}" << endl;
  file << "p_2: {";
  for(auto [id,u] : platform.users) file << "'u" + to_string(id) + "':" + to_string(u.p_2) + ",";
  file << "}" << endl;
  file << "beta: {";
  for(auto [id,u] : platform.users) file << "'u" + to_string(id) + "':" + to_string(u.beta) + ",";
  file << "}" << endl;
  file << "network_Bandwidth: {";
  for(auto [id,u] : platform.users) file << "'u" + to_string(id) + "':" + to_string(u.B) + ",";
  file << "}" << endl;
  file << "alpha: {";
  for(auto [id,u] : platform.users) file << "'u" + to_string(id) + "':" + to_string(u.alpha) + ",";
  file << "}" << endl;
  file << "data_cost: {";
  for(auto [id,u] : platform.users) file << "'u" + to_string(id) + "':" + to_string(u.zeta) + ",";
  file << "}" << endl;
  file << "time: {";
  for(auto [id,u] : platform.users) file << "'u" + to_string(id) + "':" + to_string(u.T) + ",";
  file << "}" << endl;
  file << "U: {";
  for(auto [id,u] : platform.users) file << "'u" + to_string(id) + "':" + to_string(u.U) + ",";
  file << "}" << endl;
  file << "s_1: {";
  for(auto [id,u] : platform.users){
    bool s_1 = u.app->lambda*(u.T*u.T)*u.p_1<=u.E && u.app->m_1<=u.M && u.D_1<u.app->R;
    file << "'u" + to_string(id) + "':" + to_string(s_1) + ",";
  }
  file << "}" << endl;
  file << "s_2: {";
  for(auto [id,u] : platform.users){
    bool s_2 = u.app->lambda*(u.T*u.T)*u.p_2<=u.E && u.app->m_2<=u.M && u.D_2+u.app->delta/u.B+u.app->D_e < u.app->R;
    file << "'u" + to_string(id) + "':" + to_string(s_2) + ",";
  }
  file << "}" << endl;

  file.close();
}

int
main(int argc, char **argv)
{

  //----------------------------------------------------------------------------
  // MPI SETUP
  //----------------------------------------------------------------------------
  MPI_Init(&argc, &argv);
  MPI_Comm mpi_comm = MPI_COMM_WORLD;

  int size;
  MPI_Comm_size(mpi_comm, &size);

  int rank;
  MPI_Comm_rank(mpi_comm, &rank);

  std::string filename = "parameters.json";


  // //----------------------------------------------------------------------------
  // // READ CL OPTIONS
  // //----------------------------------------------------------------------------
  // GetPot cl(argc,argv);
  //
  // if(cl.search(2,"-h","-help") && rank==0) printHelp();
  // bool verbose = cl.search("-v");
  // bool plot = cl.search("-p");
  // unsigned seed = 9999;
  // std::default_random_engine engine{seed};
  //
  // //----------------------------------------------------------------------------
  // // SYSTEM SETUP
  // //----------------------------------------------------------------------------
  // Platform platform = create_platform(filename,engine);
  // if(verbose && rank == 0) platform.print();
  // std::vector<std::string> initialOrder = {};
  // for(auto it=platform.apps.cbegin();it!=platform.apps.cend();++it)
  //   initialOrder.push_back(it->second->name);
  //
  // //writeParams(platform);
  //
  // //----------------------------------------------------------------------------
  // // DIVIDE WORKLOAD AMONG RANKS
  // //----------------------------------------------------------------------------
  // auto el_to_permute = cl.follow(initialOrder.size(),"-o"); // For optimized version
  // if(el_to_permute > initialOrder.size() && rank==0){
  //   std::cerr << " Wrong option selected, n. of element to permute must not exceed "
  //   << initialOrder.size() << ",\n please try again lowering the number after '-o'." << std::endl;
  //   std::exit(1);
  // }
  // auto orders = getPermutations(initialOrder,el_to_permute);
  // std::vector<std::vector<std::string>> localOrders={};
  // for(std::size_t i=0;i<orders.size();++i){
  //   if(i%size == rank) localOrders.push_back(orders[i]);
  // }
  //
  // //----------------------------------------------------------------------------
  // // TEST ALGORITHMS
  // //----------------------------------------------------------------------------
  // int n = 5; // size of Elite Solution
  //
  // using json=nlohmann::json;
  // std::ifstream ifile(filename);
  // json data;
  // ifile >>data;
  // // Method method = BEST_PROFIT; // Method for selection of neighbor in Tabu Search.
  // auto tabuAttribute = data["tabuAttribute"].get<int>();
  // auto tabuSize = data["tabuSize"].get<int>();
  // auto method = data["method_selection_neighbors"].get<Method>();
  //
  // auto startPerm = std::chrono::steady_clock::now(); // Start clock
  // int maxIter = data["max_iter"].get<int>();
  // // Each Rank runs the heuristic approach considering its local ordering of applications
  // auto startHeuristic = std::chrono::steady_clock::now();
  // auto eliteSol = platform.algorithm2(n,localOrders,verbose);
  // auto initialSol = platform.algorithm3(eliteSol,verbose);
  // auto startTabu = std::chrono::steady_clock::now();
  // auto sol = platform.algorithm4(initialSol,maxIter,engine,verbose,method,rank,plot,tabuAttribute,tabuSize);
  // auto end = std::chrono::steady_clock::now();
  // std::chrono::duration<double> timeHeuristic = end-startHeuristic;
  // std::chrono::duration<double> timeTabu = end-startTabu;
  //
  // double tabuImprove = (sol.profit-initialSol.profit)/initialSol.profit*100;
  // // Compare Solutions found by all ranks
  // if(rank>0) MPI_Send(&sol.profit,1,MPI_DOUBLE,0,0,mpi_comm);
  // int bestRank = 0;
  // if(rank==0){
  //   double bestProfit = sol.profit;
  //   for(int i=1;i<size;++i){
  //     double currProfit = 0;
  //     // Rank 0 receives all profits and updates the best one
  //     MPI_Recv(&currProfit,1,MPI_DOUBLE,i,0,mpi_comm,MPI_STATUS_IGNORE);
  //     if(currProfit>bestProfit) {
  //       bestRank = i;
  //       bestProfit = currProfit;
  //     }
  //   }
  // }
  // MPI_Bcast(&bestRank,1,MPI_INT,0,mpi_comm); // Broadcast best rank too all ranks.
  //
  // // ----------------------------------------------------------------------------
  // // OUTPUT RESULTS
  // // ----------------------------------------------------------------------------
  // if(rank==bestRank){
  //   auto endPerm = std::chrono::steady_clock::now(); // Stop clock
  //   std::chrono::duration<double> timePerm = endPerm-startPerm;
  //   std::cout << "==============================================================================================================" << std::endl;
  //   std::cout << " \t \t EDGE CLOUD GAME - HEURISTIC APPROACH" << std::endl;
  //   std::cout << "Number of ranks: " << size << std::endl;
  //   std::cout << "Elapsed time: " << timePerm.count() << " sec." << std::endl;
  //   std::cout << "Best Solution found by rank: " << rank << std::endl;
  //   sol.print();
  //   std::cout << "==============================================================================================================" << std::endl;
  //   std::cout << " \t \t EDGE CLOUD GAME - TABU SEARCH" << std::endl;
  //   std::string selection;
  //   std::string attribute;
  //   switch(method){
  //     case RANDOM: selection = "RANDOM"; break;
  //     case BEST_PROFIT: selection = "BEST PROFIT"; break;
  //   }
  //   switch(tabuAttribute){
  //     case 0: attribute = "n. of edge and cloud resources"; break;
  //     case 1: attribute = "n. of resources and users on cloud and edge"; break;
  //   }
  //   std::cout << "Method for selection among neighbors: " << selection << std::endl;
  //   std::cout << "Attribute for tabu list: " << attribute << std::endl;
  //   std::cout << "Size of tabu list: " << tabuSize << ", max Iter: " << maxIter << std::endl;
  //   std::cout << "Elapsed time: " << timeTabu.count() << " sec." << std::endl;
  //   std::cout << "Percentage increase in profit: " << tabuImprove << std::endl;
  //   std::cout << "------------------------------------------------------------------------------------------------------------- \n" << std::endl;
  //
  // }
  // // platform.printLong();
  //
  // MPI_Finalize();



  // // ____________________________________________________________________________
  // // ############################# TEST 1 ######################################
  // GetPot cl(argc,argv);
  //
  // if(cl.search(2,"-h","-help") && rank==0) printHelp();
  // bool verbose = cl.search("-v");
  // bool plot = cl.search("-p");
  //
  //
  //   using json=nlohmann::json;
  //   std::ifstream ifile(filename);
  //   json data = json::parse(ifile);
  //   //ifile >> data;
  //   auto min_n_users = data["min_n_users"].get<int>();
  //   auto max_n_users = data["max_n_users"].get<int>();
  //   auto min_n_app = data["min_n_app"].get<int>();
  //   auto max_n_app = data["max_n_app"].get<int>();
  //   auto user_increment = data["user_increment"].get<int>();
  //   int n_points = (max_n_users-min_n_users)/user_increment+1;
  //   ifile.close();
  //   std::vector<int> nUsers(n_points,min_n_users);
  //   for(int i=0;i<n_points;++i) nUsers[i]+=user_increment*i;
  //   std::vector<int> nApp(max_n_app-min_n_app+1,0);
  //   for(int i=0; i<nApp.size(); ++i) nApp[i] = i+min_n_app;
  //   std::vector<std::vector<std::chrono::duration<double>::rep>> vecTimesHeuristic={};
  //   vecTimesHeuristic.resize(nApp.size());
  //   std::vector<std::vector<std::chrono::duration<double>::rep>> vecTimesTabu={};
  //   vecTimesTabu.resize(nApp.size());
  //   std::vector<std::vector<double>> vecProfit={};
  //   vecProfit.resize(nApp.size());
  //
  //   int n_instances = 10;
  //
  //   // Loop for all the combinations of the number of applications and users
  //   for(int ii=0;ii<nApp.size();++ii){
  //     int n_apps = 5;
  //     std::string outputFile = "outputFiles/heuristic"+std::to_string(n_apps)+"apps.txt";
  //     std::ofstream file(outputFile);
  //     for(int jj=0;jj<nUsers.size();++jj){
  //       unsigned seed = 987654;
  //       std::chrono::duration<double>::rep countHeuristic = 0;
  //       std::chrono::duration<double>::rep countTabu = 0;
  //       double finalProfit = 0;
  //       double avgEdgeNodes = 0;
  //       double avgCloudVMs = 0;
  //       double avgUsersEdge = 0;
  //       double avgUsersCloud = 0;
  //       double avgR = 0;
  //       int n_users = nUsers[jj];
  //       for(int nn=0; nn<n_instances; ++nn){
  //         if(rank==0) std::cout << "Iter: " << ii << "-" << jj << "-" << nn << std::endl;
  //         std::default_random_engine engine{seed};
  //         ++seed;
  //         // CREATE PLATFORM
  //         Platform platform = create_platform(filename,engine,n_users,n_apps);
  //         //writeParams(platform,nn);
  //         if(verbose && rank == 0) platform.print();
  //         std::vector<std::string> initialOrder = {};
  //         for(auto it=platform.apps.cbegin();it!=platform.apps.cend();++it)
  //           initialOrder.push_back(it->second->name);
  //           // DIVIDE WORKLOAD AMONG RANKS
  //             //auto el_to_permute = cl.follow(initialOrder.size(),"-o"); // For optimized version
  //             int el_to_permute = nApp[ii];
  //             if(el_to_permute > initialOrder.size() && rank==0){
  //               std::cerr << " Wrong option selected, n. of element to permute must not exceed "
  //               << initialOrder.size() << ",\n please try again lowering the number after '-o'." << std::endl;
  //               std::exit(1);
  //             }
  //             auto orders = getPermutations(initialOrder,el_to_permute);
  //             std::vector<std::vector<std::string>> localOrders={};
  //             for(std::size_t i=0;i<orders.size();++i){
  //               if(i%size == rank) localOrders.push_back(orders[i]);
  //             }
  //             if(rank==0) std::cout << "Local order: " << localOrders.size() << std::endl;
  //
  //             // TEST ALGORITHMS
  //             int n = 10; // size of Elite Solution
  //             using json=nlohmann::json;
  //             std::ifstream ifile(filename);
  //             json data;
  //             ifile >>data;
  //             // Method method = BEST_PROFIT; // Method for selection of neighbor in Tabu Search.
  //             auto tabuAttribute = data["tabuAttribute"].get<int>();
  //             auto tabuSize = data["tabuSize"].get<int>();
  //             auto method = data["method_selection_neighbors"].get<Method>();
  //             auto maxIter = data["max_iter"].get<int>();
  //
  //             auto startHeuristic = std::chrono::steady_clock::now();
  //             Solution sol;
  //             if(localOrders.size()>0){
  //               auto eliteSol = platform.algorithm2(n,localOrders,verbose);
  //               auto initialSol = platform.algorithm3(eliteSol,verbose);
  //               auto startTabu = std::chrono::steady_clock::now();
  //             //auto sol = platform.algorithm4(initialSol,maxIter,engine,verbose,method,rank,plot,tabuAttribute,tabuSize);
  //               sol = initialSol;
  //             }
  //             auto end = std::chrono::steady_clock::now();
  //             std::chrono::duration<double> timeHeuristic = end-startHeuristic;
  //             // std::chrono::duration<double> timeTabu = end-startTabu;
  //             // auto sol = platform.heuristicApproach(n,localOrders,maxIter,verbose,engine,method,
  //             //   rank,plot,tabuAttribute,tabuSize);
  //
  //             // Compare Solutions found by all ranks
  //             if(rank>0) MPI_Send(&sol.profit,1,MPI_DOUBLE,0,0,mpi_comm);
  //             int bestRank = 0;
  //             if(rank==0){
  //               double bestProfit = sol.profit;
  //               for(int i=1;i<size;++i){
  //                 double currProfit = 0;
  //                 // Rank 0 receives all profits and updates the best one
  //                 MPI_Recv(&currProfit,1,MPI_DOUBLE,i,0,mpi_comm,MPI_STATUS_IGNORE);
  //                 if(currProfit>bestProfit) {
  //                   bestRank = i;
  //                   bestProfit = currProfit;
  //                 }
  //               }
  //               finalProfit += bestProfit/n_instances;
  //             }
  //             MPI_Bcast(&bestRank,1,MPI_INT,0,mpi_comm); // Broadcast best rank too all ranks.
  //             if(rank==bestRank){
  //               avgR += sol.r/n_instances;
  //               for(auto [name,app] : platform.apps){
  //                 avgEdgeNodes += (sol.edgeServers.find(name)->second+0.)/n_instances;
  //                 avgCloudVMs += (sol.cloudServers.find(name)->second+0.)/n_instances;
  //               }
  //               for(auto [id,u] : platform.users){
  //                 avgUsersEdge += (0.+sol.userEdgeAllocation.find(id)->second)/n_instances;
  //                 avgUsersCloud += (0.+sol.userCloudAllocation.find(id)->second)/n_instances;
  //               }
  //             }
  //             MPI_Bcast(&avgEdgeNodes,1,MPI_DOUBLE,bestRank,mpi_comm);
  //             MPI_Bcast(&avgCloudVMs,1,MPI_DOUBLE,bestRank,mpi_comm);
  //             MPI_Bcast(&avgUsersEdge,1,MPI_DOUBLE,bestRank,mpi_comm);
  //             MPI_Bcast(&avgUsersCloud,1,MPI_DOUBLE,bestRank,mpi_comm);
  //             MPI_Bcast(&avgR,1,MPI_DOUBLE,bestRank,mpi_comm);
  //             if(rank==0) std::cout << "best rank: " << bestRank << "\n" << std::endl;
  //
  //       /// UPDATE VECTOR OF TIMES FOR PERFORMANCE ANALYSIS
  //         countHeuristic += timeHeuristic.count()/n_instances;
  //         // countTabu += std::log(timeTabu.count())/n_instances;
  //
  //       }
  //       // vecTimesHeuristic[ii].push_back(countHeuristic);
  //       // vecTimesTabu[ii].push_back(countTabu);
  //       vecProfit[ii].push_back(finalProfit);
  //       // if(rank==0) {
  //       //   file << std::to_string(n_users) + " " + std::to_string(finalProfit)+ " " +
  //       //   std::to_string(countHeuristic) + " " + std::to_string(avgEdgeNodes)+" "+std::to_string(avgCloudVMs)
  //       //   +" "+std::to_string(avgUsersEdge)+" "+std::to_string(avgUsersCloud)+" "+std::to_string(avgR)<< std::endl;
  //       // }
  //     }
  //   }

    // if(rank==0){
    //   std::string title = "'Execution Time of Tabu Search'";
    //   Gnuplot gp;
    //   gp << "TITLE = " << title << std::endl;
    //   gp << "set title TITLE" << std::endl;
    //   gp << "set key top left " << std::endl;
    //   gp << "set xlabel 'n. of users'" << std::endl;
    //   gp << "set ylabel 'time (sec) - log scale'" << std::endl;
    //   gp << "plot";
    //   for(int i=0; i<vecTimesTabu.size();++i){
    //     std::string legend = "'" + std::to_string(nApp[i]) + " apps'";
    //     gp << gp.file1d(std::tie(nUsers,vecTimesTabu[i])) << "w lp lw 2 title " << legend  << ",";
    //   }
    //   gp << std::endl;
    // }
    //
    // if(rank==0){
    //   std::string title = "'Execution Time of Heuristic Approach with 6 ranks'";
    //   Gnuplot gp;
    //   gp << "TITLE = " << title << std::endl;
    //   gp << "set title TITLE" << std::endl;
    //   gp << "set key top left " << std::endl;
    //   gp << "set xlabel 'n. of users'" << std::endl;
    //   gp << "set ylabel 'time (sec) - log scale'" << std::endl;
    //   gp << "plot";
    //   for(int i=0; i<vecTimesHeuristic.size();++i){
    //     std::string legend = "'" + std::to_string(nApp[i]) + " apps'";
    //     gp << gp.file1d(std::tie(nUsers,vecTimesHeuristic[i])) << "w lp lw 2 title " << legend  << ",";
    //   }
    //   gp << std::endl;
    // }

    // if(rank==0){
    //   std::string title = "'Total Platform Profit (5 applications)'";
    //   Gnuplot gp;
    //   gp << "TITLE = " << title << std::endl;
    //   gp << "set title TITLE" << std::endl;
    //   gp << "set key top left " << std::endl;
    //   gp << "set xlabel 'n. of users'" << std::endl;
    //   gp << "set ylabel 'Platform Profit'" << std::endl;
    //   gp << "plot";
    //   for(int i=0; i<vecProfit.size();++i){
    //     std::string legend = "'" + std::to_string(nApp[i]) + "-perm profit'";
    //     gp << gp.file1d(std::tie(nUsers,vecProfit[i])) << "w lp lw 2 title " << legend  << ",";
    //   }
    //   gp << std::endl;
    // }
    //
    // MPI_Finalize();





  // ----------------------------------------------------------------------------
  // READ CL OPTIONS
  // ----------------------------------------------------------------------------
  GetPot cl(argc,argv);

  if(cl.search(2,"-h","-help") && rank==0) printHelp();
  bool verbose = cl.search("-v");
  bool plot = cl.search("-p");

  //----------------------------------------------------------------------------
  // SYSTEM SETUP
  //----------------------------------------------------------------------------
  std::vector<double> tabuVec00(9,0);
  std::vector<double> tabuVec01(9,0);
  std::vector<double> tabuVec10(9,0);
  std::vector<double> tabuVec11(9,0);
  std::vector<double> varVec(9,0);
  std::vector<std::string> varStr(9,"");
  int n_instances = 10;
      std::string outputFile = "outputFiles/tabu.txt";
      std::ofstream file(outputFile);
  for(int ii=0; ii<9; ++ii){
    unsigned seed=9999;
    double var = (ii+1.)/10;
    varVec[ii]=var*100;
    varStr[ii]=std::to_string(var*100)+"%";
    double averageTabu = 0;
    for(int j=0; j<n_instances; ++j) {
      std::default_random_engine engine{seed};
      seed++;
      if(rank==0) std::cout << "Iter: " << ii << "-" << j << std::endl;
      Platform platform = create_platform(filename,engine,var);
      if(rank==0) platform.print();
      if(verbose && rank == 0) platform.print();
      std::vector<std::string> initialOrder = {};
      for(auto it=platform.apps.cbegin();it!=platform.apps.cend();++it)
        initialOrder.push_back(it->second->name);

  //----------------------------------------------------------------------------
  // DIVIDE WORKLOAD AMONG RANKS
  //----------------------------------------------------------------------------
      auto el_to_permute = cl.follow(initialOrder.size(),"-o"); // For optimized version
      if(el_to_permute > initialOrder.size() && rank==0){
        std::cerr << " Wrong option selected, n. of element to permute must not exceed "
        << initialOrder.size() << ",\n please try again lowering the number after '-o'." << std::endl;
        std::exit(1);
      }
      auto orders = getPermutations(initialOrder,el_to_permute);
      std::vector<std::vector<std::string>> localOrders={};
      for(std::size_t i=0;i<orders.size();++i){
        if(i%size == rank) localOrders.push_back(orders[i]);
      }


  //----------------------------------------------------------------------------
  // TEST ALGORITHMS
  //----------------------------------------------------------------------------
        int n = 5; // size of Elite Solution

        using json=nlohmann::json;
        std::ifstream ifile(filename);
        json data;
        ifile >>data;
        // Method method = BEST_PROFIT; // Method for selection of neighbor in Tabu Search.
        auto tabuAttribute = data["tabuAttribute"].get<int>();
        auto tabuSize = data["tabuSize"].get<int>();
        auto method = data["method_selection_neighbors"].get<Method>();
        auto maxIter = data["max_iter"].get<int>();

        auto startPerm = std::chrono::steady_clock::now(); // Start clock
        // Each Rank runs the heuristic approach considering its local ordering of applications
        //Solution sol = platform.heuristicApproach(n,localOrders,100,verbose,engine,method,rank,plot,tabuAttribute,tabuSize);
        auto eliteSol = platform.algorithm2(n,localOrders,verbose);
        auto initialSol = platform.algorithm3(eliteSol,verbose);

        auto sol00 = platform.algorithm4(initialSol,maxIter,engine,verbose,RANDOM,rank,plot,0,tabuSize);
        auto sol01 = platform.algorithm4(initialSol,maxIter,engine,verbose,RANDOM,rank,plot,1,tabuSize);
        auto sol10 = platform.algorithm4(initialSol,maxIter,engine,verbose,BEST_PROFIT,rank,plot,0,tabuSize);
        auto sol11 = platform.algorithm4(initialSol,maxIter,engine,verbose,BEST_PROFIT,rank,plot,1,tabuSize);

        double tabuImprove00 = (sol00.profit-initialSol.profit)/initialSol.profit*100;
        double tabuImprove01 = (sol01.profit-initialSol.profit)/initialSol.profit*100;
        double tabuImprove10 = (sol10.profit-initialSol.profit)/initialSol.profit*100;
        double tabuImprove11 = (sol11.profit-initialSol.profit)/initialSol.profit*100;

        // Compare Solutions found by all ranks
        if(rank>0) MPI_Send(&initialSol.profit,1,MPI_DOUBLE,0,0,mpi_comm);
        int bestRank = 0;
        if(rank==0){
          double bestProfit = initialSol.profit;
          for(int i=1;i<size;++i){
            double currProfit = 0;
            // Rank 0 receives all profits and updates the best one
            MPI_Recv(&currProfit,1,MPI_DOUBLE,i,0,mpi_comm,MPI_STATUS_IGNORE);
            if(currProfit>bestProfit) {
              bestRank = i;
              bestProfit = currProfit;
            }
          }
        }
        MPI_Bcast(&bestRank,1,MPI_INT,0,mpi_comm); // Broadcast best rank too all ranks.

        MPI_Bcast(&tabuImprove00,1,MPI_DOUBLE,bestRank,mpi_comm);
        MPI_Bcast(&tabuImprove01,1,MPI_DOUBLE,bestRank,mpi_comm);
        MPI_Bcast(&tabuImprove10,1,MPI_DOUBLE,bestRank,mpi_comm);
        MPI_Bcast(&tabuImprove11,1,MPI_DOUBLE,bestRank,mpi_comm);

        if(rank==0) std::cout << "Tabu improve: " << tabuImprove00 << ", "
        << tabuImprove01 << ", "<< tabuImprove10 << ", "<< tabuImprove11 << "\n " <<  std::endl;

        tabuVec00[ii]+=tabuImprove00/n_instances;
        tabuVec01[ii]+=tabuImprove01/n_instances;
        tabuVec10[ii]+=tabuImprove10/n_instances;
        tabuVec11[ii]+=tabuImprove11/n_instances;

      }
      if(rank==0) {
        file << std::to_string(tabuVec00[ii]) + " " + std::to_string(tabuVec01[ii]) + " "
        + std::to_string(tabuVec10[ii]) + " " + std::to_string(tabuVec11[ii]) << std::endl;
      }
    }

    if(rank==0){
      std::string title = "'Tabu Search Profit Increase'";
      Gnuplot gp;
      gp << "TITLE = " << title << std::endl;
      gp << "set title TITLE" << std::endl;
      gp << "set key top left " << std::endl;
      gp << "set xlabel 'variance of users in %'" << std::endl;
      gp << "set ylabel 'profit % increase '" << std::endl;
      gp << "set xlabel font ',12'" << std::endl;
      gp << "set ylabel font ',12'" << std::endl;
      gp << "plot"
      << gp.file1d(std::tie(varVec,tabuVec00)) <<
      "w lp lw 2 title 'Neighbor selection type 0, tabu list attribute type 0',"
      << gp.file1d(std::tie(varVec,tabuVec01)) <<
      "w lp lw 2 title 'Neighbor selection type 0, tabu list attribute type 1',"
      << gp.file1d(std::tie(varVec,tabuVec10)) <<
      "w lp lw 2 title 'Neighbor selection type 1, tabu list attribute type 0',"
      << gp.file1d(std::tie(varVec,tabuVec11)) <<
      "w lp lw 2 title 'Neighbor selection type 1, tabu list attribute type 1'"
      << std::endl;
      gp << "set xticks (10%, 20%, 30%, 40%, 50%, 60%, 70%, 80%, 90%)" << std::endl;
    }

  // // ----------------------------------------------------------------------------
  // // OUTPUT RESULTS
  // // ----------------------------------------------------------------------------
  // if(rank==bestRank){
  //   auto endPerm = std::chrono::steady_clock::now(); // Stop clock
  //   std::chrono::duration<double> timePerm = endPerm-startPerm;
  //   std::cout << "==============================================================================================================" << std::endl;
  //   std::cout << " \t \t EDGE CLOUD GAME - HEURISTIC APPROACH" << std::endl;
  //   std::cout << "Number of ranks: " << size << std::endl;
  //   std::cout << "Elapsed time: " << timePerm.count() << " sec." << std::endl;
  //   std::cout << "Best Solution found by rank: " << rank << std::endl;
  //   sol.print();
  // }

  MPI_Finalize();



  return 0;
}
