#ifndef USER_CLASS_HPP
#define USER_CLASS_HPP

#include<iostream>
#include"Application.hpp"
#include<memory>
#include<limits>
#include<random>
#include<cstdint>
#include<cmath>
#include<algorithm>
#include<chrono>
#include<array>

class User {
public:
  unsigned id;

  /// USER PARAMETERS

  // Demanding time to run deployments on user's device
  double D_1;
  double D_2;
  // Power consumption of user's device to run deployments
  double p_1;
  double p_2;
  // Coeff to convert energy cons. to cost
  double beta;
  // Bandwidth
  double B;
  // Energy/Cost trade-off parameter
  double alpha;
  // Max energy of user's device
  double E;
  // Max memory of user's device
  double M;
  // Coefficient to convert data transfer size to monetary cost
  double zeta;
  // Time users uses app
  int T;
  // Utility
  double U;

  // pointer to application
  std::shared_ptr<Application> app = nullptr;

  // Constructors
  User(unsigned id, double D_1, double D_2, double p_1, double p_2, double beta,
        double B, double alpha, double E, double M, int T, double zeta, double U);
  User() = default;

  // Algorithm 1 - User's Deployment Decision Algorithm
  // Checks that memeory and energy constraints are satisfied and computes the best
  // deployment for the user, hence the deployment that the User is going to pick.
  //
  // INPUTS: r -> extra cost charged by the platform for running second deployment
  // OUTPUTS: two-dimensional array of boolean elements indicating whether The User
  //          chose the corresponding deployment.
  std::array<bool,2> algorithm1(double r) const;

  // Change Deployment Point
  // Computes the value for r that makes the cost for running both deployments equal,
  // meaning that values of r greater and smaller will result in a different deployment
  // chosen by the User.
  // INPUTS: -
  // OUTPUTS: change deployment point of the User.
  double changeDeplPoint() const;

  // Drop Deployment Point
  // Computes the value for r that makes the cost for running the
  // second deployment equal to the utility of the User, meaning that values of
  // r greater than the drop-value will result in the User not chosing deployment 2.
  // INPUTS: -
  // OUTPUTS: drop deployment point of the User.
  double dropDeplPoint() const;

  // Utility methods
  double getV() const;  // returns User's violation
  double getLD() const; // returns User's local delay
  void print() const {std::cout << "User id: " << id << ", user app: "<<
  app->name << std::endl;};
  void printLong() const {
    std::cout << D_1 << " " << D_2 << " "<< p_1 << " " << p_2 << " " << beta << " " << B << " " << alpha << " " << E
    << " " << M << " " << T << " " << U << std::endl;
  }


};

#endif
