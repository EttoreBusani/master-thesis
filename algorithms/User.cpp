#include"User.hpp"

// Constructors
User::User(unsigned id, double D_1, double D_2, double p_1, double p_2, double beta,
      double B, double alpha, double E, double M, int T, double zeta, double U): id(id),
      D_1(D_1), D_2(D_2), p_1(p_1), p_2(p_2), beta(beta), B(B), alpha(alpha), E(E),
      M(M), T(T), zeta(zeta), U(U) {}

// Algorithm 1 - User's Deployment Decision Algorithm
std::array<bool,2> User::algorithm1(double r) const{
  double r_2 = app->r_1+app->gamma*r; // cost of 2° deployment
  // Initialization
  double cost1 = std::numeric_limits<double>::infinity();
  double cost2 = std::numeric_limits<double>::infinity();
  // Check Energy and Memory requirements
  if(app->lambda*(T*T)*p_1<=E && app->m_1<=M && D_1<app->R){
    cost1 = alpha*app->r_1/3600*T + (1-alpha)*beta/3600000*p_1*app->lambda*T*T;
  }
  if(app->lambda*(T*T)*p_2<=E && app->m_2<=M && D_2+app->delta/B+app->D_e < app->R)
    cost2 = alpha*r_2/3600*T + (1-alpha)*beta/3600000*p_2*app->lambda*T*T + T*app->delta*app->lambda*zeta;
  // Choose best deployment for User
  bool x_1=(U/3600*T-cost1 >= 0 && cost1 < cost2);
  bool x_2=(U/3600*T-cost2 >= 0 && cost2 <= cost1);
  //std::cout << "User " << id << "-->" << "x_1 = " << x_1 << ", x_2 = " << x_2 << std::endl;
  return {x_1, x_2};
}

// Change Deployment Point
double User::changeDeplPoint() const{
  double lambda = app->lambda;
  double gamma = app->gamma;
  double delta = app->delta;
  return ((1-alpha)*beta*T*lambda*(p_1-p_2)/1000 -
  zeta*lambda*delta*3600)/(alpha*gamma);
}

// Drop Deployment Poinnt
double User::dropDeplPoint() const{
  double lambda = app->lambda;
  double gamma = app->gamma;
  double delta = app->delta;
  double r_1 = app->r_1;
  return (U-alpha*r_1 + (1-alpha)*beta*p_2*lambda*T/1000 +
  3600*zeta*delta*lambda)/(alpha*gamma);
}

// User's violation
double User::getV() const{
  return D_2 + app->delta/B;
}
// User's local delay
double User::getLD() const{ return app->R-D_2-app->delta/B;}
