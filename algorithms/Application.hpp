#ifndef APPLICATION_CLASS_HPP
#define APPLICATION_CLASS_HPP

#include<string>
#include<iostream>

class Application {
public:
  // Application Name
  std::string name;

  /// APPLICATION PARAMETERS

  // Demanding time to run second 2° partition of 2° deployment on edge servers
  double D_e;
  // Demanding time to run second 2° partition of 2° deployment on cloud VMs
  double D_c;
  // Incoming workload of users running second deployment of application
  double lambda;
  // Data transfer size of 1° partition of 1° and 2° depl from users to platform
  double delta;
  // Memory requirements for 1° and 2° deployment
  double m_1;
  double m_2;
  // Time-response limit
  double R;
  // Cost Parameters
  double r_1; // Cost of 1° depl
  double gamma; // Coeff of extra cost for 2° depl

  // Constructor
  Application(const std::string& name, double D_e, double D_c, double lambda,
  double delta, double m_1, double m_2, double R, double r_1,
  double gamma);
  // Default Constructor
  Application() = default;

  // Print
  void print() const{
    std::cout << "App Name: " << name << std::endl;
  };

};

#endif
