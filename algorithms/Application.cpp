#include"Application.hpp"

// Constructor
Application::Application(const std::string& name, double D_e, double D_c, double lambda,
  double delta, double m_1, double m_2, double R,
  double r_1,double gamma): name(name), D_e(D_e), D_c(D_c), lambda(lambda),
  delta(delta),m_1(m_1),m_2(m_2),R(R), r_1(r_1), gamma(gamma) {}
