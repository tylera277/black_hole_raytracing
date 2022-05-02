//Formulas needed during the ray tracing calculations

#include <cmath>

#include "math_formula.hpp"


double MathFormula::sigma(double radius, double kerr_metric, double theta)
{
  /* Sigma must be in radians */

  double r = radius;
  double a = kerr_metric;

  double del_ = this->del(radius, kerr_metric, theta);
  
  double first_term = pow(pow(r,2) + pow(a,2),2);
  double second_term = pow(a,2) * del_ * pow(sin(theta),2);

  double result = pow(first_term - second_term, 1.0/2.0);

  return result;
}



double MathFormula::del(double radius, double kerr_metric, double theta)
{
  double result = pow(radius,2) - 2*radius + pow(kerr_metric,2);

  return result;
}


double MathFormula::rho(double radius, double kerr_metric, double theta)
{
  double first_term = pow(radius,2) + pow(kerr_metric,2)*pow(cos(theta),2);

  double result = pow(first_term, 1.0/2.0);

  return result;


}


double MathFormula::alpha(double radius, double kerr_metric, double theta)
{
  double rho_ = this->rho(radius, kerr_metric, theta);
  double del_ = this->del(radius, kerr_metric, theta);
  double sigma_ = this->sigma(radius, kerr_metric, theta);
  
  double top = rho_ * pow(del_,1.0/2.0);
  double bottom = sigma_;

  double result = top/bottom;

  return result;
}


double MathFormula::omega(double radius, double kerr_metric, double theta)
{
  double sigma_ = this->sigma(radius, kerr_metric, theta);
  
  double top = 2 * kerr_metric * radius;
  double bottom = pow(sigma_,2);

  double result = top/bottom;

  return result;
}

double MathFormula::omega_bar(double radius, double kerr_metric, double theta)
{
  double sigma_ = this->sigma(radius, kerr_metric, theta);
  double rho_ = this->rho(radius, kerr_metric, theta);

  double top = sigma_ * sin(theta);
  double bottom = rho_;

  double result = top / bottom;

  return result;
}

double MathFormula::Omega(double radius, double kerr_metric, double theta){
  double top = 1;
  double bottom = kerr_metric + pow(radius, 3.0/2.0);

  double result = top / bottom;

  return result;
}

double MathFormula::energy_f(double radius, double kerr_metric, double theta, double phi_comp_little_n){

  double alpha_ = this->alpha(radius, kerr_metric, theta);
  double omega_ = this->omega(radius, kerr_metric, theta);
  double omega_bar_ = this->omega_bar(radius, kerr_metric, theta);

  double result = 1 / (alpha_ + omega_ * omega_bar_ * phi_comp_little_n);

  return result;


}
