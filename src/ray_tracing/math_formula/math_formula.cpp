//Formulas needed during the ray tracing calculations

#include <cmath>
#include <iostream>

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




double MathFormula::Omega(double radius, double kerr_metric, double theta)
{
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




double MathFormula::b_zero(double radius, double kerr_metric){
  
  double top = pow(radius,3) - (3*pow(radius,2)) + pow(kerr_metric,2)*radius + pow(kerr_metric,2);
  
  double bottom = kerr_metric * (radius - 1);

  double result = -top/bottom;
  return result;
}





double MathFormula::q_zero(double radius, double kerr_metric){
  double top_first_part = pow(radius,3);
  double top_second_part = pow(radius,3)-6*pow(radius,2)+9*radius-4*pow(kerr_metric,2);

  double bottom = pow(kerr_metric,2)*pow((radius - 1),2);

  double result = -(top_first_part * top_second_part) / bottom;

  return result;
}




double MathFormula::R(double radius, double a, double b, double q){
  double first_part, second_part, E;

  first_part = pow(radius,4)+pow(a,4)+pow(a,2)*pow(b,2)+2*pow(radius,2)*pow(a,2) -2*pow(radius,2)*a*b-
    2*pow(a,3)*b;

  E = pow((b-a),2)+q;
  
  second_part = -pow(radius,2)*E+2*radius*E+pow(a,2)*E;
    

  double result = first_part + second_part;
  
  return result;
}




double MathFormula::R_prime(double radius, double a, double b, double q){
  
  double E = pow((b-a),2)+q;

  double result = 4*pow(radius,3) + 4*radius*pow(a,2)-4*radius*a*b-2*radius*E + 2*E;
  
  return result;
}




double MathFormula::NewtonRaphson(double q, double a, double b,
					       double initial_guess, double termination_criteria){
  
  double current_value = initial_guess;
  double difference=10000;

  while(difference>termination_criteria){
    double R_ = this->R(current_value, a, b, q);
    double R_prime_ = this->R_prime(current_value,a,b,q);

    current_value -= R_/R_prime_;
    
    difference = R_/R_prime_;
    //std::cout << "DIFF:"<< difference<<"\n";
    //std::cout << "VALUE:" << current_value << "\n";
  }
  return current_value;
}
