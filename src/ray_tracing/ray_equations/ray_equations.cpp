
#include "ray_equations.hpp"
#include "equations_of_motion/equations_of_motion.hpp"

std::vector<double> RayEquations::integrate_ray_equations(std::vector<double>
							  initial_conditions){

  EquationsOfMotion EOM;
  
  // Unpacking all of the initial conditions
  double r = initial_conditions.at(0);
  double theta = initial_conditions.at(1);
  double phi = initial_conditions.at(2);

  double p_r = initial_conditions.at(3);
  double p_theta = initial_conditions.at(4);
  double p_phi = initial_conditions.at(5);

  double kerr_constant = initial_conditions.at(6);

  double b = initial_conditions.at(7);
  double q = initial_conditions.at(8);
  
  // Time conditions which I am integrating over
  double eta = 0;
  double eta_f = 1e5;
  
  while(eta<eta_f){
  
    double dr_deta = EOM.dr_deta(r,theta,kerr_constant,b,q,p_r,p_theta);
  
  }
  
}
