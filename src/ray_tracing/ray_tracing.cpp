
#include <Eigen/Dense>
#include <cmath>
#include <cassert>
#include <iostream>

#include "ray_tracing.hpp"
#include "math_formula/math_formula.hpp"

double RayTracing::calculate_cameras_speed(std::vector<double> components)
{

  MathFormula MF;
  
  double r_c = components.at(0);
  double theta_c = components.at(1);
  double rho_c = components.at(2);
  double angular_momentum = components.at(3);

  double kerr_metric = angular_momentum;
  
  double beta = (MF.omega_bar(r_c, kerr_metric, theta_c) * (MF.Omega(r_c, kerr_metric, theta_c) -
							    MF.omega(r_c, kerr_metric, theta_c)))
    / MF.alpha(r_c, kerr_metric, theta_c);

  return beta;
  
}

std::vector<double> RayTracing::cameras_comp_of_motion_wrt_fido(){
  std::vector<double> B;
  B.push_back(0); //radial
  B.push_back(0); //theta
  B.push_back(1); //phi

  return B;

}


std::vector<double> RayTracing::cartesian_components_inc_ray_cameras_reference(double theta_cs, double phi_cs)
{
  double N_x = sin(theta_cs) * cos(phi_cs);
  double N_y = sin(theta_cs) * sin(phi_cs);
  double N_z = cos(theta_cs);

  std::vector<double> components;
  components.push_back(N_x);
  components.push_back(N_y);
  components.push_back(N_z);
  
  return components;


}

std::vector<double> RayTracing::cartesian_components_inc_ray_fido_reference(double beta_speed,
									    std::vector<double> cart_comp_inc_ray_cameras_reference){

  std::vector<double> cart_comp = cart_comp_inc_ray_cameras_reference;
  
  double nfx_top = -pow(1-pow(beta_speed,2),1.0/2.0) * cart_comp.at(0);
  double nfx_bot = 1 - (beta_speed * cart_comp.at(1));
  double nfx = nfx_top / nfx_bot;

  double nfy_top = - cart_comp.at(1) + beta_speed;
  double nfy_bot = 1 - beta_speed * cart_comp.at(1);
  double nfy = nfy_top / nfy_bot;

  double nfz_top = -pow(1-pow(beta_speed,2),1.0/2.0) * cart_comp.at(2);
  double nfz_bot = 1 - beta_speed * cart_comp.at(1);
  double nfz = nfz_top / nfz_bot;

  std::vector<double> nf;
  nf.push_back(nfx);
  nf.push_back(nfy);
  nf.push_back(nfz);

  return nf;

}

std::vector<double> RayTracing::spherical_components_inc_ray_fido_reference(std::vector<double> B,
									    std::vector<double> nf_cartesian){

  double kappa_1 = pow(1-B.at(1),1.0/2.0);
  double kappa_2 = pow(pow(B.at(0),2) + pow(B.at(2),2),1.0/2.0);

  assert(kappa_1 == kappa_2);
  
  std::vector<double> nf_spherical;
  
  double radial_term_first_part = B.at(2)/kappa_1 * nf_cartesian.at(0);
  double radial_term_second_part =  B.at(0) * nf_cartesian.at(1);
  double radial_term_third_part = ((B.at(1)*B.at(0))/kappa_1)*nf_cartesian.at(2);
  double radial_term = radial_term_first_part + radial_term_second_part + radial_term_third_part;

  double theta_term_first_part = B.at(1) * nf_cartesian.at(1);
  double theta_term_second_part = kappa_1 * nf_cartesian.at(2);
  double theta_term = theta_term_first_part + theta_term_second_part;
    
  double phi_term_first_part = -(B.at(0)/kappa_1) * nf_cartesian.at(0);
  double phi_term_second_part = B.at(2) * nf_cartesian.at(1);
  double phi_term_third_part = (B.at(1)*B.at(2))/kappa_1 * nf_cartesian.at(2);
  double phi_term = phi_term_first_part + phi_term_second_part + phi_term_third_part;
  
  nf_spherical.push_back(radial_term);
  nf_spherical.push_back(theta_term);
  nf_spherical.push_back(phi_term);

  return nf_spherical;
    
}
