
#include <Eigen/Dense>
#include <cmath>

#include "ray_tracing.hpp"
#include "math_formula/math_formula.hpp"

double RayTracing::calculate_cameras_speed(std::vector<double> components)
{

  MathFormula MF;
  
  double r_c = components.at(0);
  double theta_c = components.at(1);
  double rho_c = components.at(2);
  double angular_momentum = components.at(3);

  
  std::vector<double> B;
  B.push_back(0); //radial
  B.push_back(0); //theta
  B.push_back(1); //phi

  double kerr_metric = angular_momentum;
  
  double beta = (MF.omega_bar(r_c, kerr_metric, theta_c) * (MF.Omega(r_c, kerr_metric, theta_c) -
							    MF.omega(r_c, kerr_metric, theta_c)))
    / MF.alpha(r_c, kerr_metric, theta_c);

  return beta;
  
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
