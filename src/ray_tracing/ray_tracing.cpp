#include "ray_tracing.hpp"
#include "math_formula/math_formula.hpp"


#include <cmath>
#include <cassert>
#include <iostream>
#include <vector>



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
  double nfz_bot = 1 - (beta_speed * cart_comp.at(1));
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
  double kappa_2 = pow(pow(B.at(0),2) + pow(B.at(2),2), 1.0/2.0);

  assert(kappa_1 == kappa_2);
  
  std::vector<double> nf_spherical;
  
  double radial_term_first_part = B.at(2)/kappa_1 * nf_cartesian.at(0);
  double radial_term_second_part =  B.at(0) * nf_cartesian.at(1);
  double radial_term_third_part = ((B.at(1)*B.at(0))/kappa_1)*nf_cartesian.at(2);
  double radial_term = radial_term_first_part + radial_term_second_part + radial_term_third_part;

  double theta_term_first_part = B.at(1) * nf_cartesian.at(1);
  double theta_term_second_part = kappa_1 * nf_cartesian.at(2);
  double theta_term = theta_term_first_part - theta_term_second_part;
    
  double phi_term_first_part = -(B.at(0)/kappa_1) * nf_cartesian.at(0);
  double phi_term_second_part = B.at(2) * nf_cartesian.at(1);
  double phi_term_third_part = (B.at(1)*B.at(2))/kappa_1 * nf_cartesian.at(2);
  double phi_term = phi_term_first_part + phi_term_second_part + phi_term_third_part;
  
  nf_spherical.push_back(radial_term);
  nf_spherical.push_back(theta_term);
  nf_spherical.push_back(phi_term);

  return nf_spherical;
    
}

std::vector<double> RayTracing::rays_canonical_momenta(std::vector<double> cameras_orbital_components,
						       std::vector<double> spherical_comp_ray_wrt_fido_reference)
{
  
  MathFormula MF;

  double r_c = cameras_orbital_components.at(0);
  double theta_c = cameras_orbital_components.at(1);
  double phi_c = cameras_orbital_components.at(2);
  double angular_momentum = cameras_orbital_components.at(3);

  double kerr_metric = angular_momentum;
  
  double energy_f = MF.energy_f(r_c,
				 kerr_metric,
				 theta_c,
				 spherical_comp_ray_wrt_fido_reference.at(2));


  double p_t = -1;
  
  double p_r = energy_f * (MF.rho(r_c, kerr_metric, theta_c)/pow(MF.del(r_c, kerr_metric, theta_c),1.0/2.0)) \
    * spherical_comp_ray_wrt_fido_reference.at(0);

			   
  double p_theta = energy_f * MF.rho(r_c, kerr_metric, theta_c) * spherical_comp_ray_wrt_fido_reference.at(1);
			
  double p_phi = energy_f * MF.omega_bar(r_c, kerr_metric, theta_c)*spherical_comp_ray_wrt_fido_reference.at(2);

  std::vector<double> ray_momenta;
  ray_momenta.push_back(p_t);
  ray_momenta.push_back(p_r);
  ray_momenta.push_back(p_theta);
  ray_momenta.push_back(p_phi);

  return ray_momenta;
}


std::vector<double> RayTracing::axial_ang_mom_carter_constant(std::vector<double> camera_orbital_comps,
							      std::vector<double> rays_canonical_momenta)
{
  double r_c = camera_orbital_comps.at(0);
  double theta_c = camera_orbital_comps.at(1);
  double phi_c = camera_orbital_comps.at(2);
  double angular_momentum = camera_orbital_comps.at(3);

  double kerr_metric = angular_momentum;
  
  
  // Axial angular momentum
  double b = rays_canonical_momenta.at(3);

  // Carter Constant
  double q = pow(rays_canonical_momenta.at(2),2) + pow(cos(theta_c),2)*(((pow(b,2)/(pow(sin(theta_c),2)))
								  - pow(kerr_metric,2)));

  std::vector<double> constants;
  constants.push_back(b);
  constants.push_back(q);

  return constants;
}

bool RayTracing::determine_location_of_beam(double b, double q, double kerr_metric,
					    std::vector<double> rays_canonical_momenta,
					    double camera_radius){

  MathFormula MF;

  double r1 = 2 * (1 + cos((2.0/3.0)*acos(-kerr_metric)));
  double r2 = 2 * (1 + cos((2.0/3.0)*acos(kerr_metric)));

  // b calculations
  double b_zero_r1 = MF.b_zero(r2, kerr_metric);
  double b_zero_r2 = MF.b_zero(r1, kerr_metric);


  // q calculations
  double q_zero_b = MF.q_zero(b, kerr_metric);

  double p_r = rays_canonical_momenta.at(1);

  /*
  std::cout << "r1: "<< r1 << "\n";
  std::cout << "r2: "<< r2 << "\n";

  std::cout << "p_r: " << rays_canonical_momenta.at(1)<<"\n";
  
  std::cout << "b: "<< b << "\n";
  std::cout << "b_zero_1: "<< b_zero_r1 << "\n";
  std::cout << "b_zero_2: "<< b_zero_r2 << "\n";
  
  
  std::cout << "q: "<< q << "\n";
  std::cout << "q_zero_b: "<< q_zero_b << "\n";
  */
  
  // Conditional checking
  if ((q<q_zero_b) && (b_zero_r1 < b < b_zero_r2))
    {
      
      if(p_r>0){
	return true;
      }
      else if(p_r<0){
	return false;
      }
      
    }
  else{
    std::cout << "ROOTING\n";
    double initial_guess = 500;
    double cease_criteria = 0.05;
    
    double root = MF.NewtonRaphson(q,
				   kerr_metric,
				   b,
				   initial_guess,
				   cease_criteria);
    
    std::cout <<"ROOT: " << root << "\n";

    if(camera_radius >= root){
      return false;
    }
    else{
      return true;
    }
    

  }

  
  

}
