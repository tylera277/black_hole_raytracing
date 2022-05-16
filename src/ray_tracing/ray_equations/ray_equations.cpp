
#include <iostream>

#include "ray_equations.hpp"
#include "equations_of_motion/equations_of_motion.hpp"

std::vector<double> RayEquations::integrate_ray_equations(std::vector<double>
							  initial_conditions){

  EquationsOfMotion EOM;

  std::vector<double> celestial_sphere_angles;
  
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
  double eta_f = -1e8;

  // Im going to try a one step size fits all approach, may need to
  // adjust in future.
  // * For the partial derivative estimates * 
  double partial_der_step_size = 0.01;

  // Step size for time integration
  double time_step_size = 100;
  

  while(eta>eta_f){

    
    // Initial Slopes
    double a1 = time_step_size*EOM.dr_deta(r,theta,kerr_constant,b,q,p_r,p_theta);
    double b1 = time_step_size*EOM.dtheta_deta(r,theta,kerr_constant,b,q,p_r,p_theta);

    double c1 = time_step_size*EOM.dphi_deta(r,theta,kerr_constant,b,q,p_r,p_theta,
					     partial_der_step_size);

    double d1 = time_step_size*EOM.dpr_deta(r,theta,kerr_constant,b,q,p_r,p_theta,
					    partial_der_step_size);

    double e1 = time_step_size*EOM.dptheta_deta(r,theta,kerr_constant,b,q,p_r,p_theta,
						partial_der_step_size);
    //////

    // Middle slopes
    double a2 = time_step_size * EOM.dr_deta(r+0.5*a1, theta+0.5*b1, kerr_constant,
					     b, q, p_r+0.5*d1, p_theta+0.5*e1);
    
    double b2 = time_step_size * EOM.dtheta_deta(r+0.5*a1, theta+0.5*b1, kerr_constant,
						 b, q, p_r+0.5*d1, p_theta+0.5*e1);

    double c2 = time_step_size * EOM.dphi_deta(r+0.5*a1, theta+0.5*b1, kerr_constant,
					     b, q, p_r+0.5*d1, p_theta+0.5*e1,
					     partial_der_step_size);
    
    double d2 = time_step_size * EOM.dpr_deta(r+0.5*a1, theta+0.5*b1, kerr_constant,
					    b, q, p_r+0.5*d1, p_theta+0.5*e1,
					    partial_der_step_size);
    double e2 = time_step_size * EOM.dptheta_deta(r+0.5*a1, theta+0.5*b1, kerr_constant,
						b, q, p_r+0.5*d1, p_theta+0.5*e1,
						partial_der_step_size);

    
    double a3 = time_step_size * EOM.dr_deta(r+0.5*a2, theta+0.5*b2, kerr_constant,
					     b, q, p_r+0.5*d2, p_theta+0.5*e2);
    
    double b3 = time_step_size * EOM.dtheta_deta(r+0.5*a2, theta+0.5*b2, kerr_constant,
					       b, q, p_r+0.5*d2, p_theta+0.5*e2);

    double c3 = time_step_size * EOM.dphi_deta(r+0.5*a2, theta+0.5*b2, kerr_constant,
					     b, q, p_r+0.5*d2, p_theta+0.5*e2,
					     partial_der_step_size);
    
    double d3 = time_step_size * EOM.dpr_deta(r+0.5*a2, theta+0.5*b2, kerr_constant,
					    b, q, p_r+0.5*d2, p_theta+0.5*e2,
					    partial_der_step_size);
    double e3 = time_step_size * EOM.dptheta_deta(r+0.5*a2, theta+0.5*b2, kerr_constant,
						b, q, p_r+0.5*d2, p_theta+0.5*e2,
						partial_der_step_size);

    // Final slopes

    double a4 = time_step_size * EOM.dr_deta(r+a3, theta+b3, kerr_constant,
					     b, q, p_r+d3, p_theta+e3);
    
    double b4 = time_step_size * EOM.dtheta_deta(r+a3, theta+b3, kerr_constant,
					       b, q, p_r+d3, p_theta+e3);

    double c4 = time_step_size * EOM.dphi_deta(r+a3, theta+b3, kerr_constant,
					     b, q, p_r+d3, p_theta+e3,
					     partial_der_step_size);
    
    double d4 = time_step_size * EOM.dpr_deta(r+a3, theta+b3, kerr_constant,
					    b, q, p_r+d3, p_theta+e3,
					    partial_der_step_size);
    double e4 = time_step_size * EOM.dptheta_deta(r+a3, theta+b3, kerr_constant,
						b, q, p_r+d3, p_theta+e3,
						partial_der_step_size);



    // Updating the values
    r += (1.0/6.0) * (a1 + 2*a2 + 2*a3 + a4);
    theta += (1.0/6.0) * (b1 + 2*b2 + 2*b3 + b4);
    phi += (1.0/6.0) * (c1 + 2*c2 + 2*c3 + c4);
    p_r += (1.0/6.0) * (d1 + 2*d2 + 2*d3 + d4);
    p_theta += (1.0/6.0) * (e1+ 2*e2 + 2*e3 + e4);

    //std::cout << "r: " << r << ";";
    //std::cout << "theta: " << theta << "; ";
    //std::cout << "phi: " << phi << "\n";
    
    //std::cout << "rad: " << (a1 + 2*a2 + 2*a3 + a4) << "; ";// <<"\n";
    //std::cout << "theta: " << (b1 + 2*b2 + 2*b3 + b4) <<"; ";
    //std::cout << "phi: " << (c1 + 2*c2 + 2*c3 + c4) << "\n";
    //std::cout << "pr: " << (d1 + 2*d2 + 2*d3 + d4) << ";";
    //std::cout << "ptheta: " <<(e1+ 2*e2 + 2*e3 + e4) << "\n";
    
    
    eta -= time_step_size;
  }

  celestial_sphere_angles.push_back(theta);
  celestial_sphere_angles.push_back(phi);

  return celestial_sphere_angles;
  
}
