// Implementing the Runge-Kutta Fehlberg variation as was done
// in the Interstellar paper

// 09May2022

#include <cmath>
#include <vector>
#include <iostream>

#include "ray_equations2.hpp"
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
  double eta_f = -1e10;

  // Im going to try a one step size fits all approach, may need to
  // adjust in future.
  // * For the partial derivative estimates * 
  double partial_der_step_size = 0.01;
  double pdss = partial_der_step_size;

  // Step size for time integration
  double time_step_size = 100;
  double dt = time_step_size;
  
  std::cout << "TIME: " << dt << "\n";
  std::cout << "INITIAL ANGLES: " << theta <<"; " << phi<<"\n";

  while(eta>eta_f){
    ///
    double a1 = dt * EOM.dr_deta(r,theta,kerr_constant,b,q,p_r,p_theta);
    double b1 = dt * EOM.dtheta_deta(r,theta,kerr_constant,b,q,p_r,p_theta);
    double c1 = dt * EOM.dphi_deta(r,theta,kerr_constant,b,q,p_r,p_theta,pdss);
    double d1 = dt * EOM.dpr_deta(r,theta,kerr_constant,b,q,p_r,p_theta,pdss);
    double e1 = dt * EOM.dptheta_deta(r,theta,kerr_constant,b,q,p_r,p_theta,pdss);
    //std::cout << "MEH: "<<a1 <<"; " <<b1<< "; "<<c1 <<"; " <<d1 <<"; " <<e1<< "\n";
    
    //
    double a2 = dt * EOM.dr_deta(r+(1.0/4.0)*a1,theta+(1.0/4.0)*b1,kerr_constant,
				 b, q,
				 p_r+(1.0/4.0)*d1,p_theta+(1.0/4.0)*e1);
    double b2 = dt * EOM.dtheta_deta(r+(1.0/4.0)*a1,theta+(1.0/4.0)*b1,kerr_constant,
				 b, q,
				 p_r+(1.0/4.0)*d1,p_theta+(1.0/4.0)*e1);
    double c2 = dt * EOM.dphi_deta(r+(1.0/4.0)*a1,theta+(1.0/4.0)*b1,kerr_constant,
				 b, q,
				   p_r+(1.0/4.0)*d1,p_theta+(1.0/4.0)*e1, pdss);
    double d2 = dt * EOM.dpr_deta(r+(1.0/4.0)*a1,theta+(1.0/4.0)*b1,kerr_constant,
				   b, q,
				   p_r+(1.0/4.0)*d1,p_theta+(1.0/4.0)*e1, pdss);
    double e2 = dt * EOM.dptheta_deta(r+(1.0/4.0)*a1,theta+(1.0/4.0)*b1,kerr_constant,
				      b, q,
				      p_r+(1.0/4.0)*d1,p_theta+(1.0/4.0)*e1, pdss);
    //std::cout << "MEH2: "<<a2 <<"; " <<b2<< "; "<<c2 <<"; " <<d2 <<"; " <<e2<< "\n";

    //
    double a3 = dt * EOM.dr_deta(r+(3.0/32.0)*a1 + (9.0/32.0)*a2,
				 theta+(3.0/32.0)*b1 + (9.0/32.0)*b2,
				 kerr_constant,
				 b,q,
				 p_r+(3.0/32.0)*d1 + (9.0/32.0)*d2,
				 p_theta+(3.0/32.0)*e1 + (9.0/32.0)*e2);
    
    double b3 = dt * EOM.dtheta_deta(r+(3.0/32.0)*a1 + (9.0/32.0)*a2,
				     theta+(3.0/32.0)*b1 + (9.0/32.0)*b2,
				     kerr_constant,
				     b,q,
				     p_r+(3.0/32.0)*d1 + (9.0/32.0)*d2,
				     p_theta+(3.0/32.0)*e1 + (9.0/32.0)*e2);

    double c3 = dt * EOM.dphi_deta(r+(3.0/32.0)*a1 + (9.0/32.0)*a2,
				   theta+(3.0/32.0)*b1 + (9.0/32.0)*b2,
				   kerr_constant,
				   b,q,
				   p_r+(3.0/32.0)*d1 + (9.0/32.0)*d2,
				   p_theta+(3.0/32.0)*e1 + (9.0/32.0)*e2,
				   pdss);
    
    double d3 = dt * EOM.dpr_deta(r+(3.0/32.0)*a1 + (9.0/32.0)*a2,
				   theta+(3.0/32.0)*b1 + (9.0/32.0)*b2,
				   kerr_constant,
				   b,q,
				   p_r+(3.0/32.0)*d1 + (9.0/32.0)*d2,
				   p_theta+(3.0/32.0)*e1 + (9.0/32.0)*e2,
				   pdss);
    
    double e3 = dt * EOM.dptheta_deta(r+(3.0/32.0)*a1 + (9.0/32.0)*a2,
				   theta+(3.0/32.0)*b1 + (9.0/32.0)*b2,
				   kerr_constant,
				   b,q,
				   p_r+(3.0/32.0)*d1 + (9.0/32.0)*d2,
				   p_theta+(3.0/32.0)*e1 + (9.0/32.0)*e2,
				   pdss);
    //std::cout << "MEH3: "<<a3 <<"; " <<b3<< "; "<<c3 <<"; " <<d3 <<"; " <<e3<< "\n";

    //
    double a4 = dt * EOM.dr_deta(r+(1932.0/2197.0)*a1 + (-7200.0/2197.0)*a2 + (7296.0/2197.0)*a3,
				 theta+(1932.0/2197.0)*b1 + (-7200.0/2197.0)*b2 + (7296.0/2197.0)*b3,
				 kerr_constant, b, q,
				 p_r+(1932.0/2197.0)*d1 + (-7200.0/2197.0)*d2 + (7296.0/2197.0)*d3,
				 p_theta+(1932.0/2197.0)*e1 + (-7200.0/2197.0)*e2 + (7296.0/2197.0)*e3);
    
    double b4 = dt * EOM.dtheta_deta(r+(1932.0/2197.0)*a1 + (-7200.0/2197.0)*a2 + (7296.0/2197.0)*a3,
				 theta+(1932.0/2197.0)*b1 + (-7200.0/2197.0)*b2 + (7296.0/2197.0)*b3,
				 kerr_constant, b, q,
				 p_r+(1932.0/2197.0)*d1 + (-7200.0/2197.0)*d2 + (7296.0/2197.0)*d3,
				 p_theta+(1932.0/2197.0)*e1 + (-7200.0/2197.0)*e2 + (7296.0/2197.0)*e3);
    
    double c4 = dt * EOM.dphi_deta(r+(1932.0/2197.0)*a1 + (-7200.0/2197.0)*a2 + (7296.0/2197.0)*a3,
				    theta+(1932.0/2197.0)*b1 + (-7200.0/2197.0)*b2 + (7296.0/2197.0)*b3,
				    kerr_constant, b, q,
				    p_r+(1932.0/2197.0)*d1 + (-7200.0/2197.0)*d2 + (7296.0/2197.0)*d3,
				   p_theta+(1932.0/2197.0)*e1 + (-7200.0/2197.0)*e2 + (7296.0/2197.0)*e3,
				    pdss);
    
    double d4 = dt * EOM.dpr_deta(r+(1932.0/2197.0)*a1 + (-7200.0/2197.0)*a2 + (7296.0/2197.0)*a3,
				   theta+(1932.0/2197.0)*b1 + (-7200.0/2197.0)*b2 + (7296.0/2197.0)*b3,
				   kerr_constant, b, q,
				   p_r+(1932.0/2197.0)*d1 + (-7200.0/2197.0)*d2 + (7296.0/2197.0)*d3,
				   p_theta+(1932.0/2197.0)*e1 + (-7200.0/2197.0)*e2 + (7296.0/2197.0)*e3,
				   pdss);
    
    double e4 = dt * EOM.dptheta_deta(r+(1932.0/2197.0)*a1 + (-7200.0/2197.0)*a2 + (7296.0/2197.0)*a3,
				      theta+(1932.0/2197.0)*b1 + (-7200.0/2197.0)*b2 + (7296.0/2197.0)*b3,
				      kerr_constant, b, q,
				      p_r+(1932.0/2197.0)*d1 + (-7200.0/2197.0)*d2 + (7296.0/2197.0)*d3,
				      p_theta+(1932.0/2197.0)*e1 + (-7200.0/2197.0)*e2 + (7296.0/2197.0)*e3,
				      pdss);
    //std::cout << "MEH4: "<<a4 <<"; " <<b4<< "; "<<c4 <<"; " <<d4 <<"; " <<e4<< "\n";
	
    /////

    double a5 = dt*EOM.dr_deta(r+(439.0/216.0)*a1 + (-8.0)*a2 + (3680.0/513.0)*a3 + (-854.0/4104.0)*a4,
			       theta+(439.0/216.0)*b1 + (-8.0)*b2 + (3680.0/513.0)*b3 + (-854.0/4104.0)*b4,
			       kerr_constant, b, q,
			       p_r+(439.0/216.0)*d1 + (-8.0)*d2 + (3680.0/513.0)*d3 + (-854.0/4104.0)*d4,
			       p_theta+(439.0/216.0)*e1 + (-8.0)*e2 + (3680.0/513.0)*e3 + (-854.0/4104.0)*e4);

    double b5 = dt*EOM.dtheta_deta(r+(439.0/216.0)*a1 + (-8.0)*a2 + (3680.0/513.0)*a3 + (-854.0/4104.0)*a4,
				   theta+(439.0/216.0)*b1 + (-8.0)*b2 + (3680.0/513.0)*b3 + (-854.0/4104.0)*b4,
				   kerr_constant, b, q,
				   p_r+(439.0/216.0)*d1 + (-8.0)*d2 + (3680.0/513.0)*d3 + (-854.0/4104.0)*d4,
				   p_theta+(439.0/216.0)*e1 + (-8.0)*e2 + (3680.0/513.0)*e3 + (-854.0/4104.0)*e4);
    
    double c5 = dt*EOM.dphi_deta(r+(439.0/216.0)*a1 + (-8.0)*a2 + (3680.0/513.0)*a3 + (-854.0/4104.0)*a4,
				 theta+(439.0/216.0)*b1 + (-8.0)*b2 + (3680.0/513.0)*b3 + (-854.0/4104.0)*b4,
				 kerr_constant, b, q,
				 p_r+(439.0/216.0)*d1 + (-8.0)*d2 + (3680.0/513.0)*d3 + (-854.0/4104.0)*d4,
				 p_theta+(439.0/216.0)*e1 + (-8.0)*e2 + (3680.0/513.0)*e3 + (-854.0/4104.0)*e4,
				 pdss);
				
     double d5 = dt*EOM.dpr_deta(r+(439.0/216.0)*a1 + (-8.0)*a2 + (3680.0/513.0)*a3 + (-854.0/4104.0)*a4,
				theta+(439.0/216.0)*b1 + (-8.0)*b2 + (3680.0/513.0)*b3 + (-854.0/4104.0)*b4,
				kerr_constant, b, q,
				p_r+(439.0/216.0)*d1 + (-8.0)*d2 + (3680.0/513.0)*d3 + (-854.0/4104.0)*d4,
				p_theta+(439.0/216.0)*e1 + (-8.0)*e2 + (3680.0/513.0)*e3 + (-854.0/4104.0)*e4,
				pdss);

     double e5 = dt*EOM.dptheta_deta(r+(439.0/216.0)*a1 + (-8.0)*a2 + (3680.0/513.0)*a3 + (-854.0/4104.0)*a4,
				     theta+(439.0/216.0)*b1 + (-8.0)*b2 + (3680.0/513.0)*b3 + (-854.0/4104.0)*b4,
				     kerr_constant, b, q,
				     p_r+(439.0/216.0)*d1 + (-8.0)*d2 + (3680.0/513.0)*d3 + (-854.0/4104.0)*d4,
				     p_theta+(439.0/216.0)*e1 + (-8.0)*e2 + (3680.0/513.0)*e3 + (-854.0/4104.0)*e4,
				     pdss);

     //std::cout << "MEH5: "<<a5 <<"; " <<b5<< "; "<<c5 <<"; " <<d5 <<"; " <<e5<< "\n";


     ///////
     double a6 = dt*EOM.dr_deta(r+(-8.0/27.0)*a1 + 2*a2 + (-3544.0/2565.0)*a3 + (1859.0/4104.0)*a4 +
				(-11.0/40.0)*a5,
				theta+(-8.0/27.0)*b1 + 2*b2 + (-3544.0/2565.0)*b3 + (1859.0/4104.0)*b4 +
				(-11.0/40.0)*b5,
				kerr_constant, b, q,
				p_r+(-8.0/27.0)*d1 + 2*d2 + (-3544.0/2565.0)*d3 + (1859.0/4104.0)*d4 +
				(-11.0/40.0)*d5,
				p_theta+(-8.0/27.0)*e1 + 2*e2 + (-3544.0/2565.0)*e3 + (1859.0/4104.0)*e4 +
				(-11.0/40.0)*e5);

     double b6 = dt*EOM.dtheta_deta(r+(-8.0/27.0)*a1 + 2*a2 + (-3544.0/2565.0)*a3 + (1859.0/4104.0)*a4 +
				    (-11.0/40.0)*a5,
				    theta+(-8.0/27.0)*b1 + 2*b2 + (-3544.0/2565.0)*b3 + (1859.0/4104.0)*b4 +
				    (-11.0/40.0)*b5,
				    kerr_constant, b, q,
				    p_r+(-8.0/27.0)*d1 + 2*d2 + (-3544.0/2565.0)*d3 + (1859.0/4104.0)*d4 +
				    (-11.0/40.0)*d5,
				    p_theta+(-8.0/27.0)*e1 + 2*e2 + (-3544.0/2565.0)*e3 + (1859.0/4104.0)*e4 +
				    (-11.0/40.0)*e5);
     
     double c6 = dt*EOM.dphi_deta(r+(-8.0/27.0)*a1 + 2*a2 + (-3544.0/2565.0)*a3 + (1859.0/4104.0)*a4 +
				  (-11.0/40.0)*a5,
				  theta+(-8.0/27.0)*b1 + 2*b2 + (-3544.0/2565.0)*b3 + (1859.0/4104.0)*b4 +
				  (-11.0/40.0)*b5,
				  kerr_constant, b, q,
				  p_r+(-8.0/27.0)*d1 + 2*d2 + (-3544.0/2565.0)*d3 + (1859.0/4104.0)*d4 +
				  (-11.0/40.0)*d5,
				  p_theta+(-8.0/27.0)*e1 + 2*e2 + (-3544.0/2565.0)*e3 + (1859.0/4104.0)*e4 +
				  (-11.0/40.0)*e5, pdss);

     double d6 = dt*EOM.dpr_deta(r+(-8.0/27.0)*a1 + 2*a2 + (-3544.0/2565.0)*a3 + (1859.0/4104.0)*a4 +
				 (-11.0/40.0)*a5,
				 theta+(-8.0/27.0)*b1 + 2*b2 + (-3544.0/2565.0)*b3 + (1859.0/4104.0)*b4 +
				 (-11.0/40.0)*b5,
				 kerr_constant, b, q,
				 p_r+(-8.0/27.0)*d1 + 2*d2 + (-3544.0/2565.0)*d3 + (1859.0/4104.0)*d4 +
				 (-11.0/40.0)*d5,
				 p_theta+(-8.0/27.0)*e1 + 2*e2 + (-3544.0/2565.0)*e3 + (1859.0/4104.0)*e4 +
				 (-11.0/40.0)*e5, pdss);
     
     double e6 = dt*EOM.dptheta_deta(r+(-8.0/27.0)*a1 + 2*a2 + (-3544.0/2565.0)*a3 + (1859.0/4104.0)*a4 +
				     (-11.0/40.0)*a5,
				     theta+(-8.0/27.0)*b1 + 2*b2 + (-3544.0/2565.0)*b3 + (1859.0/4104.0)*b4 +
				     (-11.0/40.0)*b5,
				     kerr_constant, b, q,
				     p_r+(-8.0/27.0)*d1 + 2*d2 + (-3544.0/2565.0)*d3 + (1859.0/4104.0)*d4 +
				     (-11.0/40.0)*d5,
				     p_theta+(-8.0/27.0)*e1 + 2*e2 + (-3544.0/2565.0)*e3 + (1859.0/4104.0)*e4 +
				     (-11.0/40.0)*e5, pdss);

     //std::cout << "MEH6: "<<a6 <<"; " <<b6<< "; "<<c6 <<"; " <<d6 <<"; " <<e6<< "\n";


     //Updating the values

     r += (47.0/450.0)*a1 + (0)*a2 + (12.0/25.0)*a3 + (32.0/225.0)*a4 + (1.0/30.0)*a5 + (6.0/25.0)*a6;
     theta += (47.0/450.0)*b1 + (0)*b2 + (12.0/25.0)*b3 + (32.0/225.0)*b4 + (1.0/30.0)*b5 + (6.0/25.0)*b6;
     phi += (47.0/450.0)*c1 + (0)*c2 + (12.0/25.0)*c3 + (32.0/225.0)*c4 + (1.0/30.0)*c5 + (6.0/25.0)*c6;
     p_r += (47.0/450.0)*d1 + (0)*d2 + (12.0/25.0)*d3 + (32.0/225.0)*d4 + (1.0/30.0)*d5 + (6.0/25.0)*d6;
     p_theta += (47.0/450.0)*e1 + (0)*e2 + (12.0/25.0)*e3 + (32.0/225.0)*e4 + (1.0/30.0)*e5 + (6.0/25.0)*e6;

     //std::cout << "VALU: " << r << "; " << theta << "; " << phi << "\n";
     // Maybe calculating the truncation error for this problem
     /*
     double truncation_error = abs((1.0/360.0)*(a1+b1+c1+d1+e1) + (-128.0/4275.0)*(a3+b3+c3+d3+e3)+ \
				   (-2197.0/75240.0)*(a4+b4+c4+d4+e4) + (1.0/50.0)*(a5+b5+c5+d5+e5) + \
				   (2.0/55.0)*(a6+b6+c6+d6+e6));
     */
     double truncation_error = abs((1.0/360.0)*(a1) + (-128.0/4275.0)*(a3)+ \
				   (-2197.0/75240.0)*(a4) + (1.0/50.0)*(a5) + \
				   (2.0/55.0)*(a6));
     double epsilon = 50;
     
     double new_step_size = 0.9 * dt * pow((epsilon/(truncation_error+0.001)),(1.0/5.0));
     //std::cout << "TE: " << truncation_error << "\n";
     //std::cout << "STEPPER: " << new_step_size << "\n";
     //std::cout << "TIME: " << eta_f << "\n";
					      
     if(truncation_error > epsilon)
       {
	 dt = new_step_size;
	 //std::cout << "NEW STEP!\n" ;
       }
     else if(truncation_error <= epsilon)
       {
	 eta -= dt;
	 dt = new_step_size;
       }
     //std::cout << "ANGLES: " << theta << "; " << phi << "\n";
				 
  }

  std::vector<double> angle_results;
  angle_results.push_back(theta);
  angle_results.push_back(phi);

  return angle_results;
  
}
