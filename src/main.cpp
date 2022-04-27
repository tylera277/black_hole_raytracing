
#include <vector>
#include <iostream>

#include "ray_tracing/ray_tracing.hpp"


int main(){
  
  RayTracing ray_trace;

  std::vector<double> cameras_orbital_components;

  // R, theta, phi of cameras location
  cameras_orbital_components.push_back(1e10);
  cameras_orbital_components.push_back(0);
  cameras_orbital_components.push_back(0);

  // Angular momentum of spinning black hole
  cameras_orbital_components.push_back(1);


  
  double beta_speed = ray_trace.calculate_cameras_speed(cameras_orbital_components);
  
  // The angles from which light rays will emanate w.r.t. the
  // cameras local sky.
  double phi_cs = 3.14159/4.0;
  double theta_cs = 3.14159/4.0;

  std::vector<double> cart_comp_inc_ray_cameras_reference;
  std::vector<double> cart_comp_inc_ray_fido_reference;

  
  while(phi_cs < ((2*3.14159)+3.14159/4.0))
    {
      while(theta_cs < ((2*3.14159)+3.14159/4.0))
	{
	  cart_comp_inc_ray_cameras_reference =
	    ray_trace.cartesian_components_inc_ray_cameras_reference(theta_cs, phi_cs);
	  
	  cart_comp_inc_ray_fido_reference = ray_trace.cartesian_components_inc_ray_fido_reference
	    (beta_speed, cart_comp_inc_ray_cameras_reference);

	  
	  theta_cs += 1.0;
	}
      phi_cs += 1.0;
    }


}
