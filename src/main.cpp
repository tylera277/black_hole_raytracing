
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


  //Step 1a
  /*
    Calculating the speed of the camera w.r.t. the FIDO at its location
   */
  double beta_speed = ray_trace.calculate_cameras_speed(cameras_orbital_components);
  
  // The angles from which light rays will emanate w.r.t. the
  // cameras local sky.
  double phi_cs = 3.14159/4.0;
  double theta_cs = 3.14159/4.0;

  std::vector<double> cameras_comp_relative_to_fido;
  std::vector<double> cart_comp_inc_ray_cameras_reference;
  std::vector<double> cart_comp_inc_ray_fido_reference;
  std::vector<double> spherical_comp_ray_wrt_fido_reference;
  std::vector<double> rays_canonical_momenta;
  

  //Step 1b
  /*
    Getting the cameras components of its direction of motion w.r.t. the FIDO at its location
    -Big B (r,theta,phi) in the paper
   */
  cameras_comp_relative_to_fido = ray_trace.cameras_comp_of_motion_wrt_fido(); 
  
  
  while(phi_cs < ((2*3.14159)+3.14159/4.0))
    {
      while(theta_cs < ((2*3.14159)+3.14159/4.0))
	{
	  
	  //Step 2
	  /*
	    Computing, in the cameras reference frame, the cartesian components of the unit vector
	    N that points in the direction of the incoming ray.
	   */
	  cart_comp_inc_ray_cameras_reference =
	    ray_trace.cartesian_components_inc_ray_cameras_reference(theta_cs, phi_cs);

	  
	  //Step 3a
	  /*
	    Computing the direction of motion of the incoming ray as measured by the FIDO in cartesian
	    coordinates.
	   */
	  cart_comp_inc_ray_fido_reference = ray_trace.cartesian_components_inc_ray_fido_reference
	    (beta_speed, cart_comp_inc_ray_cameras_reference);

	  //Step 3b
	  /*
	    Converting step 3a into the FIDO's spherical orthonormal basis.
	   */
	  spherical_comp_ray_wrt_fido_reference =
	    ray_trace.spherical_components_inc_ray_fido_reference(cameras_comp_relative_to_fido,
								  cart_comp_inc_ray_fido_reference);

	  
	  //Step 4a
	  /*
	    Computing the rays canonical momenta
	   */
	  rays_canonical_momenta = ray_trace.rays_canonical_momenta(cameras_orbital_components,
								    spherical_comp_ray_wrt_fido_reference);
	  
	  //Step 4b
	  /*
	    Calculating the ray's other two conserved quantities, b, the axial angular momentum and
	    q, the Carter constant.
	   */
	  constants = ray_trace.axial_ang_mom_carter_constant(camera_orbital_components,
							      rays_canonical_momenta);

	  double b = constants.at(0);
	  double q = constants.at(1);

	  //Step 5
	  /*
	    Step for determining if the ray comes from the event horizon (thus is black), or the background
	    celestial sphere.
	   */
	  //(1May2022) Will come back to do this later.

	  
	  //Step 6
	  /*
	    If the ray does come from the celestial sphere, now we compute its point of origin there.
	    -Basically shoot a ray from the camera into the scene at a pair of angles(theta_cs, phi_cs)
	    and propagate very far back in time to see where it ultimately "hits" the celestial sphere,
	    you then take the pixel readings at that location and put them into a log as the pixel readings
	    at its intensity for those pair of initial angles in the camera reference frame.

	    -A future step is to adjust the intensity of the 
	   */

	  

	  

	  
	  theta_cs += (3.14159/180.0);
	}
      phi_cs += (3.14159/180.0);
    }


}
