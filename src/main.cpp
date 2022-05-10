
#include <vector>
#include <iostream>
#include <fstream>

#include "ray_tracing/ray_tracing.hpp"
#include "ray_tracing/ray_equations/ray_equations2.hpp"
//#include "catch.hpp"




int main(){
  
  RayTracing ray_trace;
  RayEquations ray_equations;

  // Dimensions of the picture which I will place the black hole
  // in front of.
  double picture_width=4162;
  double picture_height=4101;

  // File which I am outputting the calculations to
  std::ofstream angles;
  angles.open("../angles_rkf4.csv");

  angles << "theta_cs, phi_cs, theta_prime, phi_prime\n";

  
  std::vector<double> cameras_orbital_components;
  double r_c = 5e6;
  double theta_c = 3.14159/2.0;
  double phi_c = 3.14159;

  double kerr_constant = 1;
  
  // R, theta, phi of cameras location
  cameras_orbital_components.push_back(r_c);
  cameras_orbital_components.push_back(theta_c);
  cameras_orbital_components.push_back(phi_c);

  // Angular momentum of spinning black hole
  cameras_orbital_components.push_back(kerr_constant);


  //Step 1a
  /*
    Calculating the speed of the camera w.r.t. the FIDO at its location
    # Sort of tested #
   */
  double beta_speed = ray_trace.calculate_cameras_speed(cameras_orbital_components);
  //std::cout << "SPEED:" << beta_speed << "\n";
  

  // The angles from which light rays will emanate w.r.t. the
  // cameras local sky.
  double theta_cs = 3.14159/2.0, phi_cs =3.14159/2.0;
  
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
  
  
  while(phi_cs <= ((3.14159/2.0)+200*((2*3.14159)/picture_width)))
    {
      theta_cs =(3.14159/2.0);
      while(theta_cs <= ((3.14159/2.0)+200*(3.14159)/picture_height))
  	{

	  
	  //Step 2
	  /*
	    Computing, in the cameras reference frame, the cartesian components of the unit vector
	    N that points in the direction of the incoming ray.

	    # Sort of tested #
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

	  std::vector<double> constants;
	  
	  constants = ray_trace.axial_ang_mom_carter_constant(cameras_orbital_components,
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

	    -A future step is to adjust the intensity of the light...
	   */

	  // The angles of the point at which the ray ends up on the celestial sphere
	  double theta_prime;
	  double phi_prime;
	  std::vector<double> celestial_sphere_angle;


	  // There is definitely a better way to do this than this monstrosity.
	  std::vector<double> initial_conditions;
	  initial_conditions.push_back(r_c);
	  initial_conditions.push_back(theta_cs);
	  initial_conditions.push_back(phi_cs);

	  initial_conditions.push_back(rays_canonical_momenta.at(1));
	  initial_conditions.push_back(rays_canonical_momenta.at(2));
	  initial_conditions.push_back(rays_canonical_momenta.at(3));
	  
	  /* 
	  initial_conditions.insert(initial_conditions.end(),
				    rays_canonical_momenta.begin(),
				    rays_canonical_momenta.end());
	  */
	  initial_conditions.push_back(kerr_constant);

	  initial_conditions.push_back(b);
	  initial_conditions.push_back(q);
	  
	  celestial_sphere_angle = ray_equations.integrate_ray_equations(initial_conditions);

	  // The angles on the celestial sphere where the ray coming from the camera,
	  // ultimately ended up.
	  theta_prime = celestial_sphere_angle.at(0);
	  phi_prime = celestial_sphere_angle.at(1);
	  std::cout <<"PHI: " << phi_prime << " THETA: " << theta_prime << "\n";
	  

	  // Printing out the angles to a csv file for processing in a
	  // python file
	  angles << theta_cs << "," << phi_cs << ",";
	  angles << theta_prime << "," << phi_prime << "\n";
	  

	  


	  
	  //theta_cs += (100*(3.14159/picture_width));
	  theta_cs += ((3.14159/picture_height));

	}
      //phi_cs += (100*((2*3.14159)/picture_height));
      phi_cs += (((2*3.14159)/picture_width));

     }
  
	 

  angles.close();

}
