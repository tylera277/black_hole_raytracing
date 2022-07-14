


#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <omp.h>

#include "ray_tracing/ray_tracing.hpp"
#include "ray_tracing/ray_equations/ray_equations2.hpp"
//#include "catch.hpp"





//////////////////////////////////
//     DEFINITIONS/ACRONYMS     //
/*

-FIDO: Fiducial Observer



 */


//////////////////////////////////





int main(){
  
  RayTracing ray_trace;
  RayEquations ray_equations;

  // Dimensions of the picture which the black hole will be placed
  // in front of.

  double picture_width=1201;
  double picture_height=802;

  
  // File which I am outputting the calculations to
  std::ofstream angles;
  //angles.open ("../outputs/angles_rkf4_openmp_JULY.csv");
  //angles.open ("output_redemption.txt", std::ios::out);
  angles.close();
  angles.open ("../outputs/angles_rkf4_openmp_JULY12.csv");
  if(angles){
    angles << "theta_cs, phi_cs, theta_prime, phi_prime\n";
    angles.flush();
  }
  
 
  // User-defined values
  std::vector<double> cameras_orbital_components;
  double r_c = 2;
  double theta_c = 3.14159/50000.0;
  double phi_c = 0;//3.14159;

  double kerr_constant = 0.99;



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
  double theta_cs = 3.14159/2.0, phi_cs = 3.14159;
  //double theta_cs =3.14159/2.0, phi_cs = 0;
  
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

  //while(phi_cs <= ((3.14159)+200*((2*3.14159)/picture_width)))
  
  //std::cout << "THREADS_AVAILABLE:" << omp_get_num_threads() << "\n";

  double initial_theta_cs = 3.14159/2.0;


  
#pragma omp parallel for default(none) shared(picture_width, picture_height, \
					      cameras_orbital_components, kerr_constant, \
					      cameras_comp_relative_to_fido,r_c,beta_speed, \
					      std::cout, angles, initial_theta_cs ) \
  private(cart_comp_inc_ray_fido_reference,cart_comp_inc_ray_cameras_reference,\
	  spherical_comp_ray_wrt_fido_reference,rays_canonical_momenta,\
	  ray_equations,ray_trace,theta_cs, phi_cs)
  
  






    
  for(int horizontal_pixel=0; horizontal_pixel<200; horizontal_pixel++){
    


    phi_cs = 3.14159+horizontal_pixel*((2*3.14159)/picture_width);
    
    theta_cs=3.14159/2.0;
      
      
      for(int vertical_pixel=0; vertical_pixel<200; vertical_pixel++){

	theta_cs = (initial_theta_cs)+vertical_pixel*((3.14159)/picture_height);


  	//{


	  
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
	  
	  //for(int i=0;i<3;i++){
	  // std::cout <<"nf"<<i<<": "<< spherical_comp_ray_wrt_fido_reference.at(i)<<"\n";
	  //}






	  
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

	  // Boolean which will store the status of whether the light beam comes from the
	  // horizon or not. The only other possibility at the moment is the celestial sphere.
	  bool event_horizon;
	  
	  event_horizon = ray_trace.determine_location_of_beam(b, q, kerr_constant,
							       rays_canonical_momenta,r_c);
	  std::cout << "STATUS: " << event_horizon<<"\n";


	  // I need to change the if condition later to one that only specifies
	  // it doesnt come from event horizon.
	  if (event_horizon == 0 || event_horizon == 1){


	    
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
	    if(angles){
	      angles << theta_cs << "," << phi_cs << ",";
	      angles << theta_prime << "," << phi_prime << "\n";
	      angles.flush();
	    }
	  }
	  /*
	  else if(event_horizon==1){
	    std::cout << "MAYBE???\n";
	    // If the light comes from the event horizon, that pixel needs to be black in the
	    // resulting image perceived at the cameras location
	    angles << theta_cs << "," << phi_cs << ",";
	    angles << 0 << "," << 0 << "\n";
	  }
	  */
	  


	  
	  //theta_cs += (100*(3.14159/picture_width));
	  theta_cs += ((3.14159/picture_height));

	}
      //phi_cs += (100*((2*3.14159)/picture_height));
      phi_cs += (((2*3.14159)/picture_width));

     }
  
	 

  angles.close();

}
