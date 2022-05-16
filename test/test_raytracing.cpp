#define CATCH_CONFIG_MAIN
#include "catch.hpp"


#include "../src/ray_tracing/ray_tracing.hpp"

#include <vector>



TEST_CASE("speed test", "[unit]"){

  RayTracing ray_trace;

  SECTION("testing the speed calculation"){

    // r=ANYTHING, theta=0, a=0 
    std::vector<double> camera_components;
    camera_components.push_back(1e1);
    camera_components.push_back(0);
    camera_components.push_back(0);
    camera_components.push_back(1);
  
    REQUIRE(ray_trace.calculate_cameras_speed(camera_components) == Approx(0));

    //////////////////

    //r=0, theta=pi, a=5
    std::vector<double> camera_components1;
    camera_components1.push_back(0);
    camera_components1.push_back(3.14159);
    camera_components1.push_back(0);
    camera_components1.push_back(5);
    
    REQUIRE(ray_trace.calculate_cameras_speed(camera_components1) == Approx(0).margin(1e-5));

    //////////////////

    //r=3, a=0, theta=pi/2
    double radius = 3;
    std::vector<double> camera_components2;
    camera_components2.push_back(radius);
    camera_components2.push_back(3.14159/2.0);
    camera_components2.push_back(0);
    camera_components2.push_back(0);

    double result = pow(radius,1.0/2.0) / pow((pow(radius,2)-2*radius),1.0/2.0);
    REQUIRE(ray_trace.calculate_cameras_speed(camera_components2) == Approx(result).margin(1e-5));
  }

  SECTION("components of unit vector pointing in direction of incoming ray"){

    double theta_cs=3.14159/2.0;
    double phi_cs=0;
    std::vector<double> results;

    results = ray_trace.cartesian_components_inc_ray_cameras_reference(theta_cs, phi_cs);
    
    double x_result=results.at(0), y_result=results.at(1), z_result=results.at(2);

    REQUIRE(x_result==Approx(1));
    REQUIRE(y_result==Approx(0).margin(1e-5));
    REQUIRE(z_result==Approx(0).margin(1e-5));
  }
  
}

TEST_CASE("FIDO's spherical orthonormal basis", "[unit]"){
  SECTION("testing kappas value"){
    

  }

}



TEST_CASE("zero calculations for determining origin of light beam", "[unit]"){
  /*
  RayTracing ray_trace;

  SECTION("testing the speed calculation"){
    
    

  }
  */
  
}
