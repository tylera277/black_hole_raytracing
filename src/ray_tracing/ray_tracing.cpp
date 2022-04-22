
#include <Eigen/Dense>

#include "ray_tracing.hpp"


void RayTracing::initialize_cameras_location(std::vector<double> components)
{
  double r_c = components.at(0);
  double theta_c = components.at(1);
  double rho_c = components.at(2);

  double angular_momentum = components.at(3);

  std::vector<double> B;
  B.push_back(0); //radial
  B.push_back(0); //theta
  B.push_back(1); //phi

  //double sigma = 

  
}
