#include <vector>


class RayTracing{

private:

public:


  double calculate_cameras_speed(std::vector<double>);

  std::vector<double> cameras_comp_of_motion_wrt_fido();
  
  std::vector<double> cartesian_components_inc_ray_cameras_reference(double, double);

  std::vector<double> cartesian_components_inc_ray_fido_reference(double, std::vector<double>);

  std::vector<double> spherical_components_inc_ray_fido_reference(std::vector<double>,
								  std::vector<double>);
  


};
