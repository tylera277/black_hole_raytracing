#ifndef RAY_TRACING_H
#define RAY_TRACING_H

#include <vector>


class RayTracing{

private:

public:
  RayTracing(){};

  ~RayTracing(){};

  double calculate_cameras_speed(std::vector<double>);

  std::vector<double> cameras_comp_of_motion_wrt_fido();
  
  std::vector<double> cartesian_components_inc_ray_cameras_reference(double, double);

  std::vector<double> cartesian_components_inc_ray_fido_reference(double, std::vector<double>);

  std::vector<double> spherical_components_inc_ray_fido_reference(std::vector<double>,
								  std::vector<double>);
  
  std::vector<double> rays_canonical_momenta(std::vector<double>,
					     std::vector<double>);
  std::vector<double> axial_ang_mom_carter_constant(std::vector<double>,
						    std::vector<double>);

  bool determine_location_of_beam(double, double, double,std::vector<double>,
				  double);

};

#endif
