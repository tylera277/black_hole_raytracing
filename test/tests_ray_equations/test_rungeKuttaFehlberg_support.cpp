


#include "test_rungeKuttaFehlberg_support.hpp"


#include <cmath>



double G = 6.67E-11;

double RungeKuttaFehlberg::compute_x_accel(double x1, double y1,
					   double x2, double y2,
					   double center_mass){
  double x_diff = x2 - x1;
  double y_diff = y2 - y1;

  double top = G * center_mass * x_diff;
  double bottom = pow(pow(x_diff,2) + pow(y_diff,2),3.0/2.0);

  double a_x = top/bottom;

  return a_x;


};


double RungeKuttaFehlberg::compute_y_accel(double x1, double y1,
					   double x2, double y2,
					   double center_mass){
  double x_diff = x2 - x1;
  double y_diff = y2 - y1;

  double top = G * center_mass * y_diff;
  double bottom = pow(pow(x_diff,2) + pow(y_diff,2),3.0/2.0);

  double a_y = top/bottom;

  return a_y;


};

