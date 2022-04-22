#include <cmath>

#include "math_formula.hpp"


double MathFormula::sigma(double radius, double kerr_metric,
			  double del, double theta)
{

  double r = radius;
  double a = kerr_metric;
  
  double first_term = pow(pow(r,2) + pow(a,2),2);
  double second_term = pow(a,2) * del * pow(sin(theta),2);

  double final = pow(first_term - second_term, 1.0/2.0);

  return final;
}


double del
