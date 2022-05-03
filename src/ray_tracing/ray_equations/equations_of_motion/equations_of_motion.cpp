
#incluce <cmath>

#include "equations_of_motion.hpp"
#include "../../math_formula/math_formula.hpp"



double EquationsOfMotion::dr_deta(double r, double theta, double a, double b,
				  double q, double p_r, double p_theta){
  MathFormula MF;

  double top = MF.del(r, a, theta);
  double bottom = pow(MF.rho(r,a,theta), 2.0);

  double result top/bottom * p_r;

  return result;
}

double EquationsOfMotion::dtheta_deta(double r, double theta, double a, double b,
				      double q, double p_r, double p_theta){

  double result = 1 / pow(MF.rho(r,a,theta),2) * p_theta;

  return result;

}

double EquationsOfMotion::dphi_deta(double r, double theta, double a, double b,
				    double q, double p_r, double p_theta){

  


}
