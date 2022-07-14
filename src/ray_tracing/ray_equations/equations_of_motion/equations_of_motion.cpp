
#include <cmath>

#include "equations_of_motion.hpp"
#include "../../math_formula/math_formula.hpp"



double EquationsOfMotion::dr_deta(double r, double theta, double a, double b,
				  double q, double p_r, double p_theta){
  MathFormula MF;

  double top = MF.del(r, a, theta);
  double bottom = pow(MF.rho(r,a,theta), 2.0);

  double result = top/bottom * p_r;

  return result;
}

double EquationsOfMotion::dtheta_deta(double r, double theta, double a, double b,
				      double q, double p_r, double p_theta){
  MathFormula MF;

  double result = (1 / pow(MF.rho(r,a,theta),2)) * p_theta;

  return result;

}

double EquationsOfMotion::dphi_deta_step1(double r, double theta, double a, double b,
					  double q, double p_r, double p_theta){
  MathFormula MF;

  double P = pow(r,2)+pow(a,2) - a*b;
  double del = MF.del(r, a, theta);
  double capital_theta = q - pow(cos(theta),2)*((pow(b,2)/pow(sin(theta),2))-pow(a,2));
  double capital_r = pow(P,2)- del*(pow((b-a),2) + q);
  double rho = MF.rho(r,a,theta);

  double top = capital_r + del * capital_theta;
  double bottom = 2 * del * pow(rho,2);
  double result = top/bottom;

  return result;
}

double EquationsOfMotion::dpr_deta_step1(double r, double theta, double a, double b,
					 double q, double p_r, double p_theta){

  MathFormula MF;

  double P = pow(r,2)+pow(a,2) - a*b;
  double del = MF.del(r, a, theta);
  double capital_theta = q - pow(cos(theta),2)*((pow(b,2)/pow(sin(theta),2))-pow(a,2));
  double capital_r = pow(P,2)- del*(pow((b-a),2) + q);
  double rho = MF.rho(r,a,theta);


  double first_part_top = -del * pow(p_r,2);
  double first_part_bot = 2 * pow(rho,2);
  double first_part = first_part_top / first_part_bot;

  double second_part_top = pow(p_theta,2);
  double second_part_bottom = 2 * pow(rho,2);
  double second_part = second_part_top / second_part_bottom;

  double third_part_top = capital_r + del * capital_theta;
  double third_part_bottom = 2 * del * pow(rho,2);
  double third_part = third_part_top / third_part_bottom;

  double result = first_part - second_part + third_part;

  return result;
}

double EquationsOfMotion::dptheta_deta_step1(double r, double theta, double a, double b,
					     double q, double p_r, double p_theta){

  MathFormula MF;

  double P = pow(r,2)+pow(a,2) - a*b;
  double del = MF.del(r, a, theta);
  double capital_theta = q - pow(cos(theta),2)*((pow(b,2)/pow(sin(theta),2))-pow(a,2));
  double capital_r = pow(P,2)- del*(pow((b-a),2) + q);
  double rho = MF.rho(r,a,theta);


  double first_part_top = -del * pow(p_r,2);
  double first_part_bot = 2 * pow(rho,2);
  double first_part = first_part_top / first_part_bot;

  double second_part_top = pow(p_theta,2);
  double second_part_bottom = 2 * pow(rho,2);
  double second_part = second_part_top / second_part_bottom;

  double third_part_top = capital_r + del * capital_theta;
  double third_part_bottom = 2 * del * pow(rho,2);
  double third_part = third_part_top/third_part_bottom;

  double result = first_part - second_part + third_part;

  return result;
}

double EquationsOfMotion::dphi_deta(double r, double theta, double a, double b,
				    double q, double p_r, double p_theta, double step_size){

  // Doing these 3 steps in order to approximately get the partial derivative w.r.t. b
  double dphi_deta_plus_step = this->dphi_deta_step1(r,theta,a,b+step_size,q,p_r,p_theta);
  double dphi_deta_minus_step = this->dphi_deta_step1(r,theta,a,b-step_size,q,p_r,p_theta);
  double dphi_deta = (dphi_deta_plus_step - dphi_deta_minus_step) / (2 * step_size);

  return (-dphi_deta);
}

double EquationsOfMotion::dpr_deta(double r, double theta, double a, double b,
				    double q, double p_r, double p_theta, double step_size){

// Doing these 3 steps in order to approximately get the partial derivative w.r.t. r
  double dpr_deta_plus_step = this->dpr_deta_step1(r+step_size,theta,a,b,q,p_r,p_theta);
  double dpr_deta_minus_step = this->dpr_deta_step1(r-step_size,theta,a,b,q,p_r,p_theta);
  double dpr_deta = (dpr_deta_plus_step - dpr_deta_minus_step) / (2 * step_size);

  return dpr_deta;
}

double EquationsOfMotion::dptheta_deta(double r, double theta, double a, double b,
				    double q, double p_r, double p_theta, double step_size){
  
  // Doing these 3 steps in order to approximately get the partial derivative w.r.t. theta
  
    double dptheta_deta_plus_step = this->dptheta_deta_step1(r,theta+step_size,a,b,q,p_r,p_theta);
    double dptheta_deta_minus_step = this->dptheta_deta_step1(r,theta-step_size,a,b,q,p_r,p_theta);
    double dptheta_deta = (dptheta_deta_plus_step - dptheta_deta_minus_step) / (2 * step_size);

    return dptheta_deta;


}
