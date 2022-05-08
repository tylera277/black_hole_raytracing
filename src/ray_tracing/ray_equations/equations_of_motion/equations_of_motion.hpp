
class EquationsOfMotion{
private:

public:
  
  // r, theta, a, b, p_r, p_theta
  double dr_deta(double, double, double, double, double, double, double);

  double dtheta_deta(double, double, double, double, double, double, double);

  double dphi_deta_step1(double, double, double, double, double, double, double);

  double dpr_deta_step1(double, double, double, double, double, double, double);

  double dptheta_deta_step1(double, double, double, double, double, double, double);
  
  // r, theta, a, b, p_r, p_theta, step_size

  double dphi_deta(double, double, double, double, double, double, double, double);

  double dpr_deta(double, double, double, double, double, double, double, double);

  double dptheta_deta(double, double, double, double, double, double, double, double);


};

