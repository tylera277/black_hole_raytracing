
#include <vector>


class MathFormula{

  
public:
  
  /////////////////////////////
  // Used throughout at various places
  double sigma(double, double, double);

  double del(double, double, double);

  double rho(double, double, double);

  double alpha(double, double, double);

  double omega(double, double, double);

  double omega_bar(double, double, double);

  double Omega(double, double, double);
  ///////////////////////////////

  /////////////////////////////
  // Used in calculating the rays canonical momenta
  double energy_f(double, double, double, double);
  /////////////////////////////

  /////////////////////////////
  // Used in determing if the ray ends up at the event horizon
  // or celestial sphere. Properties of the ray.
  double b_zero(double, double);
  double q_zero(double,double);

  // Used in finding the root of R(r) which again is used in
  // determining where the ray ends up in the scene.
  double R(double, double, double, double);
  double R_prime(double, double, double, double);
  double NewtonRaphson(double, double, double, double,
					   double);
  /////////////////////////////

};
