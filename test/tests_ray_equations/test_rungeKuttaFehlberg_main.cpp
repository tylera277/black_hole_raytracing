// Testing the RungeKuttaFehlberg method
// to see if I have the theory correctly understood/
// implemented in my main simulation program

// I will be testing using the motion of the
// planets, inputting initial conditions and
// then checking with actual position and velocity
// values to check for agreement.


#include <iostream>
#include <fstream>
#include <cmath>

#include "test_rungeKuttaFehlberg_support.hpp"


int main(){

  RungeKuttaFehlberg rkf;

  //Initial Conditions
  double center_mass_sun = 1.989E30;
  double center_mass_earth = 5.972E24;
  
  double x1 = 0;
  double y1 =  0;
  double vx1 = 0;
  double vy1 = 0;

  double x2 = -2.6127E10;
  double y2 =  1.4477E11;
  double vx2 = -2.9812E4;
  double vy2 = -5.4009E3;

  double dt = 3600;

  double eta_f = 86400.0*365.0;
  double eta = 0;
  //////////////////////

  // File where positions will be sent
  //std::ofstream positions_output;
  //positions_output.open("output.csv");

  //positions_output << " earth_x_pos, earth_y_pos, sun_x_pos, sun_y_pos \n";


  while(eta<eta_f){
  

    double a1_x = dt * vx1;
    double a1_y = dt * vy1;
  
    double b1_x = dt * vx2;
    double b1_y = dt * vy2;
    
    double a1_vx = dt * rkf.compute_x_accel(x1,y1,x2,y2,center_mass_earth);
    double a1_vy = dt * rkf.compute_y_accel(x1,y1,x2,y2,center_mass_earth);

    double b1_vx = dt * rkf.compute_x_accel(x2,y2,x1,y1,center_mass_sun);
    double b1_vy = dt * rkf.compute_y_accel(x2,y2,x1,y1,center_mass_sun);

    //std::cout << "A1:" << a1_vx << ", " << a1_vy << "\n";
    //std::cout << "B1:" << b1_vx << ", " << b1_vy << "\n";

    //////////////////

  
    double a2_x = dt * (vx1 + (1.0/4.0)*a1_vx);
    double a2_y = dt * (vy1 + (1.0/4.0)*a1_vy);
  
    double b2_x = dt * (vx2 + (1.0/4.0)*b1_vx);
    double b2_y = dt * (vy2 + (1.0/4.0)*b1_vy);
    
    double a2_vx = dt * rkf.compute_x_accel(x1+(1.0/4.0)*a1_x, y1+(1.0/4.0)*a1_y,
					    x2+(1.0/4.0)*b1_x, y2+(1.0/4.0)*b1_y,
					    center_mass_earth);
    double a2_vy = dt * rkf.compute_y_accel(x1+(1.0/4.0)*a1_x, y1+(1.0/4.0)*a1_y,
					    x2+(1.0/4.0)*b1_x, y2+(1.0/4.0)*b1_y,
					    center_mass_earth);

    double b2_vx = dt * rkf.compute_x_accel(x2+(1.0/4.0)*b1_x, y2+(1.0/4.0)*b1_y,
					    x1+(1.0/4.0)*a1_x, y1+(1.0/4.0)*a1_y,
					    center_mass_sun);
    double b2_vy = dt * rkf.compute_y_accel(x2+(1.0/4.0)*b1_x, y2+(1.0/4.0)*b1_y,
					    x1+(1.0/4.0)*a1_x, y2+(1.0/4.0)*a1_y,
					    center_mass_sun);

  
    ///////////////////

  
    double a3_x = dt * (vx1 + (3.0/32.0)*a1_vx + (9.0/32.0)*a2_vx);
    double a3_y = dt * (vy1 + (3.0/32.0)*a1_vy + (9.0/32.0)*a2_vy);
    
    double b3_x = dt * (vx2 + (3.0/32.0)*b1_vx + (9.0/32.0)*b2_vx);
    double b3_y = dt * (vy2 + (3.0/32.0)*b1_vy + (9.0/32.0)*b2_vy);

    double a3_vx = dt * rkf.compute_x_accel(x1+(3.0/32.0)*a1_x + (9.0/32.0)*a2_x,
					    y1+(3.0/32.0)*a1_y + (9.0/32.0)*a2_y,
					    x2+(3.0/32.0)*b1_x + (9.0/32.0)*b2_x,
					    y2+(3.0/32.0)*b1_y + (9.0/32.0)*b2_y,
					    center_mass_earth);
    double a3_vy = dt * rkf.compute_y_accel(x1+(3.0/32.0)*a1_x + (9.0/32.0)*a2_x,
					    y1+(3.0/32.0)*a1_y + (9.0/32.0)*a2_y,
					    x2+(3.0/32.0)*b1_x + (9.0/32.0)*b2_x,
					    y2+(3.0/32.0)*b1_y + (9.0/32.0)*b2_y,
					    center_mass_earth);

    double b3_vx = dt * rkf.compute_x_accel(x2+(3.0/32.0)*b1_x + (9.0/32.0)*b2_x,
					    y2+(3.0/32.0)*b1_y + (9.0/32.0)*b2_y,
					    x1+(3.0/32.0)*a1_x + (9.0/32.0)*a2_x,
					    y1+(3.0/32.0)*a1_y + (9.0/32.0)*a2_y,
					    center_mass_sun);
    double b3_vy = dt * rkf.compute_y_accel(x2+(3.0/32.0)*b1_x + (9.0/32.0)*b2_x,
					    y2+(3.0/32.0)*b1_y + (9.0/32.0)*b2_y,
					    x1+(3.0/32.0)*a1_x + (9.0/32.0)*a2_x,
					    y1+(3.0/32.0)*a1_y + (9.0/32.0)*a2_y,
					    center_mass_sun);
    
    //////////////////
    
    double a4_x = dt * (vx1 + (1932.0/2197.0)*a1_vx + (-7200.0/2197.0)*a2_vx +
			(7296.0/2197.0)*a3_vx);
    double a4_y = dt * (vy1 + (1932.0/2197.0)*a1_vy + (-7200.0/2197.0)*a2_vy +
			(7296.0/2197.0)*a3_vy);
    
    double b4_x = dt * (vx2 + (1932.0/2197.0)*b1_vx + (-7200.0/2197.0)*b2_vx +
			(7296.0/2197.0)*b3_vx);
    double b4_y = dt * (vy2 + (1932.0/2197.0)*b1_vy + (-7200.0/2197.0)*b2_vy +
			(7296.0/2197.0)*b3_vy);
    
    double a4_vx = dt * rkf.compute_x_accel(x1+(1932.0/2197.0)*a1_x +
					    (-7200.0/2197.0)*a2_x +
					    (7296.0/2197.0)*a3_x,
					    
					    y1+(1932.0/2197.0)*a1_y +
					    (-7200.0/2197.0)*a2_y +
					    (7296.0/2197.0)*a3_y,
					    
					    x2+(1932.0/2197.0)*b1_x +
					    (-7200.0/2197.0)*b2_x +
					    (7296.0/2197.0)*b3_x,
					    
					    y2+(1932.0/2197.0)*b1_y +
					    (-7200.0/2197.0)*b2_y +
					    (7296.0/2197.0)*b3_y,
					    
					    center_mass_earth);
    double a4_vy = dt * rkf.compute_y_accel(x1+(1932.0/2197.0)*a1_x +
					    (-7200.0/2197.0)*a2_x +
					    (7296.0/2197.0)*a3_x,
					    
					    y1+(1932.0/2197.0)*a1_y +
					    (-7200.0/2197.0)*a2_y +
					    (7296.0/2197.0)*a3_y,
					    
					    x2+(1932.0/2197.0)*b1_x +
					    (-7200.0/2197.0)*b2_x +
					    (7296.0/2197.0)*b3_x,
					    
					    y2+(1932.0/2197.0)*b1_y +
					    (-7200.0/2197.0)*b2_y +
					    (7296.0/2197.0)*b3_y,
					    
					    center_mass_earth);
    
    double b4_vx = dt * rkf.compute_x_accel(x2+(1932.0/2197.0)*b1_x +
					   (-7200.0/2197.0)*b2_x +
					    (7296.0/2197.0)*b3_x,
					    
					    y2+(1932.0/2197.0)*b1_y +
					    (-7200.0/2197.0)*b2_y +
					    (7296.0/2197.0)*b3_y,
					    
					    x1+(1932.0/2197.0)*a1_x +
					    (-7200.0/2197.0)*a2_x +
					    (7296.0/2197.0)*a3_x,
					    
					    y1+(1932.0/2197.0)*a1_y +
					    (-7200.0/2197.0)*a2_y +
					    (7296.0/2197.0)*a3_y,
					    
					    center_mass_sun);
    
    double b4_vy = dt * rkf.compute_y_accel(x2+(1932.0/2197.0)*b1_x +
					    (-7200.0/2197.0)*b2_x +
					    (7296.0/2197.0)*b3_x,
					    
					    y2+(1932.0/2197.0)*b1_y +
					    (-7200.0/2197.0)*b2_y +
					    (7296.0/2197.0)*b3_y,
					    
					    x1+(1932.0/2197.0)*a1_x +
					    (-7200.0/2197.0)*a2_x +
					    (7296.0/2197.0)*a3_x,
					    
					    y1+(1932.0/2197.0)*a1_y +
					    (-7200.0/2197.0)*a2_y +
					    (7296.0/2197.0)*a3_y,
					    
					    center_mass_sun);
    
    
    //////////////////
    
    
    double a5_x = dt * (vx1 + (439.0/216.0)*a1_vx + (-8.0)*a2_vx +
			(3680.0/513.0)*a3_vx + (-845.0/4104.0)*a4_vx);
    
    double a5_y = dt * (vy1 + (439.0/216.0)*a1_vy + (-8.0)*a2_vy +
			(3680.0/513.0)*a3_vy + (-845.0/4104.0)*a4_vy);
    
    double b5_x = dt * (vx2 + (439.0/216.0)*b1_vx + (-8.0)*b2_vx +
			(3680.0/513.0)*b3_vx + (-845.0/4104.0)*b4_vx);
    
    double b5_y = dt * (vy2 + (439.0/216.0)*b1_vy + (-8.0)*b2_vy +
			(3680.0/513.0)*b3_vy + (-845.0/4104.0)*b4_vy);
    
    double a5_vx = dt * rkf.compute_x_accel(x1 + (439.0/216.0)*a1_x +
					    (-8.0)*a2_x +
					    (3680.0/513.0)*a3_x +
					    (-845.0/4104.0)*a4_x,
					    y1 + (439.0/216.0)*a1_y +
					    (-8.0)*a2_y +
					    (3680.0/513.0)*a3_y +
					    (-845.0/4104.0)*a4_y,
					    
					    x2 + (439.0/216.0)*b1_x +
					    (-8.0)*b2_x +
					    (3680.0/513.0)*b3_x +
					    (-845.0/4104.0)*b4_x,
					    y2 + (439.0/216.0)*b1_y +
					    (-8.0)*b2_y +
					    (3680.0/513.0)*b3_y +
					    (-845.0/4104.0)*b4_y,
					    center_mass_earth);
    
    double a5_vy = dt * rkf.compute_y_accel(x1 + (439.0/216.0)*a1_x +
					    (-8.0)*a2_x +
					    (3680.0/513.0)*a3_x +
					    (-845.0/4104.0)*a4_x,
					    y1 + (439.0/216.0)*a1_y +
					    (-8.0)*a2_y +
					    (3680.0/513.0)*a3_y +
					    (-845.0/4104.0)*a4_y,
					    
					    x2 + (439.0/216.0)*b1_x +
					    (-8.0)*b2_x +
					    (3680.0/513.0)*b3_x +
					    (-845.0/4104.0)*b4_x,
					    y2 + (439.0/216.0)*b1_y +
					    (-8.0)*b2_y +
					    (3680.0/513.0)*b3_y +
					    (-845.0/4104.0)*b4_y,
					    center_mass_earth);
    
    double b5_vx = dt * rkf.compute_x_accel(x2 + (439.0/216.0)*b1_x +
					    (-8.0)*b2_x +
					    (3680.0/513.0)*b3_x +
					    (-845.0/4104.0)*b4_x,
					    y2 + (439.0/216.0)*b1_y +
					    (-8.0)*b2_y +
					    (3680.0/513.0)*b3_y +
					    (-845.0/4104.0)*b4_y,
					    
					    x1 + (439.0/216.0)*a1_x +
					    (-8.0)*a2_x +
					    (3680.0/513.0)*a3_x +
					    (-845.0/4104.0)*a4_x,
					    y1 + (439.0/216.0)*a1_y +
					    (-8.0)*a2_y +
					    (3680.0/513.0)*a3_y +
					    (-845.0/4104.0)*a4_y,
					    center_mass_sun);
    
    double b5_vy = dt * rkf.compute_y_accel(x2 + (439.0/216.0)*b1_x +
					    (-8.0)*b2_x +
					    (3680.0/513.0)*b3_x +
					    (-845.0/4104.0)*b4_x,
					    y2 + (439.0/216.0)*b1_y +
					    (-8.0)*b2_y +
					    (3680.0/513.0)*b3_y +
					    (-845.0/4104.0)*b4_y,
					    
					    x1 + (439.0/216.0)*a1_x +
					    (-8.0)*a2_x +
					    (3680.0/513.0)*a3_x +
					    (-845.0/4104.0)*a4_x,
					    y1 + (439.0/216.0)*a1_y +
					    (-8.0)*a2_y +
					    (3680.0/513.0)*a3_y +
					    (-845.0/4104.0)*a4_y,
					    center_mass_sun);
    
    
    
    //////////////////
    
    double a6_x = dt * (vx1 + (-8.0/27.0)*a1_vx + (2)*a2_vx +
			(-3544.0/2565.0)*a3_vx + (1859.0/4104.0)*a4_vx +
			(-11.0/40.0)*a5_vx);
    double a6_y = dt * (vy1 + (-8.0/27.0)*a1_vy + (2)*a2_vy +
			(-3544.0/2565.0)*a3_vy + (1859.0/4104.0)*a4_vy +
			(-11.0/40.0)*a5_vy);
    
    double b6_x = dt * (vx2 + (-8.0/27.0)*b1_vx + (2)*b2_vx +
			(-3544.0/2565.0)*b3_vx + (1859.0/4104.0)*b4_vx +
			(-11.0/40.0)*b5_vx);
    double b6_y = dt * (vy2 + (-8.0/27.0)*b1_vy + (2)*b2_vy +
			(-3544.0/2565.0)*b3_vy + (1859.0/4104.0)*b4_vy +
			(-11.0/40.0)*b5_vy);

    double a6_vx = dt * rkf.compute_x_accel(x1 + (-8.0/27.0)*a1_x +
					    (2)*a2_x +
					    (-3544.0/2565.0)*a3_x +
					    (1859.0/4104.0)*a4_x +
					    (-11.0/40.0)*a5_x,
					    y1 + (439.0/216.0)*a1_y +
					    (-8.0)*a2_y +
					    (3680.0/513.0)*a3_y +
					    (-845.0/4104.0)*a4_y,
					    
					    x2 + (-8.0/27.0)*b1_x +
					    (2)*b2_x +
					    (-3544.0/2565.0)*b3_x +
					    (1859.0/4104.0)*b4_x +
					    (-11.0/40.0)*b5_x,
					    y2 + (439.0/216.0)*b1_y +
					    (-8.0)*b2_y +
					    (3680.0/513.0)*b3_y +
					    (-845.0/4104.0)*b4_y +
					    (-11.0/40.0)*b5_y,
					    center_mass_earth);
    double a6_vy = dt * rkf.compute_y_accel(x1 + (-8.0/27.0)*a1_x +
					    (2)*a2_x +
					    (-3544.0/2565.0)*a3_x +
					    (1859.0/4104.0)*a4_x +
					    (-11.0/40.0)*a5_x,
					    y1 + (439.0/216.0)*a1_y +
					    (-8.0)*a2_y +
					    (3680.0/513.0)*a3_y +
					    (-845.0/4104.0)*a4_y,
					    
					    x2 + (-8.0/27.0)*b1_x +
					    (2)*b2_x +
					    (-3544.0/2565.0)*b3_x +
					    (1859.0/4104.0)*b4_x +
					    (-11.0/40.0)*b5_x,
					    y2 + (439.0/216.0)*b1_y +
					    (-8.0)*b2_y +
					    (3680.0/513.0)*b3_y +
					    (-845.0/4104.0)*b4_y +
					    (-11.0/40.0)*b5_y,
					    center_mass_earth);
    
    double b6_vx = dt * rkf.compute_x_accel(x2 + (-8.0/27.0)*b1_x +
					    (2)*b2_x +
					    (-3544.0/2565.0)*b3_x +
					    (1859.0/4104.0)*b4_x +
					    (-11.0/40.0)*b5_x,
					    y2 + (439.0/216.0)*b1_y +
					    (-8.0)*b2_y +
					    (3680.0/513.0)*b3_y +
					    (-845.0/4104.0)*b4_y +
					    (-11.0/40.0)*b5_y,
					    
					    x1 + (-8.0/27.0)*a1_x +
					    (2)*a2_x +
					    (-3544.0/2565.0)*a3_x +
					    (1859.0/4104.0)*a4_x +
					    (-11.0/40.0)*a5_x,
					    y1 + (439.0/216.0)*a1_y +
					    (-8.0)*a2_y +
					    (3680.0/513.0)*a3_y +
					    (-845.0/4104.0)*a4_y,
					    center_mass_sun);
    double b6_vy = dt * rkf.compute_y_accel(x2 + (-8.0/27.0)*b1_x +
					    (2)*b2_x +
					    (-3544.0/2565.0)*b3_x +
					    (1859.0/4104.0)*b4_x +
					    (-11.0/40.0)*b5_x,
					    y2 + (439.0/216.0)*b1_y +
					    (-8.0)*b2_y +
					    (3680.0/513.0)*b3_y +
					    (-845.0/4104.0)*b4_y +
					    (-11.0/40.0)*b5_y,
					    x1 + (-8.0/27.0)*a1_x +
					    (2)*a2_x +
					    (-3544.0/2565.0)*a3_x +
					    (1859.0/4104.0)*a4_x +
					    (-11.0/40.0)*a5_x,
					    y1 + (439.0/216.0)*a1_y +
					    (-8.0)*a2_y +
					    (3680.0/513.0)*a3_y +
					    (-845.0/4104.0)*a4_y+
					    (-11.0/40.0)*a5_y,
					    center_mass_sun);
    
    
    //////////////////
    
    
    x1 += (16.0/135.0)*a1_x + (0)*a2_x + (6656.0/12825.0)*a3_x + (28561.0/56430.0)*a4_x +
      (-9.0/50.0)*a5_x + (2.0/55.0)*a6_x;
    y1 += (16.0/135.0)*a1_y + (0)*a2_y + (6656.0/12825.0)*a3_y + (28561.0/56430.0)*a4_y +
      (-9.0/50.0)*a5_y + (2.0/55.0)*a6_y;
    
    x2 += (16.0/135.0)*b1_x + (0)*b2_x + (6656.0/12825.0)*b3_x + (28561.0/56430.0)*b4_x +
      (-9.0/50.0)*b5_x + (2.0/55.0)*b6_x;
    y2 += (16.0/135.0)*b1_y + (0)*b2_y + (6656.0/12825.0)*b3_y + (28561.0/56430.0)*b4_y +
      (-9.0/50.0)*b5_y + (2.0/55.0)*b6_y;
    
    vx1 += (16.0/135.0)*a1_vx + (0)*a2_vx + (6656.0/12825.0)*a3_vx + (28561.0/56430.0)*a4_vx +
      (-9.0/50.0)*a5_vx + (2.0/55.0)*a6_vx;
    vy1 += (16.0/135.0)*a1_vy + (0)*a2_vy + (6656.0/12825.0)*a3_vy + (28561.0/56430.0)*a4_vy +
      (-9.0/50.0)*a5_vy + (2.0/55.0)*a6_vy;
    
    vx2 += (16.0/135.0)*b1_vx + (0)*b2_vx + (6656.0/12825.0)*b3_vx + (28561.0/56430.0)*b4_vx +
      (-9.0/50.0)*b5_vx + (2.0/55.0)*b6_vx;
    vy2 += (16.0/135.0)*b1_vy + (0)*b2_vy + (6656.0/12825.0)*b3_vy + (28561.0/56430.0)*b4_vy +
      (-9.0/50.0)*b5_vy + (2.0/55.0)*b6_vy;
    
    
    double truncation_error = abs((1.0/360.0)*(a1_x+a1_y+b1_x+b1_y + a1_vx+a1_vy+b1_vx+b1_vy) +
				  
				  (-128.0/4275.0)*(a3_x+a3_y+b3_x+b3_y + a3_vx+a3_vy+b3_vx+b3_vy) +
				  
				  (-2197.0/75240.0)*(a4_x+a4_y+b4_x+b4_y + a4_vx+a4_vy+b4_vx+b4_vy) +
				  
				  (1.0/50.0)*(a5_x+a5_y+b5_x+b5_y + a5_vx+a5_vy+b5_vx+b5_vy) +
				  
				  (2.0/55.0)*( a6_x+a6_y+b6_x+b6_y + a6_vx+a6_vy+b6_vx+b6_vy));
    
    double epsilon = 0.2;
    
    double new_step_size = 0.9 * dt * pow((epsilon/(truncation_error+0.001)),(1.0/5.0));
    std::cout << "STEPPERINO: " << new_step_size << "\n";
    std::cout << "ERRORINO: " << truncation_error << "\n";
    std::cout << "P1: " << x1 << ", " << y1 << ". ";
    std::cout << "P2: " << x2 << ", " << y2 << "\n";
    if(truncation_error > epsilon)
      {
	dt = new_step_size;
	std::cout << "NEW STEP!\n" ;
      }
    else if(truncation_error <= epsilon)
      {
	eta += dt;
	dt = new_step_size;
      }

    //positions_output << x2 << "," << y2 << ",";
    //positions_output << x1 << "," << y1 << "\n";
  }
  //positions_output.close();
}
