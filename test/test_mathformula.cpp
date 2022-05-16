#define CATCH_CONFIG_MAIN
#include "catch.hpp"


#include "../src/ray_tracing/math_formula/math_formula.hpp"

#include <cmath>
#include <iostream>


TEST_CASE("math formula tests", "[unit]"){
  MathFormula MF;

  SECTION("omega_bar"){
    // Radius=0, a=1, theta=0 -> omega_bar ~= 0
    REQUIRE(MF.omega_bar(0,1,0) == Approx(0));

    double radius = 123123;
    double theta = 0;
    double a = 1;
    REQUIRE(MF.omega_bar(radius,a,theta) == Approx(0));
  }

  SECTION("sigma"){
    double radius = 0, theta = 123, a=1;
    REQUIRE(MF.omega(radius,a,theta) == Approx(0));

    double radius1 = 100, theta1=0, a1=1;
    REQUIRE(MF.omega(radius1,a1,theta1) ==
	    Approx((2*a1*radius1)/(pow(pow(radius1,2)+pow(a1,2),2))).epsilon(1e-5));
  }

  SECTION("capital omega"){
    double radius2=0, theta2=0, a2 = 1.0/5.0;
    REQUIRE(MF.Omega(radius2, a2, theta2) == Approx(1/a2));
    std::cout << "VALUE: " << 1/a2 <<"\n";
  }
}

TEST_CASE("b_zero and q_zero tests", "[unit]"){
  MathFormula MF;

  SECTION("b_zero"){
    double result = MF.b_zero(0,5);
    REQUIRE(result == Approx(5));
    
  }

}
