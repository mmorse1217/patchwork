#include "bspline_codegen.hpp"


double bspline_deg3_deriv0(double x) {
   double bspline_deg3_deriv0_result;
   if (x >= -3 && x < -2) {
      bspline_deg3_deriv0_result = x*(x*((1.0L/6.0L)*x + 3.0L/2.0L) + 9.0L/2.0L) + 9.0L/2.0L;
   }
   else if (x >= -2 && x < -1) {
      bspline_deg3_deriv0_result = x*(x*(-1.0L/2.0L*x - 5.0L/2.0L) - 7.0L/2.0L) - 5.0L/6.0L;
   }
   else if (x >= -1 && x < 0) {
      bspline_deg3_deriv0_result = x*(x*((1.0L/2.0L)*x + 1.0L/2.0L) - 1.0L/2.0L) + 1.0L/6.0L;
   }
   else if (x >= 0 && x <= 1) {
      bspline_deg3_deriv0_result = x*(x*(-1.0L/6.0L*x + 1.0L/2.0L) - 1.0L/2.0L) + 1.0L/6.0L;
   }
   else {
      bspline_deg3_deriv0_result = 0;
   }
   return bspline_deg3_deriv0_result;
}


double bspline_deg3_deriv1(double x) {
   double bspline_deg3_deriv1_result;
   if (x >= -3 && x < -2) {
      bspline_deg3_deriv1_result = x*((1.0L/2.0L)*x + 3) + 9.0L/2.0L;
   }
   else if (x >= -2 && x < -1) {
      bspline_deg3_deriv1_result = x*(-3.0L/2.0L*x - 5) - 7.0L/2.0L;
   }
   else if (x >= -1 && x < 0) {
      bspline_deg3_deriv1_result = x*((3.0L/2.0L)*x + 1) - 1.0L/2.0L;
   }
   else if (x >= 0 && x <= 1) {
      bspline_deg3_deriv1_result = x*(-1.0L/2.0L*x + 1) - 1.0L/2.0L;
   }
   else {
      bspline_deg3_deriv1_result = 0;
   }
   return bspline_deg3_deriv1_result;
}


double bspline_deg3_deriv2(double x) {
   double bspline_deg3_deriv2_result;
   if (x >= -3 && x < -2) {
      bspline_deg3_deriv2_result = x + 3;
   }
   else if (x >= -2 && x < -1) {
      bspline_deg3_deriv2_result = -3*x - 5;
   }
   else if (x >= -1 && x < 0) {
      bspline_deg3_deriv2_result = 3*x + 1;
   }
   else if (x >= 0 && x <= 1) {
      bspline_deg3_deriv2_result = -x + 1;
   }
   else {
      bspline_deg3_deriv2_result = 0;
   }
   return bspline_deg3_deriv2_result;
}


double bspline_deg4_deriv0(double x) {
   double bspline_deg4_deriv0_result;
   if (x >= -4 && x < -3) {
      bspline_deg4_deriv0_result = x*(x*(x*((1.0L/24.0L)*x + 2.0L/3.0L) + 4) + 32.0L/3.0L) + 32.0L/3.0L;
   }
   else if (x >= -3 && x < -2) {
      bspline_deg4_deriv0_result = x*(x*(x*(-1.0L/6.0L*x - 11.0L/6.0L) - 29.0L/4.0L) - 71.0L/6.0L) - 149.0L/24.0L;
   }
   else if (x >= -2 && x < -1) {
      bspline_deg4_deriv0_result = x*(x*(x*((1.0L/4.0L)*x + 3.0L/2.0L) + 11.0L/4.0L) + 3.0L/2.0L) + 11.0L/24.0L;
   }
   else if (x >= -1 && x < 0) {
      bspline_deg4_deriv0_result = x*(x*(x*(-1.0L/6.0L*x - 1.0L/6.0L) + 1.0L/4.0L) - 1.0L/6.0L) + 1.0L/24.0L;
   }
   else if (x >= 0 && x <= 1) {
      bspline_deg4_deriv0_result = x*(x*(x*((1.0L/24.0L)*x - 1.0L/6.0L) + 1.0L/4.0L) - 1.0L/6.0L) + 1.0L/24.0L;
   }
   else {
      bspline_deg4_deriv0_result = 0;
   }
   return bspline_deg4_deriv0_result;
}


double bspline_deg4_deriv1(double x) {
   double bspline_deg4_deriv1_result;
   if (x >= -4 && x < -3) {
      bspline_deg4_deriv1_result = x*(x*((1.0L/6.0L)*x + 2) + 8) + 32.0L/3.0L;
   }
   else if (x >= -3 && x < -2) {
      bspline_deg4_deriv1_result = x*(x*(-2.0L/3.0L*x - 11.0L/2.0L) - 29.0L/2.0L) - 71.0L/6.0L;
   }
   else if (x >= -2 && x < -1) {
      bspline_deg4_deriv1_result = x*(x*(x + 9.0L/2.0L) + 11.0L/2.0L) + 3.0L/2.0L;
   }
   else if (x >= -1 && x < 0) {
      bspline_deg4_deriv1_result = x*(x*(-2.0L/3.0L*x - 1.0L/2.0L) + 1.0L/2.0L) - 1.0L/6.0L;
   }
   else if (x >= 0 && x <= 1) {
      bspline_deg4_deriv1_result = x*(x*((1.0L/6.0L)*x - 1.0L/2.0L) + 1.0L/2.0L) - 1.0L/6.0L;
   }
   else {
      bspline_deg4_deriv1_result = 0;
   }
   return bspline_deg4_deriv1_result;
}


double bspline_deg4_deriv2(double x) {
   double bspline_deg4_deriv2_result;
   if (x >= -4 && x < -3) {
      bspline_deg4_deriv2_result = x*((1.0L/2.0L)*x + 4) + 8;
   }
   else if (x >= -3 && x < -2) {
      bspline_deg4_deriv2_result = x*(-2*x - 11) - 29.0L/2.0L;
   }
   else if (x >= -2 && x < -1) {
      bspline_deg4_deriv2_result = x*(3*x + 9) + 11.0L/2.0L;
   }
   else if (x >= -1 && x < 0) {
      bspline_deg4_deriv2_result = x*(-2*x - 1) + 1.0L/2.0L;
   }
   else if (x >= 0 && x <= 1) {
      bspline_deg4_deriv2_result = x*((1.0L/2.0L)*x - 1) + 1.0L/2.0L;
   }
   else {
      bspline_deg4_deriv2_result = 0;
   }
   return bspline_deg4_deriv2_result;
}


double bspline_deg5_deriv0(double x) {
   double bspline_deg5_deriv0_result;
   if (x >= -5 && x < -4) {
      bspline_deg5_deriv0_result = x*(x*(x*(x*((1.0L/120.0L)*x + 5.0L/24.0L) + 25.0L/12.0L) + 125.0L/12.0L) + 625.0L/24.0L) + 625.0L/24.0L;
   }
   else if (x >= -4 && x < -3) {
      bspline_deg5_deriv0_result = x*(x*(x*(x*(-1.0L/24.0L*x - 19.0L/24.0L) - 71.0L/12.0L) - 259.0L/12.0L) - 911.0L/24.0L) - 3019.0L/120.0L;
   }
   else if (x >= -3 && x < -2) {
      bspline_deg5_deriv0_result = x*(x*(x*(x*((1.0L/12.0L)*x + 13.0L/12.0L) + 16.0L/3.0L) + 73.0L/6.0L) + 38.0L/3.0L) + 313.0L/60.0L;
   }
   else if (x >= -2 && x < -1) {
      bspline_deg5_deriv0_result = x*(x*(x*(x*(-1.0L/12.0L*x - 7.0L/12.0L) - 4.0L/3.0L) - 7.0L/6.0L) - 2.0L/3.0L) - 7.0L/60.0L;
   }
   else if (x >= -1 && x < 0) {
      bspline_deg5_deriv0_result = x*(x*(x*(x*((1.0L/24.0L)*x + 1.0L/24.0L) - 1.0L/12.0L) + 1.0L/12.0L) - 1.0L/24.0L) + 1.0L/120.0L;
   }
   else if (x >= 0 && x <= 1) {
      bspline_deg5_deriv0_result = x*(x*(x*(x*(-1.0L/120.0L*x + 1.0L/24.0L) - 1.0L/12.0L) + 1.0L/12.0L) - 1.0L/24.0L) + 1.0L/120.0L;
   }
   else {
      bspline_deg5_deriv0_result = 0;
   }
   return bspline_deg5_deriv0_result;
}


double bspline_deg5_deriv1(double x) {
   double bspline_deg5_deriv1_result;
   if (x >= -5 && x < -4) {
      bspline_deg5_deriv1_result = x*(x*(x*((1.0L/24.0L)*x + 5.0L/6.0L) + 25.0L/4.0L) + 125.0L/6.0L) + 625.0L/24.0L;
   }
   else if (x >= -4 && x < -3) {
      bspline_deg5_deriv1_result = x*(x*(x*(-5.0L/24.0L*x - 19.0L/6.0L) - 71.0L/4.0L) - 259.0L/6.0L) - 911.0L/24.0L;
   }
   else if (x >= -3 && x < -2) {
      bspline_deg5_deriv1_result = x*(x*(x*((5.0L/12.0L)*x + 13.0L/3.0L) + 16) + 73.0L/3.0L) + 38.0L/3.0L;
   }
   else if (x >= -2 && x < -1) {
      bspline_deg5_deriv1_result = x*(x*(x*(-5.0L/12.0L*x - 7.0L/3.0L) - 4) - 7.0L/3.0L) - 2.0L/3.0L;
   }
   else if (x >= -1 && x < 0) {
      bspline_deg5_deriv1_result = x*(x*(x*((5.0L/24.0L)*x + 1.0L/6.0L) - 1.0L/4.0L) + 1.0L/6.0L) - 1.0L/24.0L;
   }
   else if (x >= 0 && x <= 1) {
      bspline_deg5_deriv1_result = x*(x*(x*(-1.0L/24.0L*x + 1.0L/6.0L) - 1.0L/4.0L) + 1.0L/6.0L) - 1.0L/24.0L;
   }
   else {
      bspline_deg5_deriv1_result = 0;
   }
   return bspline_deg5_deriv1_result;
}


double bspline_deg5_deriv2(double x) {
   double bspline_deg5_deriv2_result;
   if (x >= -5 && x < -4) {
      bspline_deg5_deriv2_result = x*(x*((1.0L/6.0L)*x + 5.0L/2.0L) + 25.0L/2.0L) + 125.0L/6.0L;
   }
   else if (x >= -4 && x < -3) {
      bspline_deg5_deriv2_result = x*(x*(-5.0L/6.0L*x - 19.0L/2.0L) - 71.0L/2.0L) - 259.0L/6.0L;
   }
   else if (x >= -3 && x < -2) {
      bspline_deg5_deriv2_result = x*(x*((5.0L/3.0L)*x + 13) + 32) + 73.0L/3.0L;
   }
   else if (x >= -2 && x < -1) {
      bspline_deg5_deriv2_result = x*(x*(-5.0L/3.0L*x - 7) - 8) - 7.0L/3.0L;
   }
   else if (x >= -1 && x < 0) {
      bspline_deg5_deriv2_result = x*(x*((5.0L/6.0L)*x + 1.0L/2.0L) - 1.0L/2.0L) + 1.0L/6.0L;
   }
   else if (x >= 0 && x <= 1) {
      bspline_deg5_deriv2_result = x*(x*(-1.0L/6.0L*x + 1.0L/2.0L) - 1.0L/2.0L) + 1.0L/6.0L;
   }
   else {
      bspline_deg5_deriv2_result = 0;
   }
   return bspline_deg5_deriv2_result;
}
