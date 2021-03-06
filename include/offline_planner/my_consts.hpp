//////////////////////////////////////////////////////////////////////////////////////////////////////////
//   Project      : OFFLINE PATH PLANNER FOR BELL 412 HELICOPTOR                                        //
//   File         : my_consts.hpp                                                                       //  
//   Description  : Constant declarations used by the functions in my_fuctions.cpp.                     //     
//   Created On   : 28/02/2021                                                                          //
//   Created By   : Awantha Jayasiri <mailto:awantha.jayasiri@nrc-cnrc.gc.ca>                           //
//////////////////////////////////////////////////////////////////////////////////////////////////////////


#ifndef MY_CONSTS_HPP
#define MY_CONSTS_HPP

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif


// FOR CONVERTING LAT, LON DIFFERENCES TO DISTENCES
//const double EARTH_RADIUS = 6372797.56085;
const double EARTH_RADIUS = 6378100.00000;


// SOME CONSTANTS ARE READ THROUGH THE CURVES STRUCT FILLED BY THE CURVE_POLY_TXT

const double EPS= 1e-9 ;
const double WIND_SPEED= 5.0;  // in m/s...14/02/2022 for adding wind
const double WIND_DIR_FRM_E= M_PI/6;        //IN RADS, FROM EAST..... 14/02/2022 for adding wind

const double DIST_TOL= 1.0;//1e-3;// 1.0
//const double DIST_TOL_FOR_STRAIGHT_LINE= 0.001; // FOR DIFFERENCE IN STAIGHT LINE END TO CURVE START

const double GRAVITY= 9.81;
// **** IMPORTANT: IF YOU CHANGE THE FOLLOWING CONSTRAINTS PLEASE REGENERATE THE CURVE COEFFICIENT FILE FOR MAKE THESE CHANGES IN EFFECT****// 19/10/2021
const double max_airspeed=50.0; // according to the curve coefficients file...starts with an enven number ()
const double min_airspeed=10.0; //20.0; this is 10 for curv_polys2.txt file
const double max_roll=((20.0/180.0)* M_PI); // converting 30 degrees to radians
const double max_roll_rate=((15.0/180.0)* M_PI); 
const double max_roll_rate_rate=((15.0/180.0)* M_PI); 
//const double MAX_ACCEL=4.0;//2.0;//1.0;//0.1*GRAVITY; IN METERS PER SECOND SQ, will be calculated based on horisontal velocity and height above ground
const double max_jerk=0.1*GRAVITY;
//const double ASSUMED_ACCEL = 0.4 *MAX_ACCEL;
const double max_vel_z=1000*0.00508; // same as defined in KiTe algorithm
////////////*********************////////////////////


const int VV_COEFF_ORDER =4;
const int GAMMA_COEFF_ORDER =4;
const int CURV_POLY_DEG =4;
const int CURV_COEFF_LENGTH= CURV_POLY_DEG+1;
const int PSI_COEFF_LENGTH= CURV_COEFF_LENGTH+1;
const int NUM_CURVS_FOR_TURN= 3;
const int NUM_BREAKS_FOR_TURN= NUM_CURVS_FOR_TURN+1;
const int TURN_VEL_COEFF_LENGTH=2; // for now keep constant speed for entire turn for the simplecity and accelration is zero

const double HEAD_CHANGE_THRESH=0.05;//SHOULD BE CHANGED TO: 0.05; //0.018 (rad)=1.031 deg.// in rad (~3.00 deg) if the waypoints has a different more than this threshold it will calculate a fly-by turn
const int SPEED_POLY_DEG =4;
const int ACCEL_POLY_DEG= 3 ;
const int NUM_SAMPLES =100;
const double TIME_WEIGHT =1.0 ;
const double CURV_SMOOTH_WEIGHT =100.0;
const int NUM_SEGMENT_SAMPLES= 25;
const double TIME_STEP =0.1;//1.0

const double WPT_DIFF_THRESH=10.0;// difference in waypoints to calculate in m

// Default corridor distances for 2D path
const double DIST_L=50.0; // 50 // 300 in meters
const double DIST_R=50.0;

const double DIST_U=50.0; // in meters up and bottom distances
const double DIST_B=50.0;

const int NUM_SIDES_FOR_CORRIDOR_DISP=2; // NUMBER OF SIDES FOR DISPLAYING CORRIDOR IN CLIENT PROGRAM
//const double CLIMB_ANG_THRESH_RAD=0.017; // (1 deg) KEEPING A THRESHOLD FOR CLIMB ANGLE CALCULATION
const double VV_FM_HV_K_0_30 =300;//300;//200; //VERTICAL VELOCITY LIMIT (IN FEETS PER min) FOR HORIZONTAL VELOCITY (KNOTS) 0-30 =+- 200
const double VV_FM_HV_K_30_45 =800;//800;//500; //VERTICAL VELOCITY LIMIT (IN FEETS PER min) FOR HORIZONTAL VELOCITY (KNOTS) 30-45 =+- 500
const double VV_FM_HV_K_45 =1500; //VERTICAL VELOCITY LIMIT (IN FEETS PER min) FOR HORIZONTAL VELOCITY (KNOTS) 45+ =+- 1500

const double LA_KS_HAG_F_0_40_LV_K_0_45 =6; //2 //LONGITUDINAL ACCELERATION MAX IN k/S FOR HEIGHT ABOVE GROUND (<40 feet) AND LONG. VELOCITY (0-45 KNOTS)
const double LA_KS_HAG_F_40_I_LV_K_0_45 =10; //3 //LONGITUDINAL ACCELERATION MAX IN k/S FOR HEIGHT ABOVE GROUND (>40 feet) AND LONG. VELOCITY (0-45 KNOTS)
const double LA_KS_HAG_F_0_I_LV_K_45_I =12; //4 //LONGITUDINAL ACCELERATION MAX IN k/S FOR HEIGHT ABOVE GROUND (any feet) AND LONG. VELOCITY (>45 KNOTS)
const double KS_MSS=0.514444444 ; // meters oer second squred equvalency of 1 knots per second.


const double MS_TO_KNOTS=1.94384; // KNOTS FOR 1M/S
const double FM_TO_MS=0.00508; //   FEETS PER MIN TO METERS PER SEC

const double CLIMB_ANG_MAX_RAD_ABOVE_45K=0.279253; //  MAXIMUM CLIMB ANGLE FOR HORIZONTAL VELOCITIES ABOVE 45 KNOTS (23.15 M/S)
const double CLIMB_ANG_THRESH_RAD=0.01; //KEEPING A THRESHOLD FOR CLIMB ANGLE CALCULATION

const double ACCEPT_RADIUS=50.0; // acceptance radius defined for wps generated by offline planner

#endif
