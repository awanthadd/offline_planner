//////////////////////////////////////////////////////////////////////////////////////////////////////////
//   Project      : OFFLINE PATH PLANNER FOR BELL 412 HELICOPTOR                                        //
//   File         : my_structs.hpp                                                                      //  
//   Description  : Structure definitions used by the functions.                                        //     
//   Created On   : 28/02/2021                                                                          //
//   Created By   : Awantha Jayasiri <mailto:awantha.jayasiri@nrc-cnrc.gc.ca>                           //
//////////////////////////////////////////////////////////////////////////////////////////////////////////


#ifndef MY_STRUCTS_HPP
#define MY_STRUCTS_HPP

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h> 
//#include <sstream>
#include <string>
#include <vector>
#include <functional> // std::divides 
#include <algorithm> // std::transform 
#include <numeric> // std::adjacent_difference 
#include <offline_planner/my_consts.hpp>


struct point_xyz  // a struct to hold (way) points
{
double x=0; double y=0; double z=0;
};
struct point_xyzvdldrfo  // a struct to hold (way) points
{
double x=0; double y=0; double z=0; double v=0; double dl=DIST_L; double dr=DIST_R; bool fowp=false; // (x,y,z,vel,x_length (default=DIST_L), y_length (default=DIST_R), flyover_wp?) 
};
struct point_LatLonAltVelFo  // a struct to hold (way) points in tat, long, alt, Vel to help parse mission items
{
double Lat=0; double Lon=0; double Alt=0; double Vel=0; bool fowp=false;
};
//Define the curves struct for holing all data
struct poly_struct{
double vel=0.0;
double ktrans=0.0;
double Sf1=0.0;
double Sf2=0.0;
std::vector <double> C1;
std::vector <double> C2;
};


struct path_struct{ // used in turn struct
std::vector <point_xyz> path_xyz;
std::vector <double> path_heading;
std::vector <double> path_velocity; // Horizontal_vel
std::vector <double> path_vvelocity; // vertical velocity
std::vector <double> path_bank;
std::vector <double> path_climb; // use only for vertical movement
std::vector <double> path_curv;
std::vector <double> path_times;
std::vector <double> path_deriv; // not used yet
std::vector <double> path_roll_rate;
};


struct RK4_struct
    {
      std::vector <double> xy; std::vector <double> z; std::vector <double> t;
      std::vector <double> x; std::vector <double> y; double heading =0.0; std::vector <double> vv;
      std::vector <double> v; std::vector <double> vh; std::vector <double> gamma;  std::vector <double> vv_coeffs; // vertical velocity coeffs   
    };

struct linear_climb_struct
    {
      std::vector <double> xy; std::vector <double> z; std::vector <double> t;
      std::vector <double> x; std::vector <double> y;   double heading =0.0;    std::vector <double> vv_coeffs; 
      double v=0.0; double vh=0.0; double vv=0.0; double gamma=0.0; double length=0.0;   double tot_t=0;
    };

struct climb_turn_struct {
  point_xyz x_start;
  point_xyz x_end;
  double straight_gamma=0.0; // climb angle for straight section
  int str_ind=0; // Straight index of the climb struct. this number is the one whit shows where this struct between the straights... 
  double heading=0.0; // Heading for the climb section for calcuation of x, y
  RK4_struct climb1;
  RK4_struct climb2;
  linear_climb_struct linear_climb;

  void set_heading(const double& head)
    {
      this->heading= head;
      this->climb1.heading=head;
      this->climb2.heading=head;
      this->linear_climb.heading=head;
    }

};

struct climb_start_end_struct{
point_xyz x_start;
point_xyz x_end;
bool matched=false;
};




struct turn_struct{
double vel=0.0; // this can  be a constant (for whole turn) or the result of distance integration of vel_poly polinomial..
bool isfowp=false; //default is a flyby waypoint 
bool isclimb=false; //default is a horizontal curve
bool is2Dcurve=false; //default is a vertical climb without a horizontal curve
bool is_real=true;
double ktrans=0.0;
//double Sf1=0.0;
//double Sf2=0.0;
//double Sf3=0.0;
double Sf_total=0.0;
bool empty =true;
double curv_heading1=0.0;
double curv_heading2=0.0;
std::vector <double> vel_poly; // Add a velocity polynomial if you want to implement a changing velocity turn in  HORIZONTAL

std::vector <double> Sfs;
std::vector <std::vector <double> > Cs;

point_xyz start; // turn start
point_xyz turn_point; // THE POINT OF THE TURN WHERE TURN HAPPEN (THE TURNING WAYPOINT) 
point_xyz end; // turn end
std::vector <point_xyz> points;

climb_turn_struct climb_turn; // This will only be used with 

std::vector <double> t_breaks;
std::vector <double> s_breaks;
std::vector <std::vector <double> > t_coeffs;
std::vector <std::vector <double> > curv_coeffs;
std::vector <std::vector <double> > psi_coeffs;
path_struct path;

double curve_max;    // From model struct 2D
double curve_rate_max;
double curve_rate_rate_max;


};



struct curv_struct
    {
      int num_vels=0; 
      int num_curv=0; 
      double min_vel=0.0; 
      double max_vel=0.0; 
      double vel_step=0.0;
      std::vector <poly_struct> curv_list;
    };

struct straight_struct
    {
      point_xyz start; // straight start
      point_xyz end; // straight end
      bool empty=true;
      bool isclimb=false; // default is a horizontal straight line
      bool use_max_accel=false;
      double vel_start=0.0;
      double vel_end=0.0;       
      double distance=0.0;
      double heading=0.0;
      double climb_angle=0.0;  // Only be used for climbs
      double time=0.0;

      double max_vel=0.0; 
      double accel=0.0;

      std::vector <double> t_breaks;
      std::vector <double> s_breaks;
      std::vector <std::vector <double> > t_coeffs;
      std::vector <std::vector <double> > vel_coeffs;
      std::vector <std::vector <double> > vvel_coeffs;
     // std::vector <point_xyz> path_xyz; Don't need this
      };



/*
struct model_struct   // initial values ar for Bell 205 from NEA
{
double max_airspeed=50.0; // starts with an enven number
double min_airspeed=20.0;
double max_roll=((15.0/180.0)* M_PI); // converting 15 degrees to radians
double max_roll_rate=((15.0/180.0)* M_PI); 
double max_roll_rate_rate=((15.0/180.0)* M_PI); 
double max_accel=0.1*GRAVITY;
double max_jerk=0.1*GRAVITY;
const double ASSUMED_ACCEL = 0.4 *max_accel;


double curve_max=GRAVITY*tan(max_roll)/(max_airspeed*max_airspeed);

double curve_rate_max=GRAVITY*(max_roll_rate)/(max_airspeed*max_airspeed*max_airspeed) - 2*curve_max*max_accel/(max_airspeed*max_airspeed);
double curve_rate_rate_max=GRAVITY*(max_roll_rate_rate)/(max_airspeed*max_airspeed*max_airspeed*max_airspeed);
double max_vel_z=1000*0.00508;

void set_max_airspeed(const double& mas)
    {
      this->max_airspeed= mas;
      if (static_cast<uint64_t>( floor(max_airspeed)) % 2 != 0 )
    {
       this->max_airspeed-=1.0;
    }
       
        this->curve_max=GRAVITY*tan(max_roll)/(max_airspeed*max_airspeed);
        this->curve_rate_max= GRAVITY*(max_roll_rate)/(max_airspeed*max_airspeed*max_airspeed) - 2*curve_max*max_accel/(max_airspeed*max_airspeed);
        this->curve_rate_rate_max= GRAVITY*(max_roll_rate_rate)/(max_airspeed*max_airspeed*max_airspeed*max_airspeed);
    }
};
*/

struct initial_Sf2_out_struct{
double rem_change=0.0;
double Sf2=0.0;
};

struct find_curv_poly_fit_out_struct{
point_xyz x_init;
point_xyz x_end;
bool matched=false;
};

struct route_struct {
  int no_of_wpts=0; 
  //int wp_len=0;   // data length 6 (x,y,z,vel,x_margin, y_margin)
  int num_turns=0;
  int num_straight_lines=0;
  std::vector <point_xyzvdldrfo> wpt_list;  

  std::vector <turn_struct> turns;
  std::vector <straight_struct> straights;
  std::vector <double> t_breaks;
  std::vector <double> s_breaks;
  std::vector <std::vector <double> > t_coeffs;
  std::vector <std::vector <double> > vel_coeffs; // horizontal vel
  std::vector <std::vector <double> > vvel_coeffs; // vertical vel
  std::vector <std::vector <double> > curv_coeffs;
  std::vector <std::vector <double> > psi_coeffs;
  std::vector <std::vector <double> > climb_ang_coeffs;
  
  double total_time=0;
  path_struct path;
};

////////////////////////////// FOR PARSING THE PLAN (JSON) FILES....
struct mission_item_struct{
double AMSLAltAboveTerrain=0;
double Altitude=0;
int AltitudeMode=0;
bool autoContinue=true;
int command=0;
int doJumpId=0;
int frame=0;
std::vector <double> params; // read mission item params in json to params vector, Check the command number for processing.
std::string type="SimpleItem";
};

struct mission_struct{
double cruiseSpeed=0;
int firmwareType=0;
double hoverSpeed=0;
std::vector <mission_item_struct> items; 
//std::vector <double> plannedHomePosition;
point_LatLonAltVelFo plannedHomePosition;
int vehicleType=0;
int version=0;
};


struct CostFunctor {
   template <typename T>
   bool operator()(const T* const x, T* residual) const {
     residual[0] = 10.0 - x[0];
     return true;
   }
};




#endif