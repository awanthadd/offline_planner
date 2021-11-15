//////////////////////////////////////////////////////////////////////////////////////////////////////////
//   Project      : OFFLINE PATH PLANNER FOR BELL 412 HELICOPTOR                                        //
//   File         : my_functions.hpp                                                                    //  
//   Description  : Functions declarations required to calculate the offline planning solution.         //     
//   Created On   : 28/02/2021                                                                          //
//   Created By   : Awantha Jayasiri <mailto:awantha.jayasiri@nrc-cnrc.gc.ca>                           //
//////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef MY_FUNCTIONS_3D_HPP
#define MY_FUNCTIONS_3D_HPP

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h> 
#include <utility> 
//#include <sstream>
#include <string>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/QR>

#include <gsl/gsl_integration.h> // GNU scintific library for numerical integration
#include <offline_planner/my_structs.hpp>
#include "offline_planner/Mission_Service.h"
#include <nlohmann/json.hpp>
#include <ceres/ceres.h>



double sign_fn(const double& temp);

bool dbl_to_bool(const double& dbl);
std::vector <point_xyz> generate_corridor(const point_xyzvdldrfo& pt1,const point_xyzvdldrfo& pt2);

struct poly_struct lookup_curve_splines(const bool& condition,const int& count, const double& vel, const curv_struct& curvs);

std::vector <double> polyint(const std::vector <double>& poly, const double& cons); // polynomial integration

double polyval(const std::vector <double>& poly, const double& x);  // polynomial evaluation

double ppval(const std::vector <std::vector <double> >& pp, const std::vector <double>& breaks, const double& x); // piesewise, polynomial, evaluation

struct initial_Sf2_out_struct get_initial_sf2(const poly_struct& C1C3Sf1Sf3k_trans,const float& head_change);

double f1 (double x, void *params);
double fcos0 (double x, void *params);
double fcos1 (double x, void *params);
double fcos2 (double x, void *params);
double fcos3 (double x, void *params);
double fsin1 (double x, void *params);
double fsin2 (double x, void *params);
double fsin3 (double x, void *params);

double time_fn (double x, void *params);

double integrater1(double (*f1) (double, void *params),const double& a, const double& b);
double integrater2(double (*f2) (double, void *params), std::vector <double> psi1, const double& a, const double& b);

double dis_frm_line_to_pt_2D(const point_xyz& pt, const point_xyz& pt1, const point_xyz& pt2);
double dis_frm_line_to_pt_3D(const point_xyz& pt, const point_xyz& pt1, const point_xyz& pt2);

point_xyz get_turn_end_point(const turn_struct& turn, const double& wind);

double dis_pt_to_pt_2D(const point_xyz& pt1, const point_xyz& pt2);

double dis_pt_to_pt_3D(const point_xyz& pt1, const point_xyz& pt2);

bool check_pt_bw_segments_2D(const point_xyz& p, const point_xyz& pt1, const point_xyz& pt2, const point_xyz& pt3);

bool check_pt_bw_segments_3D(const point_xyz& pt, const point_xyz& pt1, const point_xyz& pt2, const point_xyz& pt3);

find_curv_poly_fit_out_struct find_curv_poly_fit(const point_xyz& curv_end_pt, const point_xyz& x0,const point_xyz& curve_wpt2,const point_xyz& xlim, const double& dist_tol);

void setup_piecewice_poly_in_turn(turn_struct& turn); 

std::vector<double> vector_divide( const std::vector<double>& a, const std::vector<double>& b );

void get_path_from_turn(turn_struct& turn, const double& wind, const int& samples);

bool tunnel_check(const std::vector <point_xyz>& turn_xyz, const std::vector <point_xyz>& tnl1, const std::vector <point_xyz>& tnl2);

double set_max_airspeed(const double& mas, const double& step);

struct turn_struct generate_fly_by_turn(const route_struct& route, const int& turn_i, const curv_struct& curvs);

struct turn_struct generate_fly_by_turn_sub(const route_struct& route, const int& turn_i, const curv_struct& curvs,\
 const point_xyzvdldrfo& curve_wpt1, const point_xyzvdldrfo& curve_wpt2, const point_xyzvdldrfo& curve_wpt3,\
 const double& head_change, const double& heading1, const double& heading2);

struct turn_struct generate_climb_turn1(const route_struct& route, const int& i, const double& climb_angle, const double& heading);
struct turn_struct generate_climb_turn2(const route_struct& route, const int& i, const double& climb_angle, const double& heading, turn_struct& cl_turn);


struct turn_struct generate_fly_over_turn(const route_struct& route, const int& turn_i, const curv_struct& curvs);

struct straight_struct generate_straight_line(route_struct& route, const int& straight_i, const bool& use_max_accel);

void calculate_route_coefficients(route_struct& route);

void calculate_complete_path(route_struct& route);

void to_json(nlohmann::json& j, const mission_item_struct& mi);

void to_json(nlohmann::json& j, const mission_struct& m);

bool jexists(const nlohmann::json& j, const std::string& key);

struct mission_struct parse_plan_file(const nlohmann::json& jsonfile);

std::vector <point_LatLonAltVelFo> Generate_plan_WPs_LatLonAlt(const mission_struct& mission_in);

std::vector <point_xyzvdldrfo> Generate_plan_WPs_XYZ_from_LatLonAlt(const std::vector <point_LatLonAltVelFo>& plan_wps_LLA);
std::vector <point_LatLonAltVelFo> Generate_plan_WPs_LatLonAlt_from_XYZ(const std::vector <point_xyzvdldrfo>& plan_wps_XYZ, const mission_struct& mission);
//struct mission_struct process_mission_items(const mission_struct& mission_in);
struct mission_struct process_mission_items(const ::offline_planner::Mission_struct& mission_in);

std::vector <point_xyzvdldrfo> Generate_new_plan_wps_XYZ_from_turns(const std::vector <turn_struct>& turns);

std::vector<double> convert_mission_to_wp_vector(const mission_struct& mission_in);

std::vector<double> convert_wps_to_wp_vector(const std::vector <point_xyzvdldrfo>& wps);

double CalcGPSDistance(const double& lat1d, const double& long1d, const double& lat2d, const double& long2d); // Lat, Lon in degrees

double CalcGPSBearing(const double& lat1d, const double& long1d, const double& lat2d, const double& long2d);

double dxy_dgamma(const double& gamma, const double& v, const double& K, const double& Nz);

double dz_dgamma(const double& gamma, const double& v, const double& K, const double& Nz);

double dv_dgamma(const double& gamma,const double& v, const double& K, const double& Nz);

double dt_dgamma(const double& gamma, const double& v, const double& K, const double& Nz);

double dvh_dgamma(const double& gamma, const double& v, const double& K, const double& Nz);
double dvv_dgamma(const double& gamma, const double& v, const double& K, const double& Nz);


void RK4_out(RK4_struct& states, const double& K, const double& Nz,const double& dgamma, const int& i); // including rk4 and other onterations...
 
climb_start_end_struct find_climb_start_end(const point_xyz& climb_end, const point_xyz& x0,const point_xyz& curve_wpt2,const point_xyz& xlim, const double& dist_tol, const int& i);

void polyfit(const std::vector<double> &xv, const std::vector<double> &yv, std::vector<double> &coeff, int order);

void calculate_ceres ();



mission_struct Read_mission_from_client_req(const ::offline_planner::Mission_Service::Request& req);
::offline_planner::Mission_struct Send_mission_to_client_res(const mission_struct& mission);

mission_struct Generate_mission_from_LatLonAlt(const std::vector <point_LatLonAltVelFo>& wpsLatLonAlt, const mission_struct& mission);
double  return_max_long_accel(const double& hag, const double& lv ) ;


#endif