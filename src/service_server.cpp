//////////////////////////////////////////////////////////////////////////////////////////////////////////
//   Project      : OFFLINE PATH PLANNER FOR BELL 412 HELICOPTOR                                        //
//   File         : Service_server.cpp                                                                  //  
//   Description  : ROS service server (invoked by service client) to calculate the offline             //
//                  planning solution and generate reference points necessary for trajectory control.   //     
//   Created On   : 28/02/2021                                                                          //
//   Created By   : Awantha Jayasiri <mailto:awantha.jayasiri@nrc-cnrc.gc.ca>                           //
//////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "ros/ros.h"
#include "std_srvs/Empty.h"
#include "offline_planner/Mission_Service.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h> 
#include <utility> 
//#include <sstream>
#include <string>
#include <vector>
#include <offline_planner/my_structs.hpp>
#include <offline_planner/display.hpp>  ////For printing the data from service end for testing
//#include <offline_planner/my_functions.hpp>
#include <offline_planner/my_functions_3D.hpp>  // 18/05/2021 To include 3D waypoints with climb angle +ve


// We define the callback function of the service
bool mission_callback(offline_planner::Mission_Service::Request  &req, offline_planner::Mission_Service::Response &res)
 {
ROS_INFO("mission_callback has been called"); // We print an string whenever the Service gets called   

route_struct route;
mission_struct curr_mission = Read_mission_from_client_req(req);

std::vector <point_xyzvdldrfo>XYZ_curr_wps= Generate_plan_WPs_XYZ_from_LatLonAlt(Generate_plan_WPs_LatLonAlt(curr_mission));
route.wpt_list=XYZ_curr_wps;
route.no_of_wpts=XYZ_curr_wps.size();
route.num_turns=route.no_of_wpts-2;
route.num_straight_lines=route.no_of_wpts-1;

for (int i =0; i <XYZ_curr_wps.size(); i++)   // Filling the route with wp data sent by service client
  { 
      std::cout<<"wpt_val "<<XYZ_curr_wps.at(i).x<<", " <<XYZ_curr_wps.at(i).y<<", "<<XYZ_curr_wps.at(i).z<<", "<<XYZ_curr_wps.at(i).v<<", "<<XYZ_curr_wps.at(i).dl<<", "<<XYZ_curr_wps.at(i).dr<<", "<<XYZ_curr_wps.at(i).fowp<<", "<<std::endl;
  }
std::cin.get();
curv_struct curvs;
std::vector <std::vector <std::string> > curve_poly_data;
std::ifstream infile;
std::string line;
int no_of_velocities=0; int no_of_curves=0; double start_velocity=0.0; double end_velocity=0.0; double velocity_incr=0.0; 
//infile.open("/home/awantha/catkin_ws2/src/offline_planner/src/curv_polys_bank30_30_30.txt");// give the whole name of the file
//infile.open("/home/awantha/catkin_ws2/src/offline_planner/src/curv_polys_bank20_15_15.txt"); // 19/10/2021
infile.open("/home/awantha/catkin_ws2/src/offline_planner/src/curv_polys_bank20_10_10.txt"); // 21/10/2021


  
   while (infile.good())
  {    
    std::string s;
    if (!getline( infile, s )) break;
    std::istringstream ss( s );
     
    std::vector <std::string> record;
     while (ss)
    {
      std::string s;
      if (!getline( ss, s, ',' )) break;
      record.push_back( s );
    }

if (record.size()>0) curve_poly_data.push_back( record );
  }
  infile.close();

curvs.num_vels = std::stoi(curve_poly_data[0][0]); 
curvs.num_curv = std::stoi(curve_poly_data[0][1]);
curvs.min_vel = std::stof(curve_poly_data[0][2]); 
curvs.max_vel = std::stof(curve_poly_data[0][3]); 
curvs.vel_step = std::stof(curve_poly_data[0][4]);

std::vector< std::vector<std::string> >::iterator row_itr;
//std::vector<std::string>::iterator col_itr;

for (row_itr = curve_poly_data.begin()+1; row_itr != curve_poly_data.end(); row_itr++) { 
  poly_struct poly;
      poly.vel=std::stof(*(row_itr->begin()));     
      poly.ktrans=std::stof(*(row_itr->begin()+1));      
      poly.Sf1=std::stof(*(row_itr->begin()+2));      
      poly.Sf2=std::stof(*(row_itr->begin()+3));
      //std::cout<<poly.vel<<' '<<poly.ktrans<<' '<<poly.Sf1<<' '<<poly.Sf2;  
      poly.C1.push_back(std::stod(*(row_itr->begin()+4))); 
      poly.C1.push_back(std::stod(*(row_itr->begin()+5))); 
      poly.C1.push_back(std::stod(*(row_itr->begin()+6)));
      poly.C1.push_back(std::stod(*(row_itr->begin()+7)));
      poly.C1.push_back(std::stod(*(row_itr->begin()+8)));
      poly.C2.push_back(std::stod(*(row_itr->begin()+9))); 
      poly.C2.push_back(std::stod(*(row_itr->begin()+10))); 
      poly.C2.push_back(std::stod(*(row_itr->begin()+11)));
      poly.C2.push_back(std::stod(*(row_itr->begin()+12)));
      poly.C2.push_back(std::stod(*(row_itr->begin()+13)));

      curvs.curv_list.push_back(poly);       
}
  

bool turns_success=false; 
bool straights_success=false; 
bool straight_use_max_accel=true;  // use the maximum accelerator defined in heuristics or not

//calculate_ceres();
//std::cin.get();

// CALCULATE THE TURN SEGMENTS FIRST...
std::cout<<"Calculating turns..."<<std::endl;
for (int turn_i =0; turn_i <route.num_turns; turn_i++) {
  turn_struct new_turn;
  if (route.wpt_list.at(turn_i+1).fowp==true) {
    new_turn=generate_fly_over_turn( route, turn_i, curvs); 
  }
  else {
    new_turn=generate_fly_by_turn( route, turn_i, curvs); 
  }

  if (new_turn.empty==false) {
    route.turns.push_back(new_turn);
     std::cout<<"The turn velocity "<<turn_i<<" "<<new_turn.vel <<std::endl; 
    turns_success=true;
    std::cout<<"The turn "<<turn_i<<" is calculated successfully...Points:  "<<route.turns.at(turn_i).start.x<<","<<route.turns.at(turn_i).start.y <<","<<route.turns.at(turn_i).start.z <<",    "<<route.turns.at(turn_i).end.x<<","<<route.turns.at(turn_i).end.y <<","<<route.turns.at(turn_i).end.z <<std::endl;      
  } 
  else {
    turns_success=false;
    std::cout<<"Fail to calculate the turn "<<turn_i<< ". Please redefine the turn parameters. Existing.."<< std::endl;
    break;
  }
}

///////////////////////////////////////////////29/July/2021: THIS CODE SNIPPET IS JUST TO GENERATE THE WAYPOINT VECTOR FOR NEW PLAN FILE


// 13/09/2021 FOR TESGING THE SHIFT OF GENERATED PLAN FILE. TOGGLE BETWEEN NEW AND ORIGINAL XYZ WAYPOINTS
std::vector <point_xyzvdldrfo> new_plan_wps_XYZ= Generate_new_plan_wps_XYZ_from_turns(route.turns);
std::cin.get();
//std::vector <point_xyzvdldrfo> new_plan_wps_XYZ=XYZ_curr_wps;
// LOCATE WHERE THE ERROR IS......

std::vector <point_LatLonAltVelFo> new_plan_wps_LatLonAlt=Generate_plan_WPs_LatLonAlt_from_XYZ(new_plan_wps_XYZ, curr_mission);
std::cin.get();
//New mission that will be sent to the client. 
mission_struct updated_mission = Generate_mission_from_LatLonAlt(new_plan_wps_LatLonAlt,curr_mission);


std::cin.get();

//////////////////////////////////////////////////////////////////////////////////////////
if (turns_success==true){ // Continue only if all the turns are calculated accurately.
// NOW CALCULATE THE STRAIGHT LINE SEGMENTS...
std::cout<<"Calculating straight segments..."<<std::endl; 

for (int straight_i =0; straight_i <route.num_straight_lines; straight_i++) {

  straight_struct new_stright=generate_straight_line(route, straight_i, straight_use_max_accel);
  if (new_stright.empty==false) {   
    // Perform the calculation later as we do 3D wp now
    route.straights.push_back(new_stright);
    straights_success=true;
    std::cout<<"The straight "<<straight_i<<" is added successfully. Points:  "<<route.straights.at(straight_i).start.x<<","<<route.straights.at(straight_i).start.y<<","<<route.straights.at(straight_i).start.z <<",  "<<route.straights.at(straight_i).end.x<<","<<route.straights.at(straight_i).end.y <<","<<route.straights.at(straight_i).end.z <<std::endl;  
  } 
  else {
    straights_success=false;
    std::cout<<"Fail to calculate the straight "<<straight_i<< ". Please redefine the straight line parameters. Existing.."<< std::endl;   
    break;
  }
}

if (straights_success==true){
for (int i =0; i <route.straights.size(); i++) {
// Add new function to calculate straight coefficients...
 std::cout<<"The straight "<<i<<" is calculated successfully. Points:  "<<route.straights.at(i).start.x<<","<<route.straights.at(i).start.y <<",  "<<route.straights.at(i).end.x<<","<<route.straights.at(i).end.y <<std::endl; 

} }

}
if ( (turns_success==true) && (straights_success==true)){
// NOW CALCULATE THE PATH...
std::cout<<"Calculating straight and turn coefficients..."<<std::endl; 

  calculate_route_coefficients(route);
std::cout<<"Route total time "<< route.total_time<<std::endl; 

std::cout<<"Calculating complete path of the route..."<<std::endl; 

calculate_complete_path(route);

std::cout<<"Calculation is complete.."<<std::endl;

// Strart sending the path vector to client
std::cout<<"Sending the data to client process.."<<std::endl;
//res.path.path_heading.push_back(2000); //2
for (int i=0;i<route.path.path_xyz.size();i++){
  offline_planner::Point_XYZ cp;  // define the data type of the point as mentiond in srv file 
  cp.x=route.path.path_xyz.at(i).x;
  cp.y=route.path.path_xyz.at(i).y;
  cp.z=route.path.path_xyz.at(i).z;
  res.path.path_xyz.push_back(cp);
}
res.path.path_heading=route.path.path_heading;
res.path.path_velocity=route.path.path_velocity;
res.path.path_vvelocity=route.path.path_vvelocity;
res.path.path_bank=route.path.path_bank;
res.path.path_climb=route.path.path_climb;
res.path.path_curv=route.path.path_curv;
res.path.path_times=route.path.path_times;
res.path.path_deriv=route.path.path_deriv; // not filled yet
res.path.path_roll_rate=route.path.path_roll_rate;

// 13/ 09/2021  Change this to toggle between new and old missions
//res.new_mission= Send_mission_to_client_res(curr_mission);
res.new_mission= Send_mission_to_client_res(updated_mission);

return true;
}
else return false;
 }
  
int main(int argc, char **argv)
{
  // Initialize the ROS Node
  ros::init(argc, argv, "mission_service_server_node");  
  ros::NodeHandle nh;
  ros::ServiceServer my_service = nh.advertiseService("mission_service", mission_callback); // create the Service called "mission_service" 
  ROS_INFO("Ready to receive mission data waypoints");
  ros::spin(); // mantain the service open.​
  return 0;
}