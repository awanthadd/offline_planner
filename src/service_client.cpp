//////////////////////////////////////////////////////////////////////////////////////////////////////////
//   Project      : OFFLINE PATH PLANNER FOR BELL 412 HELICOPTOR                                        //
//   File         : Service_client.cpp                                                                  //  
//   Description  : ROS service client for invoking service server to calculate the offline             //
//                  planning solution and generate reference points necessary for trajectory control.   //     
//   Created On   : 28/02/2021                                                                          //
//   Created By   : Awantha Jayasiri <mailto:awantha.jayasiri@nrc-cnrc.gc.ca>                           //
//////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <nlohmann/json.hpp>
#include "ros/ros.h"
#include "std_srvs/Empty.h"
#include <nav_msgs/Path.h>
//#include <nav_msgs/Odometry.h>
#include <geometry_msgs/PoseStamped.h>
//#include <tf/transform_listener.h>
#include "offline_planner/Mission_Service.h"
#include <offline_planner/display.hpp>  //For printing the data from service end for testing
#include <cstdlib>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>
#include <cmath>
#include <string>
#include <std_msgs/Float64.h>
#include <offline_planner/my_structs.hpp>
#include <offline_planner/Mission_struct.h>
#include <offline_planner/Mission_item_struct.h>
#include <offline_planner/Mission_struct.h>

//#include <offline_planner/my_functions.hpp>
#include <offline_planner/my_functions_3D.hpp>  // 18/05/2021 To include 3D waypoints with climb angle +ve
////////////// For MAVROS
#include <mavros_msgs/WaypointPush.h>
#include <mavros_msgs/WaypointClear.h>
#include <mavros_msgs/CommandHome.h>
#include <mavros_msgs/Waypoint.h>

////////////////////
std::vector <geometry_msgs::Point> corridor(const geometry_msgs::Point& pt1,const geometry_msgs::Point& pt2, const double& lm, const double& rm, const double& um, const double& dm)
{ std::vector <geometry_msgs::Point> tunnel_vec;
geometry_msgs::Point p1; geometry_msgs::Point p2; geometry_msgs::Point p3; geometry_msgs::Point p4;
geometry_msgs::Point p5; geometry_msgs::Point p6; geometry_msgs::Point p7; geometry_msgs::Point p8;

 p1.x=pt1.x; p1.y=pt1.y; p1.z=pt1.z-dm; p2.x=pt2.x; p2.y=pt2.y; p2.z=pt2.z-dm; p3=p2, p4=p1;   
 p5.x=pt1.x; p5.y=pt1.y; p5.z=pt1.z+um; p6.x=pt2.x; p6.y=pt2.y; p6.z=pt2.z+um; p7=p6, p8=p5;   


  double dx=pt2.x-pt1.x;  double dy=pt2.y-pt1.y;
  double normv=sqrt((dx*dx) + (dy*dy)) ; 
  double left_vect_x=-(dy/normv)*lm; double left_vect_y=(dx/normv)*lm;
  double right_vect_x=(dy/normv)*rm; double right_vect_y=-(dx/normv)*rm;
  p1.x+=left_vect_x; p1.y+=left_vect_y; 
  p2.x+=left_vect_x; p2.y+=left_vect_y;
  p3.x+=right_vect_x; p3.y+=right_vect_y; 
  p4.x+=right_vect_x; p4.y+=right_vect_y;

  p5.x+=left_vect_x; p5.y+=left_vect_y; 
  p6.x+=left_vect_x; p6.y+=left_vect_y;
  p7.x+=right_vect_x; p7.y+=right_vect_y; 
  p8.x+=right_vect_x; p8.y+=right_vect_y;

  if (NUM_SIDES_FOR_CORRIDOR_DISP==1){tunnel_vec.push_back(p1); tunnel_vec.push_back(p2); }
  else if (NUM_SIDES_FOR_CORRIDOR_DISP==2) {tunnel_vec.push_back(p1); tunnel_vec.push_back(p2); tunnel_vec.push_back(p4); tunnel_vec.push_back(p3);}  
  else if (NUM_SIDES_FOR_CORRIDOR_DISP==3) {tunnel_vec.push_back(p1); tunnel_vec.push_back(p2); tunnel_vec.push_back(p4); tunnel_vec.push_back(p3); tunnel_vec.push_back(p5); tunnel_vec.push_back(p6);}  
  else if (NUM_SIDES_FOR_CORRIDOR_DISP==4) {tunnel_vec.push_back(p1); tunnel_vec.push_back(p2); tunnel_vec.push_back(p4); tunnel_vec.push_back(p3); tunnel_vec.push_back(p5); tunnel_vec.push_back(p6); tunnel_vec.push_back(p8); tunnel_vec.push_back(p7);}  
  else { } 
  return tunnel_vec;
}
std::vector <geometry_msgs::Point> corridor_corrected(const std::vector <geometry_msgs::Point> corridor_pts){
 std::vector <geometry_msgs::Point> new_corridor_pts;
 std::vector <geometry_msgs::Point> new_corridor_pts2;

//int num_sides=4;
 for (int j=0;j<NUM_SIDES_FOR_CORRIDOR_DISP;j++){
  bool check=false; geometry_msgs::Point px;
  for (int i=j*2; i<corridor_pts.size()-(NUM_SIDES_FOR_CORRIDOR_DISP*2+1);i+=NUM_SIDES_FOR_CORRIDOR_DISP*2){
geometry_msgs::Point p;
geometry_msgs::Point p1=corridor_pts.at(i); geometry_msgs::Point p2=corridor_pts.at(i+1);
geometry_msgs::Point p3=corridor_pts.at(i+NUM_SIDES_FOR_CORRIDOR_DISP*2); geometry_msgs::Point p4=corridor_pts.at(i+(NUM_SIDES_FOR_CORRIDOR_DISP*2+1)); px=p4;
double m1=(p2.y-p1.y)/(p2.x-p1.x); double m2=(p4.y-p3.y)/(p4.x-p3.x);
if (m1==m2) { p=p2;
} 
else {
p.x=(p3.y-p1.y + m1*p1.x-m2*p3.x)/(m1-m2);
p.y= (p.x-p1.x)*m1 +p1.y ;
p.z=p3.z; }
if (check==false) {new_corridor_pts.push_back(corridor_pts.at(i)); check=true;} 
else new_corridor_pts.push_back( new_corridor_pts.back());
new_corridor_pts.push_back(p);
 //new_corridor_pts.push_back(p); 
//if (i==corridor_pts.size()-6) {
 // new_corridor_pts.push_back(new_corridor_pts.back());
 // new_corridor_pts.push_back(p4); //else new_corridor_pts.push_back( new_corridor_pts.end());
// new_corridor_pts.push_back(p4);
//}
}
new_corridor_pts.push_back(new_corridor_pts.back());
new_corridor_pts.push_back(px);
 }
//new_corridor_pts.push_back(new_corridor_pts.back());
//new_corridor_pts.push_back(corridor_pts.back());

return new_corridor_pts;
}


int main(int argc, char **argv)
{
///////////////////////////////////////////////////////////////
//std::ifstream ifs("/home/awantha/catkin_ws2/src/offline_planner/plan_files/Ottawa_Mission.plan");
//std::ifstream ifs("/home/awantha/catkin_ws2/src/offline_planner/plan_files/MSNOL6-2.plan");
//std::ifstream ifs("/home/awantha/catkin_ws2/src/offline_planner/plan_files/MsnOL3-10.plan");
//std::ifstream ifs("/home/awantha/catkin_ws2/src/offline_planner/plan_files/MSNOL6-2_FO_FB.plan"); // FILE WITH A FO WP. 27/01/2022
//std::ifstream ifs("/home/awantha/catkin_ws2/src/offline_planner/plan_files/2022/22-XX-XX-EY_CTR_Demo_Mission-Coarse_v1.plan"); // A NEW FILE . 05/02/2022
std::ifstream ifs("/home/awantha/catkin_ws2/src/offline_planner/plan_files/2022/Coarse5.plan"); // A NEW FILE . 05/02/2022
//std::ifstream ifs("/home/awantha/catkin_ws2/src/offline_planner/plan_files/2022/Tst1.plan"); // A NEW FILE . 05/02/2022



//std::ifstream ifs("/home/awantha/catkin_ws2/src/offline_planner/plan_files/MsnOL7-2.plan");
//std::ifstream ifs("/home/awantha/catkin_ws2/src/offline_planner/plan_files/MSNOL6-2tst.plan");
//std::ifstream ifs("/home/awantha/catkin_ws2/src/offline_planner/plan_files/MSNOL6-2_mod.plan");
//std::ifstream ifs("/home/awantha/catkin_ws2/src/offline_planner/plan_files/pgncc2-EY_CCT.plan");


//std::ifstream ifs("/home/awantha/catkin_ws2/src/offline_planner/plan_files/MsnOL9_7d_Acc_Dec_Rwy_25.plan");
//std::ifstream ifs("/home/awantha/catkin_ws2/src/offline_planner/plan_files/MsnOL9_1d_Rwy_25_CCT_3_deg_App.plan");
//std::ifstream ifs("/home/awantha/catkin_ws2/src/offline_planner/plan_files/MsnOL9_1d_Rwy_25_CCT_3_deg_App_test.plan");
//std::ifstream ifs("/home/awantha/catkin_ws2/src/offline_planner/plan_files/MsnOL9-6d - Acc-Dec Rwy 25.plan");
//std::ifstream ifs("/home/awantha/catkin_ws2/src/offline_planner/plan_files/MsnOL9_2d_Rwy_25_CCT_Norm_Lin_App.plan");
//std::ifstream ifs("/home/awantha/catkin_ws2/src/offline_planner/plan_files/MsnOL9_3d_Rwy_25_CCT_2_Step_App.plan");
//std::ifstream ifs("/home/awantha/catkin_ws2/src/offline_planner/plan_files/MsnOL9_4d_EY_CCT_2_Step_App.plan");
//std::ifstream ifs("/home/awantha/catkin_ws2/src/offline_planner/plan_files/MsnOL9_5e_Taxi_to_25HoldShort.plan");
//std::ifstream ifs("/home/awantha/catkin_ws2/src/offline_planner/plan_files/MsnOL9_6d_Acc_Dec_Rwy_25.plan");

//std::ifstream ifs("/home/awantha/catkin_ws2/src/offline_planner/plan_files/21_09_24_Coarse_Mission_for_Input_to_Offline_planner_Rwy_25_CCT.plan");



nlohmann::json plan_file=nlohmann::json::parse(ifs);

mission_struct current_mission=parse_plan_file(plan_file);

//std::vector <point_LatLonAltVel>PLAN_WPs_Old= Generate_plan_WPs_LatLonAlt(current_mission);
std::vector <point_xyzvdldrfo>XYZ_wps_old= Generate_plan_WPs_XYZ_from_LatLonAlt(Generate_plan_WPs_LatLonAlt(current_mission));
//std::vector<double> convert_mission_to_wp_vector(current_mission);

//std::vector<double> convert_mission_to_wp_vector(current_mission);


::offline_planner::Mission_struct mission;
mission.cruiseSpeed= current_mission.cruiseSpeed;
mission.firmwareType= current_mission.firmwareType;
mission.hoverSpeed=current_mission.hoverSpeed;
mission.vehicleType=current_mission.vehicleType;
mission.version=current_mission.version;
mission.plannedHomePosition.Alt =current_mission.plannedHomePosition.Alt;
mission.plannedHomePosition.Lat =current_mission.plannedHomePosition.Lat;
mission.plannedHomePosition.Lon =current_mission.plannedHomePosition.Lon;
mission.plannedHomePosition.Vel =current_mission.plannedHomePosition.Vel;

for (int i=0; i<current_mission.items.size();i++){
  ::offline_planner::Mission_item_struct mission_item;
  mission_item.AMSLAltAboveTerrain=current_mission.items.at(i).AMSLAltAboveTerrain;
  mission_item.Altitude=current_mission.items.at(i).Altitude;
  mission_item.AltitudeMode=current_mission.items.at(i).AltitudeMode;
  mission_item.autoContinue=current_mission.items.at(i).autoContinue;
  mission_item.command=current_mission.items.at(i).command;
  mission_item.doJumpId=current_mission.items.at(i).doJumpId; 
  mission_item.frame=current_mission.items.at(i).frame;
  mission_item.type=current_mission.items.at(i).type;
  for (int j=0; j<current_mission.items.at(i).params.size();j++){
      mission_item.params.push_back(current_mission.items.at(i).params.at(j)); }
  mission.items.push_back( mission_item);    
}

std::vector <geometry_msgs::Point> XYZ_wps;
for (int i=0; i < XYZ_wps_old.size();i++){
  geometry_msgs::Point p1; p1.x=XYZ_wps_old.at(i).x; p1.y=XYZ_wps_old.at(i).y; p1.z=XYZ_wps_old.at(i).z;
  XYZ_wps.push_back(p1);
  std::cout<<"wway points: "<< XYZ_wps_old.at(i).x<<", "<< XYZ_wps_old.at(i).y<<", "<< XYZ_wps_old.at(i).z<<std::endl; 
 }

  // Initialize the ROS client Node
  ros::init(argc, argv, "mission_service_client_node");  // client node name
  ros::NodeHandle nh;
  //ros::NodeHandle nh2;

  ros::Publisher pub = nh.advertise<nav_msgs::Path>("/path",1000);
  ros::Publisher pub_h = nh.advertise< std_msgs::Float64>("/heading",1000);
  ros::Publisher pub_c = nh.advertise< std_msgs::Float64>("/climb",1000);
  ros::Publisher pub_b = nh.advertise< std_msgs::Float64>("/bank",1000);
  ros::Publisher pub_hv = nh.advertise< std_msgs::Float64>("/hvelocity",1000);
  ros::Publisher pub_vv = nh.advertise< std_msgs::Float64>("/vvelocity",1000);
  ros::Publisher marker_pub = nh.advertise<visualization_msgs::Marker>( "/visualization_marker", 10 );
  ros::Publisher markerArray_pub = nh.advertise<visualization_msgs::MarkerArray>("/MarkerArray", 10);
  //ros::Rate loop_rate(100);
  ros::Rate loop_rate(10);

  ros::ServiceClient mission_client = nh.serviceClient<offline_planner::Mission_Service>("mission_service"); // create the client (mission_client)
  offline_planner::Mission_Service srv;
 // std::vector<double> wp_vector2; // way point vector

 //std::vector<double> wp_vector= convert_wps_to_wp_vector(XYZ_wps_old);
 srv.request.mission=mission;
  if (mission_client.call(srv))   // Calls the service "mission_service"
  {
    ROS_INFO("Mission service server is called and successful route is achieved!!!..");

    nav_msgs::Path path;
    path.header.frame_id="/map";
    //std::vector<geometry_msgs::PoseStamped> plan;
    for(int i=0; i<srv.response.path.path_xyz.size();i++){
     geometry_msgs:: PoseStamped pose;
     pose.header.seq=i;
     pose.header.frame_id="/map";
      pose.pose.position.x=srv.response.path.path_xyz[i].x;
      //pose.pose.position.x=XYZ_wps_old.at(i).x;
      pose.pose.position.y=srv.response.path.path_xyz[i].y;
      //pose.pose.position.y=XYZ_wps_old.at(i).y;
      pose.pose.position.z=srv.response.path.path_xyz[i].z;
      //pose.pose.position.z=XYZ_wps_old.at(i).z;
      path.poses.push_back(pose);        
     
    }
    


///////////////////////////////GENERATE NEW WAY-POINTS FILE...THIS COMES LATER, AFTER CALULATING THE OFFLINE PLANNER

mission_struct new_mission= process_mission_items(srv.response.new_mission);
nlohmann::json new_plan_file;  
// Make a json object using new_mission struct.
nlohmann::json newj=new_mission;
for (auto it =plan_file.begin(); it != plan_file.end(); ++it)
{
   if(it.key() == "mission") new_plan_file[it.key()]=newj;
   else  new_plan_file[it.key()]=it.value();  
}
// Generating the new plan (json) file.
 //Sstd::ofstream file("/home/awantha/catkin_ws2/src/offline_planner/plan_files/test.json"); 
 std::ofstream file("/home/awantha/catkin_ws2/src/offline_planner/plan_files/new_plan.plan"); 
 file << std::setw(4)<<new_plan_file<<std::endl;

///////////////////////////////////////////////////////////





    visualization_msgs::MarkerArray markerArray;
  for(int f=0;f<XYZ_wps.size();f++){
  visualization_msgs::Marker  text; text.header.frame_id = "/map"; text.action= visualization_msgs::Marker::ADD; text.id = f;
   text.type= visualization_msgs::Marker::TEXT_VIEW_FACING; text.scale.x = 40; text.scale.y = 40; text.scale.z = 40;
   text.color.r = 1.0; text.color.g = 1.0; text.color.b = 1.0; text.color.a = 1.0;
   std::string s1=std::to_string(XYZ_wps_old.at(f).x); std::string s2=std::to_string(XYZ_wps_old.at(f).y); std::string s3=std::to_string(XYZ_wps_old.at(f).z); 
   std::string s4=std::to_string(XYZ_wps_old.at(f).v);
   std::string sx="\n"; std::string txt=s1+sx+ s2+sx+s3+sx+s4;
  text.text=txt;
  text.pose.position.x = XYZ_wps_old.at(f).x;
  text.pose.position.y = XYZ_wps_old.at(f).y;
  text.pose.position.z = XYZ_wps_old.at(f).z;
markerArray.markers.push_back(text); 
}

    int cnt=0;
    while (ros::ok()){
        visualization_msgs::Marker points, line_strip, line_list, text;
        
        points.header.frame_id = line_strip.header.frame_id = line_list.header.frame_id = text.header.frame_id = "/map";
        points.header.stamp = line_strip.header.stamp = line_list.header.stamp = ros::Time::now();
        points.ns = line_strip.ns = line_list.ns = text.ns = "points_and_lines";
        points.action = line_strip.action = line_list.action =text.action= visualization_msgs::Marker::ADD;
        points.pose.orientation.w = line_strip.pose.orientation.w = line_list.pose.orientation.w = text.pose.orientation.w=1.0;
        points.id = 0;line_strip.id = 1; line_list.id = 2; text.id = 3;

        points.type = visualization_msgs::Marker::POINTS;
        line_strip.type = visualization_msgs::Marker::LINE_STRIP;
        line_list.type = visualization_msgs::Marker::LINE_LIST;
        text.type= visualization_msgs::Marker::TEXT_VIEW_FACING;

        // POINTS markers use x and y scale for width/height respectively
        points.scale.x = 40; points.scale.y = 40; // change scale 17/06/2021
        text.scale.x = 200; text.scale.y = 200; text.scale.z = 200;
        // LINE_STRIP/LINE_LIST markers use only the x component of scale, for the line width
        line_strip.scale.x = 5;  // change this to 10, 17/06/2021
        line_list.scale.x = 5;

        // Way Points are blue
        points.color.b = 10.0; points.color.a = 10.0;

        // Line list is red
        line_list.color.r = 1.0; line_list.color.a = 1.0;
        text.color.r = 1.0;  text.color.g = 1.0;  text.color.b = 1.0;   text.color.a = 1.0;

        // Line strip is -blue
        line_strip.color.b = 5.0;line_strip.color.a = 5.0;    

///////////////////////////////////////////// FOR DISPLAY
std::vector <geometry_msgs::Point> way_points;
std::vector <geometry_msgs::Point> corridor_points;
std::vector <geometry_msgs::Point> corridor_points2;
for (int i=0; i< XYZ_wps.size()-1; i++) {
geometry_msgs::Point p1; geometry_msgs::Point p2; double lm=0.0; double rm=0.0;

  // p1.x =srv.request.wp_data[i*(srv.request.wp_length)];
   p1.x=XYZ_wps.at(i).x;
   //p1.y =srv.request.wp_data[i*(srv.request.wp_length)+1]; 
   p1.y=XYZ_wps.at(i).y;
   //p1.z =srv.request.wp_data[i*(srv.request.wp_length)+2]; 
   p1.z=XYZ_wps.at(i).z;
   //lm= srv.request.wp_data[i*(srv.request.wp_length)+4];
   lm=DIST_L;
   //rm= srv.request.wp_data[i*(srv.request.wp_length)+5];
   rm=DIST_R;
   //p2.x =srv.request.wp_data[(i+1)*(srv.request.wp_length)];
   p2.x=XYZ_wps.at(i+1).x;
   //p2.y =srv.request.wp_data[(i+1)*(srv.request.wp_length)+1]; 
   p2.y=XYZ_wps.at(i+1).y;
   //p2.z =srv.request.wp_data[(i+1)*(srv.request.wp_length)+2]; 
   p2.z=XYZ_wps.at(i+1).z;

std::vector <geometry_msgs::Point> cps= corridor( p1, p2, DIST_L,DIST_R,DIST_U,DIST_B); 
 //std::cout<<" corridor1 size: "<< cps.size()<<std::endl; 



corridor_points.insert(std::end(corridor_points), std::begin(cps), std::end(cps));
way_points.push_back(p1);
if (i==XYZ_wps.size()-2) way_points.push_back(p2);
}

corridor_points2=corridor_corrected(corridor_points); //removed in 17/06/2021

corridor_points.clear(); //corridor_points2.clear();

std_msgs::Float64 h_msg; h_msg.data=srv.response.path.path_heading.at(cnt);
std_msgs::Float64 c_msg; c_msg.data=srv.response.path.path_climb.at(cnt);
std_msgs::Float64 b_msg; b_msg.data=srv.response.path.path_bank.at(cnt);
std_msgs::Float64 hv_msg; hv_msg.data=srv.response.path.path_velocity.at(cnt);
std_msgs::Float64 vv_msg; vv_msg.data=srv.response.path.path_vvelocity.at(cnt);

//std::cout<<"cnt "<<cnt<<" "<<"vel "<<srv.response.path.path_velocity.at(cnt)<<std::endl;
//std::cout<<"cnt "<<cnt<<" "<<"x,y "<< srv.response.path.path_xyz[cnt].x<<", "<< srv.response.path.path_xyz[cnt].y<<std::endl;
if (cnt == srv.response.path.path_heading.size()-1) break;
cnt++;
//points.points=way_points;
points.points=XYZ_wps;

line_list.points=corridor_points2;
corridor_points2.clear();
line_strip.points=XYZ_wps;
marker_pub.publish(points);
// PUBLISH FOLLOWING TO GET THE PATH
marker_pub.publish(line_strip);


//if (cnt<10) 
marker_pub.publish(text);
markerArray_pub.publish(markerArray);
marker_pub.publish(line_list);
pub.publish(path);
pub_h.publish(h_msg); pub_c.publish(c_msg); pub_b.publish(b_msg); pub_hv.publish(hv_msg); pub_vv.publish(vv_msg);
ros::spinOnce();
loop_rate.sleep();


    }
  }
  else
  {
    ROS_ERROR("Mission service server is called but route is unsuccessful. Please redefine the parameters ");
    return 1;
  }





  return 0;
}
