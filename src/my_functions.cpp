//////////////////////////////////////////////////////////////////////////////////////////////////////////
//   Project      : OFFLINE PATH PLANNER FOR BELL 412 HELICOPTOR                                        //
//   File         : my_functions.cpp                                                                    //  
//   Description  : Functions definitions required to calculate the offline planning solution.          //     
//   Created On   : 28/02/2021                                                                          //
//   Created By   : Awantha Jayasiri <mailto:awantha.jayasiri@nrc-cnrc.gc.ca>                           //
//////////////////////////////////////////////////////////////////////////////////////////////////////////


#include <offline_planner/my_functions.hpp>
#include <offline_planner/display.hpp>

bool dbl_to_bool(const double& dbl){
return ((dbl== 1.0)? true : false);
}

std::vector <point_xyz> generate_corridor(const point_xyzvdldrfo& pt1,const point_xyzvdldrfo& pt2)
{ std::vector <point_xyz> tunnel_vec;
  point_xyz p1; point_xyz p2; point_xyz p3; point_xyz p4;
  p1.x=pt1.x; p1.y=pt1.y; p1.z=pt1.z; p2.x=pt2.x; p2.y=pt2.y; p2.z=pt2.z;
  p3=p2, p4=p1;
  double dx=pt2.x-pt1.x;
  double dy=pt2.y-pt1.y;
  double normv=sqrt((dx*dx) + (dy*dy)) ;
  double left_vect_x=-(dy/normv)*pt1.dl; 
  double left_vect_y=(dx/normv)*pt2.dl;
  double right_vect_x=(dy/normv)*pt1.dr; 
  double right_vect_y=-(dx/normv)*pt2.dr;

  p1.x+=left_vect_x; p1.y+=left_vect_y;
  p2.x+=left_vect_x; p2.y+=left_vect_y;
  p3.x+=right_vect_x; p3.y+=right_vect_y;  
  p4.x+=right_vect_x; p4.y+=right_vect_y;
  tunnel_vec.push_back(p1); tunnel_vec.push_back(p2); tunnel_vec.push_back(p3); tunnel_vec.push_back(p4);
  return tunnel_vec;
}

struct poly_struct lookup_curve_splines(const bool& condition,const int& count, const double& vel, const curv_struct& curvs){
poly_struct return_poly;

int vel_idx=int(((vel-curvs.min_vel)/curvs.vel_step)*curvs.num_curv);
int idx=0;
if (condition == 1) {
  idx=vel_idx+count;
}
else {
  idx=vel_idx+ (int)(curvs.num_curv/2) +count;
}
return_poly=curvs.curv_list.at(idx);
return return_poly;
}


std::vector <double> polyint(const std::vector <double>& poly, const double& cons){
int poly_length=poly.size();
std::vector <double> intg_coeff; 
for (int i=0; i<poly_length; i++){
  intg_coeff.push_back(poly.at(i)/(poly_length-i));  
}
intg_coeff.push_back(cons);
return intg_coeff;
}

double polyval(const std::vector <double>& poly, const double& x){
int poly_length=poly.size();
double val=0.0; 
for (int i=0; i<poly_length; i++){
val += poly.at(i)* pow(x,(double)(poly_length-1-i));  
}
return val;
}

double ppval(const std::vector <std::vector <double> >& pp, const std::vector <double>& breaks, const double& x) // piesewise, polynomial, evaluation
{
   int breaks_length=breaks.size();   double ppvalue=0.0;
  for (int i=0; i<breaks_length-1; i++){
    if ( (x>=breaks.at(i)) && (x<=breaks.at(i+1))  ){      
      ppvalue=polyval(pp.at(i),(x-breaks.at(i)));
      break;
    }else
    {
      continue;
    }    
  }
 return ppvalue;
}

struct initial_Sf2_out_struct get_initial_sf2(const poly_struct& C1C3Sf1Sf3k_trans,const float& head_change){
  initial_Sf2_out_struct st1;
  double Sf2=0.0;
  std::vector <double> psi1; 
  std::vector <double> psi3;
  psi1= polyint(C1C3Sf1Sf3k_trans.C1, 0.0);
  psi3= polyint(C1C3Sf1Sf3k_trans.C2, 0.0);
  double heading_change=polyval(psi1,C1C3Sf1Sf3k_trans.Sf1) + polyval(psi3,C1C3Sf1Sf3k_trans.Sf2);

  st1.rem_change=head_change-abs(heading_change);
    if (st1.rem_change<0){
      st1.Sf2=EPS;
    }
    else
    {
      st1.Sf2=abs(st1.rem_change/ C1C3Sf1Sf3k_trans.ktrans);
    }       
  return st1;
}
//////////////////////////////////////////////////////////////////
 double f1 (double x, void *params)
 {
 double n = *(double *) params;
 //double f = exp (-x) * pow (x, n);
 double f1 = log(n*x)/sqrt(x);
  return f1;
}
double integrater1(double (*f1) (double, void *params), const double& a, const double& b)
{
// a and b are the limits of integration
//double a = 0., b = 1.; // limits of integration
double abserr = 0.0; 
double relerr = 1e-7; // requested errors
double result; // the integral value
double error; // the error estimate
double n = 1;
size_t np = 1000; // work area size
gsl_integration_workspace *w = gsl_integration_workspace_alloc (np);
gsl_function F;
 //F.function = &ff;
 F.function = f1; // already a pointer
 F.params = &n;
gsl_integration_qag (&F, a, b, abserr, relerr, np, GSL_INTEG_GAUSS15, w, &result, &error);

 //printf ("result = % .18f\n", result);
 //printf ("estimated error = % .18f\n", error);
//printf ("intervals = %zu\n", w->size);

gsl_integration_workspace_free (w);
return result;
}

double fcos0 (double x, void *params)
 {
 double n = *(double *) params; 
 double f=cos( n*x );
 return f;
 }

double fcos1 (double x, void *params)
 {
  double a5 = ((double *)params)[0]; double a4 = ((double *)params)[1]; double a3 = ((double *)params)[2]; 
  double a2 = ((double *)params)[3]; double a1 = ((double *)params)[4]; double a0 = ((double *)params)[5];
 //double n = *(double *) params;
  double f=cos( a5*pow(x,5) + a4*pow(x,4) + a3*pow(x,3) + a2*pow(x,2) + a1*pow(x,1) + a0*pow(x,0) );
 return f;
 }

double fcos2 (double x, void *params) {
  std::vector <double> psi;
  for (int i=0;i<PSI_COEFF_LENGTH;i++){
  psi.push_back(((double *)params)[i]);
  }
   double f=cos (polyval(psi, x));
  return f;
 }

double fcos3 (double x, void *params) {
  std::vector <double> bks;
  std::vector <std::vector <double> > psi_cofs;
// first and second memory address refers to # of breaks and # of curver of the turn/ route respectively
  int num_breaks =(int)((double *)params)[0];
  int num_curves =(int)((double *)params)[1];

    //for (int i=0;i<NUM_BREAKS_FOR_TURN;i++)
    for (int i=0;i<num_breaks;i++)
    {
    bks.push_back(((double *)params)[i+2]);
    
  }
 
  for (int i=0;i<num_curves;i++)
  {
    std::vector <double> psi_cof;
    for(int j=0;j<PSI_COEFF_LENGTH;j++){
     psi_cof.push_back(((double *)params)[(num_breaks+2)+j+(i*PSI_COEFF_LENGTH)]);        
     }     
     psi_cofs.push_back( psi_cof);
  }
 
  double f=cos (ppval(psi_cofs, bks, x));
  return f;
 }

 double fsin1 (double x, void *params)
 {
  double a5 = ((double *)params)[0]; double a4 = ((double *)params)[1]; double a3 = ((double *)params)[2]; 
  double a2 = ((double *)params)[3]; double a1 = ((double *)params)[4]; double a0 = ((double *)params)[5];
 //double n = *(double *) params;
 double f=sin( a5*pow(x,5) + a4*pow(x,4) + a3*pow(x,3) + a2*pow(x,2) + a1*pow(x,1) + a0*pow(x,0) );
 return f;
 }

double fsin2 (double x, void *params) {
  std::vector <double> psi;
  for (int i=0;i<PSI_COEFF_LENGTH;i++){
  psi.push_back(((double *)params)[i]);
  }
  double f=sin (polyval(psi, x));
  return f;
 }
double fsin3 (double x, void *params) {
  std::vector <double> bks;
  std::vector <std::vector <double> > psi_cofs;
// first and second memory address refers to # of breaks and # of curver of the turn/ route respectively
  int num_breaks =(int)((double *)params)[0];
  int num_curves =(int)((double *)params)[1];

  for (int i=0;i<num_breaks;i++) {
    bks.push_back(((double *)params)[i+2]);
    }
  for (int i=0;i<num_curves;i++)
  {
    std::vector <double> psi_cof;
    for(int j=0;j<PSI_COEFF_LENGTH;j++){
     psi_cof.push_back(((double *)params)[(num_breaks+2)+j+(i*PSI_COEFF_LENGTH)]);         
     }
     psi_cofs.push_back( psi_cof);
  }
  
  double f=sin (ppval(psi_cofs, bks, x));
   return f;
 }

double time_fn (double x, void *params) {
  std::vector <double> vel_poly;
  for (int i=0;i<TURN_VEL_COEFF_LENGTH;i++){
  vel_poly.push_back(((double *)params)[i]);
  }
  double f= 1.0/(polyval(vel_poly, x));
  return f;
 }

double integrater2(double (*my_f) (double, void *params), std::vector <double> psi_cof, const double& a, const double& b)
{
// a and b are the limits of integration
//double a = 0., b = 1.; // limits of integration
//double abserr = 0.0; 
double abserr = 1e-7;
double relerr = 1e-7; // requested errors
double result; // the integral value
double error; // the error estimate

size_t np = 10000; // work area size
gsl_integration_workspace *w = gsl_integration_workspace_alloc (np);
gsl_function F;

F.function = my_f; // already a pointer
double* prm=psi_cof.data();
F.params=prm;
 
gsl_integration_qag (&F, a, b, abserr, relerr, np, GSL_INTEG_GAUSS15, w, &result, &error);
//gsl_integration_qags (&F, a, b, abserr, relerr, np, w, &result, &error);
gsl_integration_workspace_free (w);
return result;
}

point_xyz get_turn_end_point(const turn_struct& turn, const double& wind){
  //We consider wind as 0.0 for now...

  point_xyz x_end; double x_integral=0.0; double y_integral=0.0;
  std::vector <double> psi1; std::vector <double> psi2;  std::vector <double> psi3;

  psi1= polyint(turn.Cs.at(0), turn.curv_heading1);   
  x_integral = integrater2(fcos2, psi1, 0.0, turn.Sfs.at(0)); // Sf1
  y_integral = integrater2(fsin2, psi1, 0.0, turn.Sfs.at(0));

  psi2= polyint(turn.Cs.at(1), polyval(psi1,turn.Sfs.at(0)));//Sf1
  x_integral += integrater2(fcos2, psi2,0.0, turn.Sfs.at(1)); //Sf2
  y_integral += integrater2(fsin2, psi2, 0.0, turn.Sfs.at(1));

  psi3= polyint(turn.Cs.at(2), polyval(psi2,turn.Sfs.at(1))); //Sf2
  x_integral += integrater2(fcos2, psi3, 0.0, turn.Sfs.at(2));//Sf3
  y_integral += integrater2(fsin2, psi3, 0.0, turn.Sfs.at(2));

  x_end.x=x_integral; x_end.y=y_integral; x_end.z=0.0;

return x_end;
  }
/////////////////////////////////////////////////////////////////////

double dis_pt_to_pt_2D(const point_xyz& pt1, const point_xyz& pt2){ 
 return sqrt(pow((pt2.x-pt1.x),2.0) + pow((pt2.y-pt1.y),2.0))  ;
}

double dis_pt_to_pt_3D(const point_xyz& pt1, const point_xyz& pt2){ 
 return sqrt(pow((pt2.x-pt1.x),2.0) + pow((pt2.y-pt1.y),2.0) + pow((pt2.z-pt1.z),2.0)) ;
}

double dis_frm_line_to_pt(const point_xyz& pt, const point_xyz& pt1, const point_xyz& pt2){
  double a=(pt2.y-pt1.y); double b=-(pt2.x-pt1.x); double c= (pt1.y*(pt2.x-pt1.x) - pt1.x*(pt2.y-pt1.y));   
  return abs(a*pt.x + b*pt.y +c)/sqrt(pow(a,2.0)+pow(b,2.0));
}

bool check_pt_bw_segments(const point_xyz& p, const point_xyz& pt1, const point_xyz& pt2, const point_xyz& pt3){
// Check wehther point p is on the same side of the line segments defiend by pt1, pt2, and pt3
  int sgn1=0; int sgn2=0; bool flag=false;

  sgn1 = (((pt2.x-pt1.x)*(pt3.y-pt2.y) -(pt2.y-pt1.y)*(pt3.x-pt2.x))<0.0) ? -1 : 1;
  sgn2= (((pt3.x-pt2.x)*(p.y-pt2.y) -(pt3.y-pt2.y)*(p.x-pt2.x))<0.0)? -1: 1;

  flag = (sgn1==sgn2) ? true : false ;
  return flag; 
}

find_curv_poly_fit_out_struct find_curv_poly_fit(const point_xyz& curv_end_pt, const point_xyz& x0,const point_xyz& curve_wpt2,const point_xyz& xlim, const double& dist_tol)
{ // For using only fly by way points
const int max_iter=5000;  float low =0.0; float high=1.0; float t=0.5; int num_iter=0;
find_curv_poly_fit_out_struct st1;
point_xyz x_start; point_xyz x_end; point_xyz wpt;
wpt.x=curve_wpt2.x; wpt.y=curve_wpt2.y; wpt.z=curve_wpt2.z;

while ((t < 1.0) && (num_iter<max_iter) ) {
  t=(low+high)/2;
  x_start.x= x0.x + t*(curve_wpt2.x-x0.x);
  x_start.y= x0.y + t*(curve_wpt2.y-x0.y);
  x_end.x= x_start.x + curv_end_pt.x;
  x_end.y= x_start.y + curv_end_pt.y;

// check if we are converged
if (dis_frm_line_to_pt(x_end,wpt,xlim)<dist_tol) {
//  std::cout<<"dist: "<< dis_frm_line_to_pt(x_end,wpt,xlim)<< std::endl;
  if (dis_pt_to_pt_2D(x_end,wpt)> dis_pt_to_pt_2D(xlim,wpt) ){
    std::cout<<"Curve poly NOT matched** "<< std::endl;
    st1.matched=false;
    } 
  else
    {
      std::cout<<"Curve poly matched** "<< std::endl;
      st1.x_init=x_start; st1.x_end=x_end; st1.matched=true;
    }  
  break; 
}
else {  // Shift start point otherwise
  if( check_pt_bw_segments(x_end,x0,wpt,xlim)){
    low = t;
  } else {
    high = t;
  }
num_iter+=1;
} 
}
return st1;
//std::cin.get();
}

void setup_piecewice_poly_in_turn(turn_struct& turn, const point_xyz& x_init, const point_xyz& x_end){
turn.start= x_init;
turn.end= x_end;
turn.s_breaks.clear(); turn.s_breaks.push_back(0);
turn.t_breaks.clear();  turn.t_breaks.push_back(0);
turn.curv_coeffs.clear(); turn.psi_coeffs.clear(); turn.t_coeffs.clear();
double prev_psi=turn.curv_heading1; 
//double prev_S=0.0;

if (turn.Sfs.at(0)>EPS){    //Sf1
  turn.s_breaks.push_back((turn.s_breaks.back()+turn.Sfs.at(0)));
  turn.curv_coeffs.push_back(turn.Cs.at(0));
  turn.psi_coeffs.push_back(polyint(turn.Cs.at(0),prev_psi)); // polynomial integration
  prev_psi=polyval(turn.psi_coeffs.back(),turn.Sfs.at(0)); // polynomial evaluation
 // prev_S +=turn.Sfs.at(0);
  
  turn.t_breaks.push_back( (turn.Sfs.at(0))/turn.vel_poly.back() ); // turn.vel_poly.back() has the const. airspeed
  std::vector <double> t_cof; t_cof.push_back(0); t_cof.push_back(turn.vel_poly.back()); t_cof.push_back(0); 
  turn.t_coeffs.push_back(t_cof);
}

if (turn.Sfs.at(1)>EPS){ //Sf2
  turn.s_breaks.push_back(turn.Sfs.at(1)+turn.s_breaks.back());
  turn.curv_coeffs.push_back(turn.Cs.at(1));
  turn.psi_coeffs.push_back(polyint(turn.Cs.at(1),prev_psi)); // polynomial integration
  prev_psi=polyval(turn.psi_coeffs.back(),turn.Sfs.at(1)); // polynomial evaluation
  //prev_S +=turn.Sfs.at(1);  

  turn.t_breaks.push_back( ((turn.Sfs.at(0))/turn.vel_poly.back()) + ((turn.Sfs.at(1))/turn.vel_poly.back()) ); // turn.vel_poly.back() has the const. airspeed
  std::vector <double> t_cof; t_cof.push_back(0); t_cof.push_back(turn.vel_poly.back()); t_cof.push_back(turn.Sfs.at(0)); 
  turn.t_coeffs.push_back(t_cof);  
}

if (turn.Sfs.at(2)>EPS){ //Sf3
  turn.s_breaks.push_back(turn.Sfs.at(2)+turn.s_breaks.back());
  turn.curv_coeffs.push_back(turn.Cs.at(2));
  turn.psi_coeffs.push_back(polyint(turn.Cs.at(2),prev_psi)); // polynomial integration

  turn.t_breaks.push_back( ((turn.Sfs.at(0))/turn.vel_poly.back()) + ((turn.Sfs.at(1))/turn.vel_poly.back()) + ((turn.Sfs.at(2))/turn.vel_poly.back()) ); // turn.vel_poly.back() has the const. airspeed
  std::vector <double> t_cof; t_cof.push_back(0); t_cof.push_back(turn.vel_poly.back()); t_cof.push_back(turn.Sfs.at(0)+ turn.Sfs.at(1)); 
  turn.t_coeffs.push_back(t_cof);  
}
}

std::vector<double> vector_divide( const std::vector<double>& a, const std::vector<double>& b )
{
    std::vector<double> result ;
    const std::size_t n = std::min( a.size(), b.size() ) ;
    std::transform( std::begin(a), std::begin(a)+n, std::begin(b),std::back_inserter(result), std::divides<>{} ) ;
    return result ;
}

void get_path_from_turn(turn_struct& turn, const double& wind, const int& num_samples){

// get heading function  
 std::vector <double> breaks_and_coffs; // making one long vector of brks and psi cofs to pass into integrator 3

 breaks_and_coffs.push_back(NUM_BREAKS_FOR_TURN); breaks_and_coffs.push_back(NUM_CURVS_FOR_TURN);

 breaks_and_coffs.insert(breaks_and_coffs.end(), turn.s_breaks.begin(), turn.s_breaks.end());
 for(int i=0;i<NUM_CURVS_FOR_TURN;i++){ // NUM_CURVS_FOR_TURN is same as turn.psi_coeffs.size() here
    breaks_and_coffs.insert(breaks_and_coffs.end(), turn.psi_coeffs[i].begin(), turn.psi_coeffs[i].end());
   }
   
 double sample_step= turn.Sf_total/(num_samples-1);

 turn.path={}; //reset before filing again
for (int i=0; i<num_samples; i++){
  double ith_sample_step=sample_step*i; //double x_integral=0.0; double y_integral=0.0;
  point_xyz pt_p;
 
  turn.path.path_velocity.push_back(polyval(turn.vel_poly, ith_sample_step));
  turn.path.path_heading.push_back( ppval(turn.psi_coeffs,turn.s_breaks,ith_sample_step));
  turn.path.path_curv.push_back( ppval(turn.curv_coeffs,turn.s_breaks,ith_sample_step));
  turn.path.path_bank.push_back( atan2(((pow(turn.path.path_velocity.back(),2.0))*turn.path.path_curv.back()),GRAVITY));
  turn.path.path_times.push_back(integrater2(time_fn, turn.vel_poly, 0.0, ith_sample_step));
  pt_p.x = integrater2(fcos3, breaks_and_coffs, 0.0, ith_sample_step) + turn.start.x; //+ x_direction of wind_speed* time(s);
  pt_p.y = integrater2(fsin3, breaks_and_coffs, 0.0, ith_sample_step) + turn.start.y; //+ y_direction of wind_speed* time(s);
  turn.path.path_xyz.push_back(pt_p);
}
std::vector <double> path_bank_diff=turn.path.path_bank; std::vector <double> path_times_diff=turn.path.path_times; 
std::adjacent_difference(path_bank_diff.begin(), path_bank_diff.end(), path_bank_diff.begin());
std::adjacent_difference( path_times_diff.begin(),  path_times_diff.end(), path_times_diff.begin());
path_bank_diff.erase(path_bank_diff.begin()); path_times_diff.erase(path_times_diff.begin());
turn.path.path_roll_rate=vector_divide(path_bank_diff,path_times_diff);
}

bool tunnel_check(const std::vector <point_xyz>& turn_xyz, const std::vector <point_xyz>& tnl1, const std::vector <point_xyz>& tnl2){
  point_xyz A1 = tnl1.at(0); point_xyz B1 = tnl1.at(1); point_xyz C1 = tnl1.at(2);
  point_xyz A2 = tnl2.at(0); point_xyz B2 = tnl2.at(1); point_xyz C2 = tnl2.at(2);

  double dot_A1_B1_A1_B1=  (B1.x-A1.x)*(B1.x-A1.x) + (B1.y-A1.y)*(B1.y-A1.y);
  double dot_B1_C1_B1_C1=  (C1.x-B1.x)*(C1.x-B1.x) + (C1.y-B1.y)*(C1.y-B1.y);

  double dot_A2_B2_A2_B2=  (B2.x-A2.x)*(B2.x-A2.x) + (B2.y-A2.y)*(B2.y-A2.y);
  double dot_B2_C2_B2_C2=  (C2.x-B2.x)*(C2.x-B2.x) + (C2.y-B2.y)*(C2.y-B2.y);
  
  std::vector <int> tnl_chk1;
  for (int i=0;i<turn_xyz.size();i++){
    point_xyz p=turn_xyz.at(i);
    double dot_B1A1_ptA1= (p.x-A1.x)*(B1.x-A1.x) + (p.y-A1.y)*(B1.y-A1.y);  
    double dot_C1B1_ptB1= (p.x-B1.x)*(C1.x-B1.x) + (p.y-B1.y)*(C1.y-B1.y);

    double dot_B2A2_ptA2= (p.x-A2.x)*(B2.x-A2.x) + (p.y-A2.y)*(B2.y-A2.y);  
    double dot_C2B2_ptB2= (p.x-B2.x)*(C2.x-B2.x) + (p.y-B2.y)*(C2.y-B2.y);

  if( (dot_B1A1_ptA1<0) || (dot_C1B1_ptB1<0) || (dot_B1A1_ptA1>dot_A1_B1_A1_B1) || (dot_C1B1_ptB1>dot_B1_C1_B1_C1) ){
    if((dot_B2A2_ptA2<0) || (dot_C2B2_ptB2<0) || (dot_B2A2_ptA2>dot_A2_B2_A2_B2) || (dot_C2B2_ptB2>dot_B2_C2_B2_C2) ) {
          tnl_chk1.push_back(0) ; 
       }else {
        tnl_chk1.push_back(1) ;  
       } 
  }else{
    tnl_chk1.push_back(1) ;
  }
  }
return (std::all_of(tnl_chk1.begin(), tnl_chk1.end(), [](int i) { return i==1; })) ;
}

double set_max_airspeed(const double& mas, const double& step)    {
    double max_airspeed= floor(mas);
    int mod_result=(static_cast<uint64_t> (max_airspeed)) % int(step);
    for (int i=0; i<mod_result; i++){
      max_airspeed-=1.0;
      }
    return max_airspeed;
  }

struct turn_struct generate_fly_by_turn(const route_struct& route, const int& turn_i, const curv_struct& curvs){

  bool VALID_TURN = false; // To test the turn is valid and move  forward..
  point_xyzvdldrfo curve_wpt1; point_xyzvdldrfo curve_wpt2; point_xyzvdldrfo curve_wpt3;
  // PUSHING X, Y, Z, AND VELOCITY FOR EACH CURVE WAYPOINT
  curve_wpt1=route.wpt_list[turn_i]; curve_wpt2=route.wpt_list[turn_i+1]; curve_wpt3=route.wpt_list[turn_i+2];

  double heading1=atan2((curve_wpt2.y-curve_wpt1.y),(curve_wpt2.x-curve_wpt1.x));
  double heading2=atan2((curve_wpt3.y-curve_wpt2.y),(curve_wpt3.x-curve_wpt2.x));
  heading2=(heading2 >= 0.0 ? heading2 : (2*M_PI + heading2)); // radian value, wrapTo2Pi
  heading1=(heading1 >= 0.0 ? heading1 : (2*M_PI + heading1));

  double head_change=abs(heading2-heading1); 
  head_change=(head_change > M_PI ? (2*M_PI-head_change) : head_change);
  double temp=(cos(heading1)*sin(heading2)-sin(heading1)*cos(heading2));
  double sign=(temp > 0) ? 1 : ((temp < 0) ? -1 : 0);
  head_change=sign*head_change;
  turn_struct turn;

  //if (head_change !=0.0)
  if (abs(head_change)>HEAD_CHANGE_THRESH) // else consider head change=0
{
    // GENERATE CORRIDORS.. GET TUNNELS....
  std::vector <point_xyz> tunnel1; // holds 4 points for tunnel 
  std::vector <point_xyz> tunnel2;
  tunnel1=generate_corridor(curve_wpt1,curve_wpt2);
  tunnel2=generate_corridor(curve_wpt2,curve_wpt3);

  // Iterate over speeds until we find a feasible turn...
  double max_airspeed= set_max_airspeed(curve_wpt2.v, curvs.vel_step) ; 
  point_xyz x0; point_xyz xlim;
  x0.x=curve_wpt1.x +0.5*(curve_wpt2.x-curve_wpt1.x); x0.y=curve_wpt1.y +0.5*(curve_wpt2.y-curve_wpt1.y); 
  xlim.x=curve_wpt2.x +0.75*(curve_wpt3.x-curve_wpt2.x); xlim.y=curve_wpt2.y +0.75*(curve_wpt3.y-curve_wpt2.y); 
    
  while (max_airspeed>=curvs.min_vel)
  {
        
  // model struct is already updated with the curve speed...
    int count1=0; // As the  poly_struct poly strarts from 0 instead of 1 as in the curve_poly.txt
    bool flag=false;
    while (count1<static_cast<int>(curvs.num_curv/2)) // 5 of them are - curvatures and 5 of them are + curvatures
    {
       std::cout<<" Curve count number "<<count1<< std::endl;
       std::cout<<" head_change "<<head_change<< std::endl;
      flag=false;  
       //std::cin.get(); 
    poly_struct C1C3Sf1Sf3k_trans; // struct for holding retun values from following lookup function
    C1C3Sf1Sf3k_trans=lookup_curve_splines(((head_change>0)? 1:0),count1, max_airspeed,curvs);

    // FIRST WE CONTINUE WITHOUT FLY OVER WAY POINTS 11/01/2021
       initial_Sf2_out_struct Sf2_rem_change= get_initial_sf2(C1C3Sf1Sf3k_trans,head_change);

  if(Sf2_rem_change.rem_change>=0.0){
    turn.Sfs.clear();   turn.Cs.clear();    // Clear the vectors before refilling    
    turn.Sfs.push_back(C1C3Sf1Sf3k_trans.Sf1); turn.Sfs.push_back(Sf2_rem_change.Sf2); turn.Sfs.push_back(C1C3Sf1Sf3k_trans.Sf2);
    turn.Sf_total=C1C3Sf1Sf3k_trans.Sf1+Sf2_rem_change.Sf2+C1C3Sf1Sf3k_trans.Sf2;
    turn.Cs.push_back(C1C3Sf1Sf3k_trans.C1);
    std::vector <double> C2(CURV_POLY_DEG,0.0); 
    C2.push_back(C1C3Sf1Sf3k_trans.ktrans);  turn.Cs.push_back(C2);
    turn.Cs.push_back(C1C3Sf1Sf3k_trans.C2);
    turn.ktrans=C1C3Sf1Sf3k_trans.ktrans;
    turn.curv_heading1=heading1;
    turn.curv_heading2=heading2;
    turn.empty=false;
    std::vector <double> vp; vp.push_back(0); vp.push_back(max_airspeed);
    turn.vel_poly=vp; // Assuming constant speed from start to end of the turn. It can be changed by adding a polynomial (now accelraiton =0)
    flag=true;

    
  }
  
    //if (flag == false) continue;
    if (flag == false) {count1++; continue;} //MODIFIED ON 29/04.2021

     
    point_xyz x_end=get_turn_end_point(turn, WIND_SPEED);  // correct
    point_xyz curv_wpt2_xyz; curv_wpt2_xyz.x=curve_wpt2.x; curv_wpt2_xyz.y=curve_wpt2.y; curv_wpt2_xyz.z=curve_wpt2.z;
  // Following is assuming no fly-over way point, if it exists use find_curv_poly_fit2 to release the corridor check 
  find_curv_poly_fit_out_struct poly_fit_out =find_curv_poly_fit(x_end, x0,curv_wpt2_xyz,xlim,DIST_TOL);
   if (poly_fit_out.matched==true){ // Continue as fly-by waypoint only first

    setup_piecewice_poly_in_turn(turn, poly_fit_out.x_init, poly_fit_out.x_end);     
    get_path_from_turn(turn, 0.0, 100); // inputs: turn struct, wind, number of samples

    if(tunnel_check(turn.path.path_xyz, tunnel1, tunnel2)){
      turn.vel= max_airspeed; // Assuming a constant velocity throughout the turn.(= current max_airspeed) 
      turn.curve_max=GRAVITY*tan(max_roll)/(max_airspeed*max_airspeed);
      turn.curve_rate_max= GRAVITY*(max_roll_rate)/(max_airspeed*max_airspeed*max_airspeed) - 2*turn.curve_max*max_accel/(max_airspeed*max_airspeed);
      turn.curve_rate_rate_max= GRAVITY*(max_roll_rate_rate)/(max_airspeed*max_airspeed*max_airspeed*max_airspeed);
      turn.empty =false;
      VALID_TURN = true;
      std::cout<<" A valid turn is found for fly-by curve: "<<turn_i<< std::endl;
      break;}
    else {
     VALID_TURN=false;
     std::cout<<" A valid turn is NOT found for fly-by curve: "<<turn_i<< std::endl;
     }  
  }
  count1++;
      }
  if (VALID_TURN == true) break;    
  max_airspeed -= curvs.vel_step;
  if (max_airspeed<curvs.min_vel) turn.empty=true;
  }
}
else
{ // If head change =0, a turn with no head change
 turn.vel=curve_wpt2.v;
 turn.Sf_total=0.0;
 turn.empty =false;
 turn.start.x=curve_wpt2.x; turn.start.y=curve_wpt2.y; turn.start.z=curve_wpt2.z;
 turn.end=turn.start;
}
//turn.empty = (VALID_TURN == true) ? (false) : (true); //already assigned above
return turn;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct turn_struct generate_fly_over_turn(const route_struct& route, const int& turn_i, const curv_struct& curvs){
  double fowp_ang= M_PI/4;//M_PI/6;//M_PI/4; 
  double fowp_theta= 0.0;
  turn_struct turn1;  turn_struct turn2;  turn_struct fowp_turn; // turn to represent both turn1 and turn2 and the straight line between
  fowp_turn.isfowp=true;

  bool VALID_TURN1 = false; bool VALID_TURN2 = false; // To test the turn is valid and move  forward..
  point_xyzvdldrfo curve_wpt1; point_xyzvdldrfo curve_wpt2; point_xyzvdldrfo curve_wpt3;
  // PUSHING X, Y, Z, AND VELOCITY FOR EACH CURVE WAYPOINT
  curve_wpt1=route.wpt_list[turn_i]; curve_wpt2=route.wpt_list[turn_i+1]; curve_wpt3=route.wpt_list[turn_i+2];

  double heading1=atan2((curve_wpt2.y-curve_wpt1.y),(curve_wpt2.x-curve_wpt1.x));
  double heading2=atan2((curve_wpt3.y-curve_wpt2.y),(curve_wpt3.x-curve_wpt2.x));
  heading2=(heading2 >= 0.0 ? heading2 : (2*M_PI + heading2)); // radian value, wrapTo2Pi
  heading1=(heading1 >= 0.0 ? heading1 : (2*M_PI + heading1));

  double head_change=abs(heading2-heading1); 
  head_change=(head_change > M_PI ? (2*M_PI-head_change) : head_change);
  double temp=(cos(heading1)*sin(heading2)-sin(heading1)*cos(heading2));
  double sign=(temp > 0) ? 1 : ((temp < 0) ? -1 : 0);
  head_change=sign*head_change;

  if(head_change>0){
    head_change += fowp_ang;
    fowp_theta =heading2+fowp_ang;
  }
  else if(head_change<0){
    head_change -= fowp_ang;
    fowp_theta =heading2-fowp_ang;
  }
  else{  
    // DO NOTHING//
  }  
//////////////////////////Check corridor for for turn1 only!.. 
  std::vector <point_xyz> tunnel1=generate_corridor(curve_wpt1,curve_wpt2);  // holds 4 points for tunnel
  std::vector <point_xyz> tunnel2=generate_corridor(curve_wpt2,curve_wpt3);
  // Iterate over speeds until we find a feasible turn

  double max_airspeed= set_max_airspeed(curve_wpt2.v, curvs.vel_step) ; 
  point_xyz x0; point_xyz xlim;
  x0.x=curve_wpt1.x +0.5*(curve_wpt2.x-curve_wpt1.x); x0.y=curve_wpt1.y +0.5*(curve_wpt2.y-curve_wpt1.y); 
  xlim.x=curve_wpt2.x +0.5*(curve_wpt3.x-curve_wpt2.x); xlim.y=curve_wpt2.y +0.5*(curve_wpt3.y-curve_wpt2.y); 
    
  while (max_airspeed>=curvs.min_vel)
  {
    std::cout<<" Max_airspeed: "<<max_airspeed<< std::endl;    
  // model struct is already updated with the curve speed...
    int count1=0; // As the  poly_struct poly strarts from 0 instead of 1 as in the curve_poly.txt
    bool flag1=false;
    while (count1<static_cast<int>(curvs.num_curv/2)) // 5 of them are - curvatures and 5 of them are + curvatures
    {
      std::cout<<" Curve count number for count 1: "<<count1<< std::endl;
      flag1=false;           
   // struct for holding retun values from following lookup function
    poly_struct C1C3Sf1Sf3k_trans1=lookup_curve_splines(((head_change>0)? 1:0),count1, max_airspeed,curvs);
    initial_Sf2_out_struct Sf2_rem_change1= get_initial_sf2(C1C3Sf1Sf3k_trans1,head_change);

  if(Sf2_rem_change1.rem_change>=0.0){
    turn1.Sfs.clear(); turn1.Cs.clear(); // Clear the vector before refilling.
    turn1.Sfs.push_back(C1C3Sf1Sf3k_trans1.Sf1); turn1.Sfs.push_back(Sf2_rem_change1.Sf2); turn1.Sfs.push_back(C1C3Sf1Sf3k_trans1.Sf2);
    turn1.Sf_total=C1C3Sf1Sf3k_trans1.Sf1+Sf2_rem_change1.Sf2+C1C3Sf1Sf3k_trans1.Sf2;
    turn1.Cs.push_back(C1C3Sf1Sf3k_trans1.C1);
    std::vector <double> C2(CURV_POLY_DEG,0.0); 
    C2.push_back(C1C3Sf1Sf3k_trans1.ktrans);  turn1.Cs.push_back(C2);
    turn1.Cs.push_back(C1C3Sf1Sf3k_trans1.C2);
    turn1.ktrans=C1C3Sf1Sf3k_trans1.ktrans;
    turn1.curv_heading1=heading1;
    turn1.curv_heading2=heading2;
    turn1.empty=false;
    std::vector <double> vp; vp.push_back(0); vp.push_back(max_airspeed);
    turn1.vel_poly=vp; // Assuming constant speed from start to end of the turn. It can be changed by adding a polynomial (now accelraiton =0)
    turn1.vel=max_airspeed;
    flag1=true;
  }
    //if (flag1 == false) continue;
    if (flag1 == false) {count1++; continue;} //MODIFIED ON 29/04.2021
      
      point_xyz x_end1=get_turn_end_point(turn1, WIND_SPEED);  // correct
      point_xyz x_start1; x_start1.x=curve_wpt2.x; x_start1.y=curve_wpt2.y; x_start1.z=curve_wpt2.z; // Consider start point as the fowp
      x_end1.x+= x_start1.x; x_end1.y+= x_start1.y; x_end1.z+= x_start1.z;
    // Add turn1 to breaks..
      setup_piecewice_poly_in_turn(turn1, x_start1, x_end1);

    //////////////////////////////////////////////////////////////// TUNNEL CHECK FOR FOWP (FOR ONLY TURN 1 SECTION)
    get_path_from_turn(turn1, 0.0, 100); // inputs: turn struct, wind, number of samples
    if(tunnel_check(turn1.path.path_xyz, tunnel1, tunnel2)){
      turn1.vel= max_airspeed; // Assuming a constant velocity throughout the turn.(= current max_airspeed) 
      turn1.curve_max=GRAVITY*tan(max_roll)/(max_airspeed*max_airspeed);
      turn1.curve_rate_max= GRAVITY*(max_roll_rate)/(max_airspeed*max_airspeed*max_airspeed) - 2*turn1.curve_max*max_accel/(max_airspeed*max_airspeed);
      turn1.curve_rate_rate_max= GRAVITY*(max_roll_rate_rate)/(max_airspeed*max_airspeed*max_airspeed*max_airspeed);
      turn1.empty =false;
      VALID_TURN1 = true;
      std::cout<<" A valid turn 1 is found for fly-over curve: "<<turn_i<< std::endl;
      break;}
    else {
     VALID_TURN1=false;
     std::cout<<" A valid turn 1 is NOT found for fly-over curve: "<<turn_i<< std::endl;
     }  
  
  count1++;
      }
  if (VALID_TURN1 == true) break;    
  max_airspeed -= curvs.vel_step;
  }
  if (VALID_TURN1 == false) return fowp_turn;

  // First add turn 1 coefficients to fowp turn, still need to add straight line and turn 2 parameters...
fowp_turn.t_breaks=turn1.t_breaks; fowp_turn.s_breaks=turn1.s_breaks; fowp_turn.t_coeffs=turn1.t_coeffs;  
fowp_turn.curv_coeffs=turn1.curv_coeffs; fowp_turn.psi_coeffs=turn1.psi_coeffs; fowp_turn.Sfs=turn1.Sfs;

// NOW PERFORM THE SECOND TURN
  turn2.vel_poly= turn1.vel_poly; turn2.vel= turn1.vel;
  // Get the eqn. of the straight line  fowpxy=turn1.end
  double c1=turn1.end.y- tan(fowp_theta)*turn1.end.x;
  point_xyz xy_intersect;
  double c2= curve_wpt2.y-tan(heading2)*curve_wpt2.x;
  xy_intersect.x=(c2-c1)/(tan(fowp_theta)-tan(heading2)); // (c2-c1)/(m1-m2)
  xy_intersect.y= tan(fowp_theta)*xy_intersect.x +c1;

  // Calculate new head change due to the intersection point:
  double nheading2=atan2((curve_wpt3.y-xy_intersect.y), (curve_wpt3.x-xy_intersect.x)); 
  double nheading1=atan2((xy_intersect.y-turn1.end.y), (xy_intersect.x-turn1.end.x));
  double nhead_change=abs(nheading2-nheading1);
  nhead_change=(nhead_change > M_PI ? (2*M_PI-nhead_change) : nhead_change);
  
   int count2=0; // As the  poly_struct poly strarts from 0 instead of 1 as in the curve_poly.txt
    bool flag2=false;
    while (count2<static_cast<int>(curvs.num_curv/2)) // 5 of them are - curvatures and 5 of them are + curvatures
    {
       std::cout<<" Curve count number for count 2: "<<count2<< std::endl;
      flag2=false;  

      poly_struct C1C3Sf1Sf3k_trans2=lookup_curve_splines(((head_change>0)? 0:1),count2, turn2.vel,curvs); // change the direction of second curve
      initial_Sf2_out_struct Sf2_rem_change2= get_initial_sf2(C1C3Sf1Sf3k_trans2,nhead_change);

      if(Sf2_rem_change2.rem_change>=0.0){         
        turn2.Sfs.clear(); turn2.Cs.clear(); // Clear the vector before refilling.
        turn2.Sfs.push_back(C1C3Sf1Sf3k_trans2.Sf1); turn2.Sfs.push_back(Sf2_rem_change2.Sf2); turn2.Sfs.push_back(C1C3Sf1Sf3k_trans2.Sf2);
        turn2.Sf_total=C1C3Sf1Sf3k_trans2.Sf1+Sf2_rem_change2.Sf2+C1C3Sf1Sf3k_trans2.Sf2;
        turn2.Cs.push_back(C1C3Sf1Sf3k_trans2.C1);        
        std::vector <double> C2(CURV_POLY_DEG,0.0); 
        C2.push_back(C1C3Sf1Sf3k_trans2.ktrans); turn2.Cs.push_back(C2);
        turn2.Cs.push_back(C1C3Sf1Sf3k_trans2.C2);
        turn2.ktrans=C1C3Sf1Sf3k_trans2.ktrans;
        turn2.curv_heading1=nheading1;
        turn2.curv_heading2=nheading2;
        turn2.empty=false;
        flag2=true;
       }
      if (flag2 == false) continue;  

      point_xyz x_end_fowp=get_turn_end_point(turn2, WIND_SPEED);  // Get the endpoint of the second turn of fowp
      point_xyz new_x0; new_x0.x=  (turn1.end.x+xy_intersect.x)/2; new_x0.y=  (turn1.end.y+xy_intersect.y)/2;
      point_xyz new_xlim; new_xlim.x=(curve_wpt3.x+xy_intersect.x)/2; new_xlim.y=(curve_wpt3.y+xy_intersect.y)/2;       
      find_curv_poly_fit_out_struct poly_fit_out2 =find_curv_poly_fit(x_end_fowp, new_x0,xy_intersect,new_xlim,DIST_TOL);

      if (poly_fit_out2.matched==true){ //     
        setup_piecewice_poly_in_turn(turn2, poly_fit_out2.x_init,poly_fit_out2.x_end);
        double distance_traveled=turn1.Sf_total;
        double dist_t2_end_t1_start= dis_pt_to_pt_2D(turn2.start, turn1.end); //(turn2.start-turn1.end) distance  has to follow as a straight line
        double prev_psi= polyval(turn1.psi_coeffs.back(),turn1.Sfs.back());
   
          if(dist_t2_end_t1_start>DIST_TOL){ // Add straight line coeffs to fowp turn
            fowp_turn.Sfs.push_back(dist_t2_end_t1_start);    //turn2.vel is the velocity of the total fowp turn including 2 turns and straight line between
            fowp_turn.s_breaks.push_back(fowp_turn.s_breaks.back()+dist_t2_end_t1_start);
            fowp_turn.t_breaks.push_back(fowp_turn.t_breaks.back()+ (dist_t2_end_t1_start/turn2.vel));
            std::vector <double> v1; v1.push_back(0); v1.push_back(turn2.vel); v1.push_back(distance_traveled); fowp_turn.t_coeffs.push_back(v1); // the total dist. traveled before the striahgt line
            distance_traveled+=dist_t2_end_t1_start;
            std::vector <double> c_c;for (int i = 0; i < CURV_COEFF_LENGTH; i++) c_c.push_back(0);fowp_turn.curv_coeffs.push_back(c_c);
            std::vector <double> p_c; for (int j = 0; j < PSI_COEFF_LENGTH-1; j++) p_c.push_back(0); p_c.push_back(prev_psi); fowp_turn.psi_coeffs.push_back(p_c);
        }
      // Now add turn 2 coeffients to fowp turn...
            for (int k = 1; k < turn2.s_breaks.size(); k++) fowp_turn.s_breaks.push_back(fowp_turn.s_breaks.back()+turn2.s_breaks.at(k)-turn2.s_breaks.at(k-1));
            for (int k = 1; k < turn2.t_breaks.size(); k++) fowp_turn.t_breaks.push_back(fowp_turn.t_breaks.back()+turn2.t_breaks.at(k)-turn2.t_breaks.at(k-1));      
            for (int k = 0; k < turn2.curv_coeffs.size(); k++) fowp_turn.curv_coeffs.push_back(turn2.curv_coeffs.at(k));        
            for (int k=0; k< turn2.t_coeffs.size();k++){
              turn2.t_coeffs.at(k).back()+=distance_traveled; 
              fowp_turn.t_coeffs.push_back(turn2.t_coeffs.at(k));
              }
            for (int k=0; k< turn2.Sfs.size();k++){
              fowp_turn.psi_coeffs.push_back(polyint(turn2.Cs.at(k),prev_psi));
              prev_psi= polyval(fowp_turn.psi_coeffs.back(), turn2.Sfs.at(k));
              fowp_turn.Sfs.push_back(turn2.Sfs.at(k));
            }
            //Add other parameters to fowp_turn...
            fowp_turn.start=turn1.start; fowp_turn.end=turn2.end; fowp_turn.vel=turn2.vel;
            fowp_turn.Sf_total= turn1.Sf_total+dist_t2_end_t1_start+ turn2.Sf_total;  fowp_turn.vel_poly=turn2.vel_poly;  
            flag2 = true; fowp_turn.empty=false; 
            std::cout<<" A valid turn 2 is found for the fly-over curve.. "<< std::endl;
            break;
   }   
    count2++;
    }
    if (flag2 == false) {
      std::cout<<"A valid turn 2 is NOT found for the fly-over curve.. "<< std::endl;
      fowp_turn.empty=true;
    } 
return fowp_turn;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct straight_struct generate_straight_line(const route_struct& route, const int& straight_i, const bool& use_max_accel)
{
  straight_struct straight;  
 if (straight_i==0) {  
  point_xyz pt; 
  pt.x=route.wpt_list.front().x; pt.y=route.wpt_list.front().y; pt.z=route.wpt_list.front().z;
  straight.start=pt;
  straight.end=route.turns.front().start;
  straight.vel_start=route.wpt_list.front().v;
  straight.vel_end=route.turns.front().vel;
  
 }
 else if (straight_i==route.num_straight_lines-1) {
  point_xyz pt;  
  pt.x=route.wpt_list.back().x; pt.y=route.wpt_list.back().y; pt.z=route.wpt_list.back().z;
  straight.end=pt;
  straight.start=route.turns.back().end;
  straight.vel_start=route.turns.back().vel;
  straight.vel_end=route.wpt_list.back().v; 
}
else {
  straight.start=route.turns.at(straight_i-1).end;
  straight.end=route.turns.at(straight_i).start;  
  straight.vel_start=route.turns.at(straight_i-1).vel;
  straight.vel_end=route.turns.at(straight_i).vel;  
}
straight.heading=atan2((straight.end.y-straight.start.y), (straight.end.x-straight.start.x));
straight.heading=(straight.heading <0.0 ? (2*M_PI+straight.heading) : straight.heading);
straight.distance= dis_pt_to_pt_2D(straight.start, straight.end);
  double new_accel=(pow(straight.vel_end,2.0)- pow(straight.vel_start,2.0))/(2*straight.distance);
  
  if(new_accel <=  max_accel){
    straight.empty=false;    
        
    if(straight.vel_end>straight.vel_start) // Go on acceleraton t1 tmes and go on v2 velocity t2 times
      {
      straight.accel=(use_max_accel==true)? max_accel : new_accel;
      double t1 =   (straight.vel_end-straight.vel_start)/straight.accel;
      double s1= (pow(straight.vel_end,2.0)- pow(straight.vel_start,2.0))/(2*straight.accel);
      double s2= straight.distance-s1;

        if((0<=s2)  && (s2<=DIST_TOL)) //DIST_TOL need to be changed accordingly to a larger value
        { // v_end= turn_v is achieved at the turn begining
            straight.time=t1; 
            straight.t_breaks.push_back(0); straight.t_breaks.push_back(t1);
            std::vector <double> tcof;
            tcof.push_back(0.5*straight.accel); tcof.push_back(straight.vel_start); tcof.push_back(0);  
            straight.t_coeffs.push_back(tcof);
            std::vector <double> vcof;
            vcof.push_back(straight.accel); vcof.push_back(straight.vel_start);
            straight.vel_coeffs.push_back(vcof);
            straight.s_breaks.push_back(0); straight.s_breaks.push_back(s1); 
        }
        else if (s2>DIST_TOL)
        {// v_end= turn_v is achieved before the turn begining, when accel=max_accel
          double t2=s2/straight.vel_end;
          straight.time=t1+t2;
          straight.t_breaks.push_back(0); straight.t_breaks.push_back(t1); straight.t_breaks.push_back(t1+t2);
          
          std::vector <double> tcof1;
          tcof1.push_back(0.5*straight.accel); tcof1.push_back(straight.vel_start); tcof1.push_back(0);  
          straight.t_coeffs.push_back(tcof1);
          std::vector <double> tcof2;          
          tcof2.push_back(0); tcof2.push_back(straight.vel_end); tcof2.push_back(s1);  
          straight.t_coeffs.push_back(tcof2);

          std::vector <double> vcof1;
          vcof1.push_back(straight.accel); vcof1.push_back(straight.vel_start);
          straight.vel_coeffs.push_back(vcof1);
          std::vector <double> vcof2;
          vcof2.push_back(0); vcof2.push_back(straight.vel_end);
          straight.vel_coeffs.push_back(vcof2);

          straight.s_breaks.push_back(0); straight.s_breaks.push_back(s1); straight.s_breaks.push_back(s1+s2); 
        }
        else
        {
          /* Do nothing*/
        }
      }
    else if(straight.vel_end<straight.vel_start) // Go on decelaration 
      {
        straight.accel= new_accel; // first consider vel_low=vel_end at (distance =straight.start- straight.end) as we don't know the min decelerationS
        double t1 =   (straight.vel_end-straight.vel_start)/straight.accel;
        straight.time=t1;
        straight.t_breaks.push_back(0); straight.t_breaks.push_back(t1);
        std::vector <double> tcof;
        tcof.push_back(0.5*straight.accel); tcof.push_back(straight.vel_start); tcof.push_back(0);  
        straight.t_coeffs.push_back(tcof);
        std::vector <double> vcof;
        vcof.push_back(straight.accel); vcof.push_back(straight.vel_start);
        straight.vel_coeffs.push_back(vcof);
        straight.s_breaks.push_back(0); straight.s_breaks.push_back(straight.distance); 
        
      }
    else    //straight.vel_end = straight.vel_start
    {
        double t1 = straight.distance/straight.vel_end;
        straight.time=t1;
        straight.t_breaks.push_back(0); straight.t_breaks.push_back(t1);
        std::vector <double> tcof;
        tcof.push_back(0); tcof.push_back(straight.vel_end); tcof.push_back(0);  
        straight.t_coeffs.push_back(tcof);
        std::vector <double> vcof;
        vcof.push_back(0); vcof.push_back(straight.vel_end);
        straight.vel_coeffs.push_back(vcof);
        straight.s_breaks.push_back(0); straight.s_breaks.push_back(straight.distance);
    }  
  }
else{
  std::cout<<"The required acceleation is higher than the max accelration for the straight segment "<<straight_i<<std::endl; 
  straight.empty=true;
}
  return straight;
}

void calculate_route_coefficients(route_struct& route){

  route.s_breaks.push_back(0);
  route.t_breaks.push_back(0);
  double prev_s=0; double prev_t=0;
  for (int i=0; i <route.num_straight_lines;i++)
  {
    //////////////////////////////////////FOR ANGLE COEFFICIENTS
   // double heading= atan2((route.wpt_list.at(k+1).y-route.wpt_list.at(k).y), (route.wpt_list.at(k+1).x-route.wpt_list.at(k).x));
    std::vector <double> psi_cof;
    for (int j = 0; j < PSI_COEFF_LENGTH-1; j++) psi_cof.push_back(0); psi_cof.push_back(route.straights.at(i).heading);
    route.psi_coeffs.push_back(psi_cof);

    std::vector <double> curv_cof;
    for (int j = 0; j < CURV_COEFF_LENGTH; j++) curv_cof.push_back(0);
    route.curv_coeffs.push_back(curv_cof);

    for (int j = 0; j < route.straights.at(i).vel_coeffs.size(); j++) route.vel_coeffs.push_back(route.straights.at(i).vel_coeffs.at(j));
    // adding the previous s for the straight begins (t_coeffs ends)
    for (int k = 0; k < route.straights.at(i).t_coeffs.size(); k++){
      route.straights.at(i).t_coeffs.at(k).back()+=prev_s;
      route.t_coeffs.push_back(route.straights.at(i).t_coeffs.at(k));     
    }
      // For straight lines s breaks only used for curve coefficients.. not need to add two sections for acceleration changes as curve is 0   
        route.s_breaks.push_back(route.straights.at(i).s_breaks.back()+prev_s);  

        for (int k = 0; k < route.straights.at(i).t_breaks.size(); k++) {
          route.straights.at(i).t_breaks.at(k)+=prev_t;
          if (k>0) route.t_breaks.push_back(route.straights.at(i).t_breaks.at(k));
        }
    prev_t  =  route.t_breaks.back();
    prev_s = route.s_breaks.back();

     // Integrate turns
    if(i<(route.num_straight_lines-1)){
      if (route.turns.at(i).Sf_total>0){ // Only the real curves (|heading_dif|>0) taking into account

      for (int j = 0; j < route.turns.at(i).curv_coeffs.size(); j++) route.curv_coeffs.push_back(route.turns.at(i).curv_coeffs.at(j));        

      for (int j = 0; j < route.turns.at(i).psi_coeffs.size(); j++) route.psi_coeffs.push_back(route.turns.at(i).psi_coeffs.at(j));
       
      for (int j = 0; j < route.turns.at(i).curv_coeffs.size(); j++) route.vel_coeffs.push_back(route.turns.at(i).vel_poly);///???
         
      //The above need to be adjusted as: vel_coeffs.at(j)
      for (int k = 0; k < route.turns.at(i).t_breaks.size(); k++) {
        route.turns.at(i).t_breaks.at(k)+=prev_t;
        if (k>0) route.t_breaks.push_back(route.turns.at(i).t_breaks.at(k));
      }  
      for (int k = 0; k < route.turns.at(i).t_coeffs.size(); k++){
        route.turns.at(i).t_coeffs.at(k).back()+=prev_s;
        route.t_coeffs.push_back(route.turns.at(i).t_coeffs.at(k));     
      }
      for (int k = 0; k < route.turns.at(i).s_breaks.size(); k++){
        route.turns.at(i).s_breaks.at(k) +=prev_s;
        if (k>0) route.s_breaks.push_back(route.turns.at(i).s_breaks.at(k));     
      }
    prev_t  =  route.t_breaks.back();
    prev_s = route.s_breaks.back();
    }
    }
  }
  route.total_time= route.t_breaks.back();
}

void calculate_complete_path(route_struct& route){
// get heading function  
// No need to calculate route_s_breaks_and_psi_coffs altogether as they were calcuated seperately as follows.
int tot_t_samples = (int)(ceil(route.total_time/TIME_STEP)); 
point_xyz pt_p; pt_p.x=route.wpt_list.at(0).x; pt_p.y=route.wpt_list.at(0).y; pt_p.z=route.wpt_list.at(0).z;
route.path.path_xyz.push_back(pt_p);
for (int i=0; i<tot_t_samples; i++){
  double ith_t_sample=TIME_STEP*i; //double x_integral=0.0; double y_integral=0.0;
  route.path.path_times.push_back(ith_t_sample);
  double ith_s_sample= ppval(route.t_coeffs,route.t_breaks,ith_t_sample);
  double ith_heading=ppval(route.psi_coeffs,route.s_breaks,ith_s_sample);
  ith_heading= (ith_heading<0)? (ith_heading+=2*M_PI):ith_heading;  
  route.path.path_heading.push_back(ith_heading);
  route.path.path_velocity.push_back(ppval(route.vel_coeffs,route.t_breaks,ith_t_sample));
  route.path.path_curv.push_back(ppval(route.curv_coeffs,route.s_breaks,ith_s_sample));
 // route.path.path_xyz.push_back(pt_p);
  pt_p.x +=   (double)TIME_STEP*route.path.path_velocity.back()*cos( ith_heading) ;  //+route.wpt_list.at(0).x; + x_direction of wind_speed* time(s);
  pt_p.y +=   (double)TIME_STEP*route.path.path_velocity.back()*sin(ith_heading) ; //+route.wpt_list.at(0).y; //+ y_direction of wind_speed* time(s);
  route.path.path_xyz.push_back(pt_p);
  route.path.path_bank.push_back( atan2(((pow(route.path.path_velocity.back(),2.0))*route.path.path_curv.back()),GRAVITY));
}
std::vector <double> path_bank_diff=route.path.path_bank; std::vector <double> path_times_diff=route.path.path_times; 
std::adjacent_difference(path_bank_diff.begin(), path_bank_diff.end(), path_bank_diff.begin());
std::adjacent_difference( path_times_diff.begin(),  path_times_diff.end(), path_times_diff.begin());
path_bank_diff.erase(path_bank_diff.begin()); path_times_diff.erase(path_times_diff.begin());
route.path.path_roll_rate=vector_divide(path_bank_diff,path_times_diff);
  }

void to_json(nlohmann::json& j, const mission_item_struct& mi){
j["AMSLAltAboveTerrain"]=mi.AMSLAltAboveTerrain;
j["Altitude"]=mi.Altitude;
j["AltitudeMode"]=mi.AltitudeMode;
j["autoContinue"]=mi.autoContinue;
j["command"]=mi.command;
j["doJumpId"]=mi.doJumpId;
j["frame"]=mi.frame;
j["params"]=mi.params;
j["type"]= mi.type;
}

void to_json(nlohmann::json& j, const mission_struct& m)
{   nlohmann::json nj=m.items;
    j = nlohmann::json{
            { "cruiseSpeed", m.cruiseSpeed},
            { "firmwareType", m.firmwareType},
            { "hoverSpeed", m.hoverSpeed},
            { "items", nj},
            { "plannedHomePosition", {m.plannedHomePosition.Lat, m.plannedHomePosition.Lon, m.plannedHomePosition.Alt}},
           // { "plannedHomePosition", m.plannedHomePosition},
            { "vehicleType", m.vehicleType},
            { "version", m.version}       
          };
}

bool jexists(const nlohmann::json& j, const std::string& key)
{
    return j.find(key) != j.end();
}

struct mission_struct parse_plan_file(const nlohmann::json& jsonfile){
  mission_struct mission;
  mission.cruiseSpeed=jsonfile["mission"]["cruiseSpeed"].get<double>();
  mission.firmwareType=jsonfile["mission"]["firmwareType"].get<int>();
  mission.hoverSpeed=jsonfile["mission"]["hoverSpeed"].get<double>();

  mission.plannedHomePosition.Lat=jsonfile["mission"]["plannedHomePosition"][0].get<double>();
  mission.plannedHomePosition.Lon=jsonfile["mission"]["plannedHomePosition"][1].get<double>();
  mission.plannedHomePosition.Alt=jsonfile["mission"]["plannedHomePosition"][2].get<double>(); 

  mission.vehicleType=jsonfile["mission"]["vehicleType"].get<int>();
  mission.version=jsonfile["mission"]["version"].get<int>();


  for (int i=0; i<jsonfile["mission"]["items"].size(); i++){
    mission_item_struct mission_item;

    if (jexists(jsonfile["mission"]["items"][i],"AMSLAltAboveTerrain")){
      if (jsonfile["mission"]["items"][i]["AMSLAltAboveTerrain"].is_null()) mission_item.AMSLAltAboveTerrain=NAN;
      else mission_item.AMSLAltAboveTerrain=jsonfile["mission"]["items"][i]["AMSLAltAboveTerrain"].get<double>(); }     
    else  mission_item.AMSLAltAboveTerrain=NAN;


    if (jexists(jsonfile["mission"]["items"][i],"Altitude"))
      mission_item.Altitude=jsonfile["mission"]["items"][i]["Altitude"].get<double>();
    else  mission_item.Altitude=NAN;
     
    if (jexists(jsonfile["mission"]["items"][i],"AltitudeMode"))
      mission_item.AltitudeMode=jsonfile["mission"]["items"][i]["AltitudeMode"].get<double>();
     // else  mission_item.AltitudeMode=NAN;
    else  mission_item.AltitudeMode=0; // Keep current altitude
      
    mission_item.autoContinue=jsonfile["mission"]["items"][i]["autoContinue"].get<bool>();
    mission_item.command = jsonfile["mission"]["items"][i]["command"].get<int>();
    mission_item.doJumpId = jsonfile["mission"]["items"][i]["doJumpId"].get<int>();
    mission_item.frame = jsonfile["mission"]["items"][i]["frame"].get<int>(); 
       
    for (int j=0; j<jsonfile["mission"]["items"][i]["params"].size(); j++){
      if (jsonfile["mission"]["items"][i]["params"][j].empty()) mission_item.params.push_back(NAN); 
      else mission_item.params.push_back(jsonfile["mission"]["items"][i]["params"][j].get<double>());
      } 
      mission_item.type=jsonfile["mission"]["items"][i]["type"].get<std::string>();     
     // std::cout<< mission_item.type<<std::endl;  //print_vec_float(mission_item.params);
     
    mission.items.push_back(mission_item);
    }   

  return mission;
}

struct mission_struct process_mission_items(const mission_struct& mission_in){

 mission_struct mission_out=mission_in;
return mission_out;  
}

std::vector <point_LatLonAltVel> Generate_plan_WPs_LatLonAlt(const mission_struct& mission_in){  
 // Process mission items and add new mission items depends on the vehicle constraints 
   //AltitudeMode == 2 --> absolute height above mean sea level, ltitudeMode == 1 --> Altitude sepecified relative to the home (launch) altitude 
  std::vector <point_LatLonAltVel> planWPs; 
  //planWPs.push_back(mission_in.plannedHomePosition); // Add the launch position as the 1st point

  for (int i=0; i < mission_in.items.size();i++){
    if (mission_in.items.at(i).command != 178) {  // Do not count speed change mission items as WPs 
    point_LatLonAltVel pt;
    if(mission_in.items.at(i).AltitudeMode == 2) pt.Alt= mission_in.items.at(i).params.at(6); // abs Altitude 
    else if(mission_in.items.at(i).AltitudeMode == 1) pt.Alt= mission_in.items.at(i).params.at(6)+mission_in.plannedHomePosition.Alt; // relative Altitude  
    else { }
    pt.Lat= mission_in.items.at(i).params.at(4);
    pt.Lon= mission_in.items.at(i).params.at(5);  

    if (mission_in.items.at(i).command == 16) { // Navigate to waypoint CMD
      if(mission_in.items.at(i+1).command == 178) pt.Vel=mission_in.items.at(i+1).params.at(1);
      else pt.Vel=planWPs.back().Vel;
  }
    planWPs.push_back(pt);
    }
  }
  return planWPs;
}


std::vector <point_xyzvdldrfo> Generate_plan_WPs_XYZ(const std::vector <point_LatLonAltVel>& plan_wps_LLA){  
// Assume the home location as (0,0,0)
point_xyzvdldrfo home_pt;
std::vector <point_xyzvdldrfo>  xyz_wps;
 xyz_wps.push_back(home_pt);
for (int i=1;i<plan_wps_LLA.size();i++){

  point_xyzvdldrfo pt;
  double dist_to_wpi= CalcGPSDistance(plan_wps_LLA.at(0).Lat, plan_wps_LLA.at(0).Lon,plan_wps_LLA.at(i).Lat, plan_wps_LLA.at(i).Lon);
  double bearing_to_wpi= -CalcGPSBearing(plan_wps_LLA.at(0).Lat, plan_wps_LLA.at(0).Lon,plan_wps_LLA.at(i).Lat, plan_wps_LLA.at(i).Lon);
  pt.x=dist_to_wpi*cos(bearing_to_wpi);
  pt.y=dist_to_wpi*sin(bearing_to_wpi);
  pt.z=plan_wps_LLA.at(i).Alt;
  pt.v=plan_wps_LLA.at(i).Vel;
  xyz_wps.push_back(pt);
}
return xyz_wps;
}

std::vector<double> convert_wps_to_wp_vector(const std::vector <point_xyzvdldrfo>& wps ){
  const double WPT_DIFF_THRESH=10.0;
  point_xyzvdldrfo pt;
  std::vector<double> wp_vec;
  for (int i=0;i<wps.size();i++){
     pt=wps.at(i); 
    if (i<wps.size()-1) {
     if (sqrt(pow((wps.at(i).x-wps.at(i+1).x),2.0) + pow((wps.at(i).y-wps.at(i+1).y),2.0)) <WPT_DIFF_THRESH) {
       pt=wps.at(i+1);
        i++; 
        }           
    }
    std::cout<<"way points: "<<pt.x<<", "<< pt.y<<", "<< pt.z<<", "<< pt.v<<", "<< pt.dl<<", "<< pt.dr<<", "<< pt.fowp<<std::endl; 
    wp_vec.push_back(pt.x); wp_vec.push_back(pt.y); wp_vec.push_back(pt.z); wp_vec.push_back(pt.v); wp_vec.push_back(pt.dl); wp_vec.push_back(pt.dr); wp_vec.push_back(pt.fowp);
  }
  return wp_vec;
}

std::vector<double> convert_mission_to_wp_vector(const mission_struct& mission_in){ 
 return convert_wps_to_wp_vector(Generate_plan_WPs_XYZ(Generate_plan_WPs_LatLonAlt(mission_in)));
 }

double CalcGPSDistance(const double& lat1d, const double& long1d, const double& lat2d, const double& long2d){
  // Convert to radians
  double lat1r=lat1d*M_PI/180.0; double long1r=long1d*M_PI/180.0; double lat2r=lat2d*M_PI/180.0; double long2r=long2d*M_PI/180.0;
  // Use Harversine formula
  double haversine = (pow(sin((1.0 / 2.0) * (lat2r - lat1r)), 2.0)) + ((cos(lat1r))*(cos(lat2r))*(pow(sin((1.0 / 2.0)*(long2r-long1r)), 2.0)));
  return EARTH_RADIUS* 2.0 * asin(std::min(1.0, sqrt(haversine)));
}

double CalcGPSBearing(const double& lat1d, const double& long1d, const double& lat2d, const double& long2d){
  // Convert to radians
  double lat1r=lat1d*M_PI/180.0; double long1r=long1d*M_PI/180.0; double lat2r=lat2d*M_PI/180.0; double long2r=long2d*M_PI/180.0; 

  return atan2(cos(lat2r)*sin(long2r-long1r),cos(lat1r)*sin(lat2r)- sin(lat1r)*cos(lat2r)*cos(long2r-long1r));
}

