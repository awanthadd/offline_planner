
#include <offline_planner/display.hpp>

////For printing the data from service end for testing
void print_vec_float(std::vector<double> const &input)
{
    for (int i = 0; i < input.size(); i++) {
        std::cout << i<<' '<<input.at(i) <<std::endl;
       // std::cout << ' '<<input.at(i);
    }
    std::cout<< std::endl;
}
void print_vec_vec_float(std::vector <std::vector <double> > const &input)
{
    for (int i = 0; i < input.size(); i++) { std::cout<<i<<' ';
      for (int j = 0; j < input[i].size(); j++) {
        std::cout<< input[i][j] << ' '; 
      }
        std::cout<< std::endl;    }
    std::cout<< std::endl;
    }
void print_vec_string(std::vector<std::string> const &input)
{
    for (int i = 0; i < input.size(); i++) {
        std::cout << input.at(i) << ' ';
    }
    std::cout<< std::endl;
}
/////////////////////////////////////////////////////////////////////
