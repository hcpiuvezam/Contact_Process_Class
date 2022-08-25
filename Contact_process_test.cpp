#include "Contact_process_class.hpp"
#include <fstream>
#include <sstream>

using namespace std;

int main(){
    int s = 14;
    double lambda = 1.0,
            t = 0.0;

    Contact_Process CP_object(s, lambda, 0, 4);

    cout << CP_object.t << '\n';
    
    for(int i = 0; i < CP_object.L; ++i){
            CP_object.cells[i] = 1;
            CP_object.active_indexes.push_back(i);
    }
    CP_object.active_num = CP_object.L;
    cout << CP_object.active_num << '\n';

    do{
        CP_object.advance_time();
        CP_object.simulation();
        cout << CP_object.t << ' ' << CP_object.s <<'\n';
    }while(CP_object.t < 100.0);

    return(0);
}