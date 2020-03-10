#include "analytic_evaluator.hpp"
#include <string>

AnalyticEvaluator::AnalyticEvaluator(string mesh_filename, 
        void (*func)(Constants, Point2, int, vector<Point3>&)): _func(func) {
    ifstream mesh(mesh_filename.c_str());

    _mesh =  new GpMesh(); 
    _mesh->setup(mesh, 0, Point3(1.));
}

AnalyticEvaluator::~AnalyticEvaluator(){
}
AnalyticFunction AnalyticEvaluator::select_function(int analytic_func){
    if(analytic_func == 0)
        return &sin3;
    else if(analytic_func == 1)
        return &sin_plus_sin;
    else if(analytic_func == 2)
        return &exp_sin_cos;
    else if(analytic_func == 3)
        return &torus;
    else if(analytic_func == 4)
        return &gaussian;
}

vector<Point3> AnalyticEvaluator::evaluate(int , int evaluation_flags, Point2 uv_coords){
    // patch_id is unused; assumes that we're evaluating a predefined analytic
    // function on a single patch
    int return_value_size = 0;
    if(evaluation_flags & EVAL_VALUE){
        return_value_size += 1;
    } if (evaluation_flags & EVAL_1ST_DERIV){
        return_value_size += 2;
    } if (evaluation_flags & EVAL_2ND_DERIV){
        return_value_size += 3;
    } 
    
    vector<Point3> positions_and_derivs(return_value_size, Point3(0.));
    Constants c;
    c.mean = Point2(0, 0);
    c.std_dev = Point2(1,1)*.075;
    c.amp = .5;
    //c.std_dev = Point2(1,1)*.05;

    _func(c, uv_coords,
         evaluation_flags,
         positions_and_derivs);
    return positions_and_derivs;
}
void AnalyticEvaluator::setup(){
    assert(0);
}

void AnalyticEvaluator::bounding_box(Point3& min, Point3& max){
    _mesh->bbox(min, max);
}
void gaussian(Constants constants, Point2 uv, int flags, vector<Point3>& values){
    Point2 mean = constants.mean;
    Point2 std_dev = constants.std_dev;
    double amp = constants.amp;
    double x = 1*uv.x()-.5;
    double y = 1*uv.y()-.5;

    double mean_x = mean.x();
    double mean_y = mean.y();

    double std_dev_x2 = std_dev.x()*std_dev.x();
    double std_dev_y2 = std_dev.y()*std_dev.y();
    if(flags & EVAL_VALUE){
        assert(values.size() >= 1);
        Point3 function_value(x, y,
                amp*exp(-(
                        (x - mean_x)*(x - mean_x)/(2.*std_dev_x2)
                        + (y - mean_y)*(y - mean_y)/(2.*std_dev_y2)
                     )
                   ));
        values[0] = function_value;

    } 
    if(flags & EVAL_1ST_DERIV){
        assert(flags & EVAL_VALUE);
        Point3 df_du =  Point3(1, 0,  // TODO check this is correct
                -(x - mean_x)/std_dev_x2*values[0].z());
        
        Point3 df_dv =  Point3(0, 1,  // TODO check this is correct
                -(y - mean_y)/std_dev_y2*values[0].z());

        values[1] = df_du; 
        values[2] = df_dv; 

    } 
    if(flags & EVAL_2ND_DERIV){
        assert(flags & EVAL_VALUE);
        Point3 d2f_du2(0, 0,  // TODO check this is correct
                //1./std_dev_x2*(1./(std_dev_x2)*(x - mean_x)*(x - mean_x) -1.)*values[0].z());
                1./pow(std_dev_x2,2)*pow(x - mean_x, 2)*values[0].z() - values[0].z()/std_dev_x2);
        
        Point3 d2f_dudv(0, 0,  // TODO check this is correct
                (y - mean_y)*(x - mean_x)/(std_dev_x2*std_dev_y2)*values[0].z());

        Point3 d2f_dv2(0, 0,  // TODO check this is correct
                //1./std_dev_y2*(1./(std_dev_y2)*(y - mean_y)*(y - mean_y) -1.)*values[0].z());
                1./pow(std_dev_y2,2)*pow(y - mean_y, 2)*values[0].z() - values[0].z()/std_dev_y2);
                //((y - mean_y)*(y - mean_y)/(std_dev_y2) -1.)*values[0].z()/std_dev_y2);

        values[3] = d2f_du2; 
        values[4] = d2f_dudv; 
        values[5] = d2f_dv2; 

    } 
}

void sin3(Constants constants, Point2 uv, int flags, vector<Point3>& values){
    Point2 mean = constants.mean;
    Point2 std_dev = constants.std_dev;
    double x = 2*uv.x()-1.;
    double y = 2*uv.y()-1.;

    double amp = .5;
    double freq_x = 10./(2*M_PI);
    double freq_y = 4/(2*M_PI);
    if(flags & EVAL_VALUE){
        assert(values.size() >= 1);
        Point3 function_value(x, y, amp*pow(sin(freq_x*x),2)*pow(cos(freq_y*y), 2));
        values[0] = function_value;

    } 
    if(flags & EVAL_1ST_DERIV){
        assert(flags & EVAL_VALUE);
        Point3 df_du =  Point3(1, 0, amp*freq_x*2.*sin(x*freq_x)*cos(x*freq_x)*pow(cos(y*freq_y), 2));  // TODO check this is correct
        
        Point3 df_dv =  Point3(0, 1, -amp*freq_y*2.*pow(sin(x*freq_x),2)*sin(freq_y*y)*cos(freq_y*y));  // TODO check this is correct

        values[1] = df_du; 
        values[2] = df_dv; 

    } 
    if(flags & EVAL_2ND_DERIV){
        assert(flags & EVAL_VALUE);
        Point3 d2f_du2(0, 0, 2*amp*freq_x*freq_x*cos(2*freq_x*x)*pow(cos(freq_y*y),2)); // TODO check this is correct
        
        Point3 d2f_dudv(0, 0, -4*amp*freq_x*freq_y*sin(freq_x*x)*cos(freq_x*x)*sin(freq_y*y)*cos(freq_y*y));  // TODO check this is correct

        Point3 d2f_dv2(0, 0, -2*amp*freq_y*freq_y*cos(2*freq_y*y)*pow(sin(freq_x*x),2));  // TODO check this is correct
                //1./std_dev_y2*(1./(std_dev_y2)*(y - mean_y)*(y - mean_y) -1.)*values[0].z());
                //((y - mean_y)*(y - mean_y)/(std_dev_y2) -1.)*values[0].z()/std_dev_y2);
        values[3] = d2f_du2; 
        values[4] = d2f_dudv; 
        values[5] = d2f_dv2; 

    } 
}

void sin_plus_sin(Constants constants, Point2 uv, int flags, vector<Point3>& values){
    Point2 mean = constants.mean;
    Point2 std_dev = constants.std_dev;
    double x = 2*uv.x()-1.;
    double y = 2*uv.y()-1.;

    double amp = .5;
    double freq_x = 1.;
    double freq_y = 10.;
    if(flags & EVAL_VALUE){
        assert(values.size() >= 1);
        Point3 function_value(x, y, amp*(pow(sin(freq_x*x), 2) + sin(freq_y*y)) );
        values[0] = function_value;

    } 
    if(flags & EVAL_1ST_DERIV){
        assert(flags & EVAL_VALUE);
        Point3 df_du =  Point3(1, 0,amp*2.*freq_x*sin(freq_x*x)*cos(freq_x*x) );  // TODO check this is correct
        
        Point3 df_dv =  Point3(0, 1, amp*freq_y*cos(freq_y*y));  // TODO check this is correct

        values[1] = df_du; 
        values[2] = df_dv; 

    } 
    if(flags & EVAL_2ND_DERIV){
        assert(flags & EVAL_VALUE);
        Point3 d2f_du2(0, 0, 2*amp*freq_x*freq_x*cos(2*freq_x*x)); // TODO check this is correct
        
        Point3 d2f_dudv(0, 0, 0);  // TODO check this is correct

        Point3 d2f_dv2(0, 0, -amp*freq_y*freq_y*sin(freq_y*y));  // TODO check this is correct

        values[3] = d2f_du2; 
        values[4] = d2f_dudv; 
        values[5] = d2f_dv2; 

    } 
}

void exp_sin_cos(Constants constants, Point2 uv, int flags, vector<Point3>& values){
    Point2 mean = constants.mean;
    Point2 std_dev = constants.std_dev;
    double x = 2*uv.x()-1.;
    double y = 2*uv.y()-1.;

    double amp = .1;
    double freq_x = 2*M_PI;
    double freq_y = 2*M_PI;
    if(flags & EVAL_VALUE){
        assert(values.size() >= 1);
        Point3 function_value(x, y, amp*exp(sin(freq_x*x)*cos(freq_y*y)) );
        values[0] = function_value;

    } 
    if(flags & EVAL_1ST_DERIV){
        assert(flags & EVAL_VALUE);
        Point3 df_du =  Point3(1, 0, amp*freq_x*cos(freq_x*x)*cos(freq_y*y)*exp(sin(freq_x*x)*cos(freq_y*y)) );  // TODO check this is correct
        
        Point3 df_dv =  Point3(0, 1, -amp*freq_y*sin(freq_x*x)*sin(freq_y*y)*exp(sin(freq_x*x)*cos(freq_y*y)));  // TODO check this is correct

        values[1] = df_du; 
        values[2] = df_dv; 

    } 
    if(flags & EVAL_2ND_DERIV){
        assert(flags & EVAL_VALUE);
        Point3 d2f_du2(0, 0, 
               amp*exp(sin(freq_x*x)*cos(freq_y*y))*(
                   freq_x*freq_x*pow(cos(freq_x*x),2)*pow(cos(freq_y*y),2) - 
                   freq_x*freq_x*sin(freq_x*x)*cos(freq_y*y)
                       )
                ); // TODO check this is correct
        
        Point3 d2f_dudv(0, 0, 
               -amp*freq_x*freq_y*cos(freq_x*x)*sin(freq_y*y)*
               exp(sin(freq_x*x)*cos(freq_y*y))*(
                   sin(freq_x*x)*cos(freq_y*y) + 1
                       )
                );  // TODO check this is correct

        Point3 d2f_dv2(0, 0, 
               amp*exp(sin(freq_x*x)*cos(freq_y*y))*(
                   freq_y*freq_y*pow(sin(freq_x*x),2)*pow(sin(freq_y*y),2) - 
                   freq_y*freq_y*sin(freq_x*x)*cos(freq_y*y)
                       )
                    );  // TODO check this is correct

        values[3] = d2f_du2; 
        values[4] = d2f_dudv; 
        values[5] = d2f_dv2; 

    } 
}

void torus(Constants constants, Point2 uv, int flags, vector<Point3>& values){
    Point2 mean = constants.mean;
    Point2 std_dev = constants.std_dev;
    //double x = 2*uv.x()-1.;
    //double y = 2*uv.y()-1.;
    double x = uv.x();
    double y = uv.y();

    double freq_x = 2*M_PI;
    double freq_y = 2*M_PI;
    double radius = .5;
    Point3 center(1.,1.,0.);
    if(flags & EVAL_VALUE){
        assert(values.size() >= 1);
        Point3 function_value( 
                (center.x() + radius*cos(freq_y*y))*cos(freq_x*x),
                (center.y() + radius*cos(freq_y*y))*sin(freq_x*x),
                radius*sin(freq_y*y));
        values[0] = function_value;

    } 
    if(flags & EVAL_1ST_DERIV){
        assert(flags & EVAL_VALUE);
        Point3 df_du =  Point3( 
                (center.x() + radius*cos(freq_y*y))*(-freq_x*sin(freq_x*x)),
                (center.y() + radius*cos(freq_y*y))*(freq_x*cos(freq_x*x)),
                0.);
        
        Point3 df_dv =  Point3(
            radius*cos(freq_x*x)*(-freq_y*sin(freq_y*y)),
            radius*sin(freq_x*x)*(-freq_y*sin(freq_y*y)),
            freq_y*radius*cos(freq_y*y));

        values[1] = df_du; 
        values[2] = df_dv; 

    } 
    if(flags & EVAL_2ND_DERIV){
        assert(flags & EVAL_VALUE);
        assert(0);
        Point3 d2f_du2(0, 0, 0);
        
        Point3 d2f_dudv(0, 0, 0);

        Point3 d2f_dv2(0, 0, 0); 

        values[3] = d2f_du2; 
        values[4] = d2f_dudv; 
        values[5] = d2f_dv2; 

    } 
}




