#include "diff_geom.hpp"


vector<double> Differential::first_fundamental_form(Function X_u, Function X_v, Point2 uv_coords){
    NumVector uv(2, false, uv_coords.array()); 
    // evaluate the partial derivatives at the (u,v) coordinates
    NumVector X_u_vec = X_u(uv);
    NumVector X_v_vec = X_v(uv);
    
    assert(X_u_vec.length() == 3);
    assert(X_v_vec.length() == 3);

    Point3 X_u_at_uv(X_u_vec.data());
    Point3 X_v_at_uv(X_v_vec.data());

    vector<double> form_coefficients(3,0);
    // only return E, F, G since matrix 
    // | E F | 
    // | F G |
    // is symmetric
    
    form_coefficients[0] = dot(X_u_at_uv, X_u_at_uv); // E 
    form_coefficients[1] = dot(X_u_at_uv, X_v_at_uv); // F
    form_coefficients[2] = dot(X_v_at_uv, X_v_at_uv); // G
    return form_coefficients;
}

vector<double> Differential::first_fundamental_form(Point3 X_u, Point3 X_v){

    vector<double> form_coefficients(3,0);
    form_coefficients[0] = dot(X_u, X_u); // E 
    form_coefficients[1] = dot(X_u, X_v); // F
    form_coefficients[2] = dot(X_v, X_v); // G

    // only return E, F, G since matrix 
    // | E F | 
    // | F G |
    // is symmetric
    assert(form_coefficients.size() == 3);
    return form_coefficients;
}

vector<double> Differential::second_fundamental_form(Function X_u, Function X_v, 
        Function X_uu, Function X_uv, Function X_vv,
        Point2 uv_coords){
    NumVector uv(2, false, uv_coords.array()); 
    // evaluate the first and second derivatives at the (u,v) coordinates
    NumVector X_u_vec = X_u(uv);
    NumVector X_v_vec = X_v(uv);
    NumVector X_uu_vec = X_uu(uv);
    NumVector X_uv_vec = X_uv(uv);
    NumVector X_vv_vec = X_vv(uv);
    
    assert(X_u_vec.length() == 3);
    assert(X_v_vec.length() == 3);
    assert(X_uu_vec.length() == 3);
    assert(X_uv_vec.length() == 3);
    assert(X_vv_vec.length() == 3);

    Point3 X_u_at_uv(X_u_vec.data());
    Point3 X_v_at_uv(X_v_vec.data());
    Point3 X_uu_at_uv(X_uu_vec.data());
    Point3 X_uv_at_uv(X_uv_vec.data());
    Point3 X_vv_at_uv(X_vv_vec.data());

    Point3 normal = cross(X_u_at_uv, X_v_at_uv).dir();

    vector<double> form_coefficients(3,0);
    // only return E, F, G since matrix 
    // | E F | 
    // | F G |
    // is symmetric
    
    form_coefficients[0] = dot(X_uu_at_uv, normal); // E 
    form_coefficients[1] = dot(X_uv_at_uv, normal); // F
    form_coefficients[2] = dot(X_vv_at_uv, normal); // G
    return form_coefficients;
}


vector<double> Differential::second_fundamental_form(Point3 X_u, Point3 X_v, 
        Point3 X_uu, Point3 X_uv, Point3 X_vv){
    Point3 normal = cross(X_u, X_v).dir();

    vector<double> form_coefficients(3,0);
    form_coefficients[0] = dot(X_uu, normal); // L 
    form_coefficients[1] = dot(X_uv, normal); // M 
    form_coefficients[2] = dot(X_vv, normal); // N

    // only return L, M, N since matrix 
    // | L M | 
    // | M N |
    // is symmetric
    assert(form_coefficients.size() == 3);
    return form_coefficients;
}

vector<double> Differential::first_fundamental_form_partial_derivatives(
        Function X_u, Function X_v, 
        Function X_uu, Function X_uv, Function X_vv,
        Point2 uv_coords){

    NumVector uv(2, false, uv_coords.array()); 
    // evaluate the first and second derivatives at the (u,v) coordinates
    NumVector X_u_vec = X_u(uv);
    NumVector X_v_vec = X_v(uv);
    NumVector X_uu_vec = X_uu(uv);
    NumVector X_uv_vec = X_uv(uv);
    NumVector X_vv_vec = X_vv(uv);
    
    assert(X_u_vec.length() == 3);
    assert(X_v_vec.length() == 3);
    assert(X_uu_vec.length() == 3);
    assert(X_uv_vec.length() == 3);
    assert(X_vv_vec.length() == 3);

    Point3 X_u_at_uv(X_u_vec.data());
    Point3 X_v_at_uv(X_v_vec.data());
    Point3 X_uu_at_uv(X_uu_vec.data());
    Point3 X_uv_at_uv(X_uv_vec.data());
    Point3 X_vv_at_uv(X_vv_vec.data());

    Point3 normal = cross(X_u_at_uv, X_v_at_uv).dir();


    vector<double> form_coeffs_derivatives(6,0);

    form_coeffs_derivatives[0] = 2*dot(X_u_at_uv, X_uu_at_uv); // E_u
    form_coeffs_derivatives[1] = 2*dot(X_u_at_uv, X_uv_at_uv); // E_v

    form_coeffs_derivatives[2] = dot(X_u_at_uv, X_uv_at_uv) + dot(X_uu_at_uv, X_v_at_uv); // F_u
    form_coeffs_derivatives[3] = dot(X_u_at_uv, X_vv_at_uv) + dot(X_uv_at_uv, X_v_at_uv); // F_v

    form_coeffs_derivatives[4] = 2*dot(X_v_at_uv, X_uv_at_uv); // G_u
    form_coeffs_derivatives[5] = 2*dot(X_v_at_uv, X_vv_at_uv); // G_v
    
    assert(form_coeffs_derivatives.size() == 6); 
    return form_coeffs_derivatives;
}

vector<double> Differential::first_fundamental_form_partial_derivatives(
        Point3 X_u, Point3 X_v, Point3 X_uu, Point3 X_uv, Point3 X_vv){
    
    vector<double> form_coeffs_derivatives(6,0);

    form_coeffs_derivatives[0] = 2*dot(X_u, X_uu); // E_u
    form_coeffs_derivatives[1] = 2*dot(X_u, X_uv); // E_v

    form_coeffs_derivatives[2] = dot(X_u, X_uv) + dot(X_uu, X_v); // F_u
    form_coeffs_derivatives[3] = dot(X_u, X_vv) + dot(X_uv, X_v); // F_v

    form_coeffs_derivatives[4] = 2*dot(X_v, X_uv); // G_u
    form_coeffs_derivatives[5] = 2*dot(X_v, X_vv); // G_v
    
    assert(form_coeffs_derivatives.size() == 6); 
    return form_coeffs_derivatives;
}



double Differential::determinant_of_first_fundamental_form(vector<double> first_fundamental_form_coefs){
    assert(first_fundamental_form_coefs.size() == 3);
    double E = first_fundamental_form_coefs[0];
    double F = first_fundamental_form_coefs[1];
    double G = first_fundamental_form_coefs[2];
    return sqrt(E*G-F*F);
}

vector<double> Differential::derivatives_of_det_of_first_form(
        vector<double> first_form_coefs,
        vector<double> first_form_derivatives){
    
    assert(first_form_coefs.size() == 3);
    assert(first_form_derivatives.size() == 6);

    vector<double> derivatives(2,0);
    double E = first_form_coefs[0];
    double F = first_form_coefs[1];
    double G = first_form_coefs[2];
    
    double E_u = first_form_derivatives[0];
    double E_v = first_form_derivatives[1];

    double F_u = first_form_derivatives[2];
    double F_v = first_form_derivatives[3];

    double G_u = first_form_derivatives[4];
    double G_v = first_form_derivatives[5];

    // W_u
    derivatives[0] = .5*pow(E*G - F*F, .5)*(E*G_u + G*E_u  - 2*F*F_u);

    // W_v
    derivatives[1] = .5*pow(E*G - F*F, .5)*(E*G_v + G*E_v  - 2*F*F_v);
    return derivatives;
}

double Differential::gaussian_curvature(Point3 X_u, Point3 X_v, 
        Point3 X_uu, Point3 X_uv, Point3 X_vv){
    vector<double> first_form_coefs = 
        Differential::first_fundamental_form(X_u, X_v);
    vector<double> second_form_coefs = 
        Differential::second_fundamental_form(X_u, X_v, X_uu, X_uv, X_vv);

    //double det_of_first_form = determinant_of_first_fundamental_form(first_form_coefs);
    double E = first_form_coefs[0];
    double F = first_form_coefs[1];
    double G = first_form_coefs[2];

    double L = second_form_coefs[0];
    double M = second_form_coefs[1];
    double N = second_form_coefs[2];
    return  (L*N - M*M)/(E*G - F*F);
}

double Differential::mean_curvature(Function X_u, Function X_v, 
        Function X_uu, Function X_uv, Function X_vv, Point2 uv){
    vector<double> first_form_coefs = 
        Differential::first_fundamental_form(X_u, X_v, uv);
    vector<double> second_form_coefs = 
        Differential::second_fundamental_form(X_u, X_v, X_uu, X_uv, X_vv, uv);

    double E = first_form_coefs[0];
    double F = first_form_coefs[1];
    double G = first_form_coefs[2];

    double L = second_form_coefs[0];
    double M = second_form_coefs[1];
    double N = second_form_coefs[2];
    return .5*(E*N - 2*F*M + G*L)/(E*G-F*F);
}


double Differential::mean_curvature(Point3 X_u, Point3 X_v, 
        Point3 X_uu, Point3 X_uv, Point3 X_vv){
    vector<double> first_form_coefs = 
        Differential::first_fundamental_form(X_u, X_v);
    vector<double> second_form_coefs = 
        Differential::second_fundamental_form(X_u, X_v, X_uu, X_uv, X_vv);

    double E = first_form_coefs[0];
    double F = first_form_coefs[1];
    double G = first_form_coefs[2];

    double L = second_form_coefs[0];
    double M = second_form_coefs[1];
    double N = second_form_coefs[2];
    return .5*(E*N - 2*F*M + G*L)/(E*G-F*F);
}

Point3 Differential::surface_gradient(Function X_u, Function X_v, 
        Function f_u, Function f_v, Point2 uv_coords){
    vector<double> first_form_coefs = 
        Differential::first_fundamental_form(X_u, X_v, uv_coords);
    NumVector uv(2, false, uv_coords.array()); 
    
    // evaluate the first derivatives of the scalar function at the (u,v) coordinates
    NumVector f_u_vec = f_u(uv);
    NumVector f_v_vec = f_v(uv);
    NumVector X_u_vec = X_u(uv);
    NumVector X_v_vec = X_v(uv);
    

    // only take gradient of scalars
    assert(f_u_vec.length() == 1);
    assert(f_v_vec.length() == 1);
    assert(X_u_vec.length() == 3);
    assert(X_v_vec.length() == 3);

    double f_u_at_uv = f_u_vec(0);
    double f_v_at_uv = f_v_vec(0);
    Point3 X_u_at_uv(X_u_vec.data());
    Point3 X_v_at_uv(X_v_vec.data());


    double E = first_form_coefs[0];
    double F = first_form_coefs[1];
    double G = first_form_coefs[2];
    double W = Differential::determinant_of_first_fundamental_form(first_form_coefs);
    return (G*X_u_at_uv - F*X_v_at_uv)/(W*W)*f_u_at_uv 
        + (E*X_v_at_uv - F*X_u_at_uv)/(W*W)*f_v_at_uv;
}

Point3 Differential::surface_gradient(Point3 X_u, Point3 X_v, 
        double f_u, double f_v){
    vector<double> first_form_coefs = 
        Differential::first_fundamental_form(X_u, X_v);

    double E = first_form_coefs[0];
    double F = first_form_coefs[1];
    double G = first_form_coefs[2];
    double W = Differential::determinant_of_first_fundamental_form(first_form_coefs);

    return (G*X_u - F*X_v)/(W*W)*f_u + (E*X_v - F*X_u)/(W*W)*f_v;
}

NumMatrix Differential::surface_gradient(Point3 X_u, Point3 X_v, 
        Point3 V_u, Point3 V_v){
    
    const int dim = 3;
    NumMatrix surface_grad(dim, dim);
    for(int i = 0 ; i < dim; i++){
        Point3 ith_surface_grad = 
            Differential::surface_gradient(X_u, X_v, V_u(i), V_v(i));
        for(int j =0; j < dim; j++){
            surface_grad(i,j) = ith_surface_grad(i);
        }
    }
   return surface_grad;
}

double Differential::surface_divergence(Function X_u, Function X_v, 
        Function V_u, Function V_v, Point2 uv_coords){
    NumVector uv(2, false, uv_coords.array()); 
    
    // evaluate the first derivatives of the scalar function at the (u,v) coordinates
    NumVector V_u_vec = V_u(uv);
    NumVector V_v_vec = V_v(uv);
    NumVector X_u_vec = X_u(uv);
    NumVector X_v_vec = X_v(uv);
    

    // only take gradient of scalars
    assert(V_u_vec.length() == 3);
    assert(V_v_vec.length() == 3);
    assert(X_u_vec.length() == 3);
    assert(X_v_vec.length() == 3);

    Point3 V_u_at_uv(V_u_vec.data());
    Point3 V_v_at_uv(V_v_vec.data());
    Point3 X_u_at_uv(X_u_vec.data());
    Point3 X_v_at_uv(X_v_vec.data());

    vector<double> first_form_coefs = 
        Differential::first_fundamental_form(X_u, X_v, uv_coords);

    double E = first_form_coefs[0];
    double F = first_form_coefs[1];
    double G = first_form_coefs[2];
    double W = Differential::determinant_of_first_fundamental_form(first_form_coefs);

    return dot((G*X_u_at_uv - F*X_v_at_uv)/(W*W), V_u_at_uv) + 
        dot((E*X_v_at_uv - F*X_u_at_uv)/(W*W), V_v_at_uv);

}


double Differential::surface_divergence(Point3 X_u, Point3 X_v, 
        Point3 V_u, Point3 V_v){
    vector<double> first_form_coefs = 
        Differential::first_fundamental_form(X_u, X_v);

    double E = first_form_coefs[0];
    double F = first_form_coefs[1];
    double G = first_form_coefs[2];
    double W = Differential::determinant_of_first_fundamental_form(first_form_coefs);

    return dot((G*X_u - F*X_v)/(W*W), V_u) + dot((E*X_v - F*X_u)/(W*W), V_v);
}

double Differential::laplace_beltrami(
       vector<Function> X_derivatives,
       vector<Function> f_derivatives,
       Point2 uv_coords){
    NumVector uv(2, false, uv_coords.array()); 
    assert(X_derivatives.size() == f_derivatives.size());
    int num_derivs = X_derivatives.size();

    vector<Point3> X_deriv_evalauted(num_derivs, Point3());
    vector<double> f_deriv_evalauted(num_derivs, 0);
    
    // make sure everything has the right domain/range dimensions...
    for(int i = 0; i < num_derivs; i++){
        Function dX_dx = X_derivatives[i];
        Function df_dx = f_derivatives[i];
        
        assert(dX_dx.domain_dim() == 2);
        assert(df_dx.domain_dim() == 2);
        assert(dX_dx.range_dim() == 3);
        assert(df_dx.range_dim() == 1);
    }


    for(int i = 0; i < num_derivs; i++){
        Function dX_dx = X_derivatives[i];
        Function df_dx = f_derivatives[i];

        // evaluate derivatives
        NumVector dX_dx_value = dX_dx(uv);
        NumVector df_dx_value = df_dx(uv);
        
        // store values
        Point3 dX_dx_at_uv(dX_dx_value.data());
        double df_dx_at_uv = df_dx_value(0);
        
        X_deriv_evalauted[i] = dX_dx_at_uv;
        f_deriv_evalauted[i] = df_dx_at_uv;
    }

    Point3 X_u =  X_deriv_evalauted[0];
    Point3 X_v =  X_deriv_evalauted[1];
    Point3 X_uu = X_deriv_evalauted[2];
    Point3 X_uv = X_deriv_evalauted[3];
    Point3 X_vv = X_deriv_evalauted[4];
    
    double f_u =  f_deriv_evalauted[0];
    double f_v =  f_deriv_evalauted[1];
    double f_uu = f_deriv_evalauted[2];
    double f_uv = f_deriv_evalauted[3];
    double f_vv = f_deriv_evalauted[4];


    vector<double> first_form_coefs = 
        Differential::first_fundamental_form(
                X_derivatives[0], 
                X_derivatives[1], 
                uv_coords);

    vector<double> first_form_derivatives = 
        Differential::first_fundamental_form_partial_derivatives(
                X_derivatives[0], 
                X_derivatives[1], 
                X_derivatives[2], 
                X_derivatives[3], 
                X_derivatives[4], 
                uv_coords);

    double E = first_form_coefs[0];
    double F = first_form_coefs[1];
    double G = first_form_coefs[2];
    
    //double E_u = first_form_derivatives[0];
    double E_v = first_form_derivatives[1];

    double F_u = first_form_derivatives[2];
    double F_v = first_form_derivatives[3];

    double G_u = first_form_derivatives[4];
    //double G_v = first_form_derivatives[5];

    double W = 
        Differential::determinant_of_first_fundamental_form(first_form_coefs);
    vector<double> det_of_first_form_derivs = 
        Differential::derivatives_of_det_of_first_form(first_form_coefs, first_form_derivatives);
    
    double W_u = det_of_first_form_derivs[0];
    double W_v = det_of_first_form_derivs[1];


    // see http://www.math.lsa.umich.edu/~shravan/papers/ves3d.pdf
    return 1./W*(
            (G*f_uu + f_u*G_u - F*f_uv - F_u*f_v)/W - W_u/(W*W)*(G*f_u - F*f_v)
          + (E*f_vv + f_v*E_v - F*f_uv - F_v*f_u)/W - W_v/(W*W)*(E*f_v - F*f_u));


    return 0.;
}
double Differential::laplace_beltrami(Point3 X_u, Point3 X_v, 
        Point3 X_uu, Point3 X_uv, Point3 X_vv,
        double f_u, double f_v,
        double f_uu, double f_uv, double f_vv){

    vector<double> first_form_coefs = 
        Differential::first_fundamental_form(X_u, X_v);

    vector<double> first_form_derivatives = 
        Differential::first_fundamental_form_partial_derivatives(
                X_u,  X_v,  X_uu,  X_uv,  X_vv);

    double E = first_form_coefs[0];
    double F = first_form_coefs[1];
    double G = first_form_coefs[2];
    
    //double E_u = first_form_derivatives[0];
    double E_v = first_form_derivatives[1];

    double F_u = first_form_derivatives[2];
    double F_v = first_form_derivatives[3];

    double G_u = first_form_derivatives[4];
    //double G_v = first_form_derivatives[5];

    double W = 
        Differential::determinant_of_first_fundamental_form(first_form_coefs);
    vector<double> det_of_first_form_derivs = 
        Differential::derivatives_of_det_of_first_form(first_form_coefs, first_form_derivatives);
    
    double W_u = det_of_first_form_derivs[0];
    double W_v = det_of_first_form_derivs[1];


    //Point3 surface_grad = Differential::surface_gradient(X_u, X_v, f_u, f_v);
    //return Differential::surface_divergence(X_u, X_v, surface_grad
    return 1./W*(
            (G*f_uu + f_u*G_u - F*f_uv - F_u*f_v)/W - W_u/(W*W)*(G*f_u - F*f_v)
          + (E*f_vv + f_v*E_v - F*f_uv - F_v*f_u)/W - W_v/(W*W)*(E*f_v - F*f_u));
}


Point3 Differential::laplace_beltrami(
        Point3 X_u, Point3 X_v, Point3 X_uu, Point3 X_uv, Point3 X_vv,
        Point3 V_u, Point3 V_v, Point3 V_uu, Point3 V_uv, Point3 V_vv){

    const int dim = 3;
    Point3 laplace_beltrami;
    for(int i = 0 ; i < dim; i++){
        double ith_laplace_beltrami = 
            Differential::laplace_beltrami( X_u,  X_v, X_uu,  X_uv,  X_vv,
                            V_u(i),  V_v(i), V_uu(i),  V_uv(i),  V_vv(i));
        laplace_beltrami(i) = ith_laplace_beltrami;
    }
   return laplace_beltrami;


}




double Differential::gaussian_curvature(Point2 uv, FaceMapSurf* surface, int patch_id){
    assert(patch_id >= 0);
    Point3 position_and_derivs[6];
    surface->eval(EVAL_VALUE|EVAL_1ST_DERIV|EVAL_2ND_DERIV,
            patch_id, uv, position_and_derivs);

    return Differential::gaussian_curvature(
            position_and_derivs[1], // X_u
            position_and_derivs[2], // X_v
            position_and_derivs[3], // X_uu
            position_and_derivs[4], // X_uv
            position_and_derivs[5]); // X_vv

}


double Differential::mean_curvature(Point2 uv, FaceMapSurf* surface, int patch_id){
    assert(patch_id >= 0);
    Point3 position_and_derivs[6];
   surface->eval(EVAL_VALUE|EVAL_1ST_DERIV|EVAL_2ND_DERIV,
            patch_id, uv, position_and_derivs);

    return Differential::mean_curvature(
            position_and_derivs[1], // X_u
            position_and_derivs[2], // X_v
            position_and_derivs[3], // X_uu
            position_and_derivs[4], // X_uv
            position_and_derivs[5]); // X_vv

}

Point3 Differential::surface_gradient(Point2 uv, double f_u, double f_v,
        FaceMapSurf* surface, int patch_id){
    assert(patch_id >= 0);
    Point3 position_and_derivs[3];
   surface->eval(EVAL_VALUE|EVAL_1ST_DERIV,
            patch_id, uv, position_and_derivs);

   return Differential::surface_gradient(
           position_and_derivs[1],
           position_and_derivs[2],
           f_u,
           f_v);
}

double Differential::surface_divergence(Point2 uv, Point3 V_u, Point3 V_v,
        FaceMapSurf* surface, int patch_id){
    assert(patch_id >= 0);
    Point3 position_and_derivs[3];
   surface->eval(EVAL_VALUE|EVAL_1ST_DERIV,
            patch_id, uv, position_and_derivs);
   return Differential::surface_divergence(
           position_and_derivs[1],
           position_and_derivs[2],
           V_u,
           V_v);
}

double Differential::laplace_beltrami(Point2 uv, double f_u, double f_v,
        double f_uu, double f_uv, double f_vv,
        FaceMapSurf* surface, int patch_id){
    assert(patch_id >= 0);
    Point3 position_and_derivs[6];
   surface->eval(EVAL_VALUE|EVAL_1ST_DERIV|EVAL_2ND_DERIV,
            patch_id, uv, position_and_derivs);
   return Differential::laplace_beltrami(
            position_and_derivs[1], // X_u
            position_and_derivs[2], // X_v
            position_and_derivs[3], // X_uu
            position_and_derivs[4], // X_uv
            position_and_derivs[5], // X_vv
            f_u, f_v, f_uu, f_uv, f_vv);
}

Point3 Differential::laplace_beltrami(Point2 uv, Point3 V_u, Point3 V_v,
        Point3 V_uu, Point3 V_uv, Point3 V_vv,
        FaceMapSurf* surface, int patch_id){
    assert(patch_id >= 0);
    Point3 position_and_derivs[6];
   surface->eval(EVAL_VALUE|EVAL_1ST_DERIV|EVAL_2ND_DERIV,
            patch_id, uv, position_and_derivs);
   return Differential::laplace_beltrami(
            position_and_derivs[1], // X_u
            position_and_derivs[2], // X_v
            position_and_derivs[3], // X_uu
            position_and_derivs[4], // X_uv
            position_and_derivs[5], // X_vv
            V_u, V_v, V_uu, V_uv, V_vv);
}




