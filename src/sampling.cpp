#include "sampling.hpp"

double Sampling::affine_coordinate_transform(double u, 
        Interval source_domain, 
        Interval target_domain){
    // source_domain bounds [a,b] 
    double a = source_domain.first;
    double b = source_domain.second;
    
    // target_domain bounds [c,d] 
    double c = target_domain.first;
    double d = target_domain.second;

    
    return c + (d-c)/(b-a)*(u - a);
}
Point2 Sampling::affine_coordinate_transform(Point2 uv, Rectangle source_domain, 
        Rectangle target_domain){
    
    return Point2(
            Sampling::affine_coordinate_transform(uv.x(), 
                source_domain.first, target_domain.first),
            Sampling::affine_coordinate_transform(uv.y(), 
                source_domain.second, target_domain.second)
            );
}

void Sampling::remap_samples(Interval source_domain, Interval target_domain,
        NumVector& samples){
    
    int num_samples = samples.length();
    for(int i = 0; i < num_samples; i++){
        samples(i) = affine_coordinate_transform(
                samples(i),
                source_domain,
                target_domain);
    }
}


//----------------------------------------------------------------------------
// Sampling patterns
//----------------------------------------------------------------------------
// Equispaced
//----------------------------------------------------------------------------

NumVector Sampling::equispaced(int num_samples, Interval domain){
    NumVector samples(num_samples);
    double step_size = 1./double(num_samples-1);
    
    for(int i = 0; i < num_samples; i++){
        samples(i) = i*step_size;
    }
    
    remap_samples(base_domain.first, domain, samples);
    return samples;
}


//----------------------------------------------------------------------------
// Chebyshev
//----------------------------------------------------------------------------

NumVector Sampling::chebyshev1(int num_samples, Interval domain){
    NumVector samples(num_samples);
    double step_size = 1./double(num_samples);
    
    // generate the usual chebyshev grid on [-1,1]
    for(int i = 1; i < num_samples+1; i++){
        samples(i-1) = cos((2*i-1)/2.*M_PI*step_size);
    }

    // move it to [0,1] like we want.
    remap_samples(Interval(-1., 1.), Interval(0.,1.), samples);
    // then move it to the target domain.
    remap_samples(base_domain.first, domain, samples);
    return samples;
}

NumVector Sampling::chebyshev2(int num_samples, Interval domain){
    NumVector samples(num_samples);
    double step_size = 1./double(num_samples-1);
    
    // generate the usual chebyshev grid on [-1,1]
    for(int i = 0; i < num_samples; i++){
        samples(i) = cos(i*M_PI*step_size);
    }

    // move it to [0,1] like we want.
    remap_samples(Interval(-1., 1.), Interval(0.,1.), samples);
    // then move it to the target domain.
    remap_samples(base_domain.first, domain, samples);
    return samples;
}
