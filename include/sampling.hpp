#ifndef __SAMPLING_HPP__
#define __SAMPLING_HPP__
#include <nummatrix.hpp>
#include <numvector.hpp>
#include "common.hpp"
#include <vec2t.hpp>
#include <vec3t.hpp>
/**
 * Various sampling patterns in a two dimensional parameter space [0,1]^2. 
 * All functions return a 2 x num_samples^2 length vector, indexed by two
 * variables i,j with i*num_samples + j. interval endpoints are always included.
 */
namespace Sampling {
    const int parameter_dim = 2;
    const Rectangle base_domain(Interval(0.,1.),Interval(0.,1.));
    
    /**
     * Suppose source_domain = [a,b] and target_domain = [c,d] 
     * and u \in source_domain. The affine map that sends u to its
     * image in target_domain is
     *
     * u \to  c  + (d-c)/(b-a)*(u - a)
     *
     * assuming source_domain is non-degenerate (a\not=b). This
     * function computes this map
     *
     * Recall that Interval =  pair<double, double>
     *
     * @param   double       u                  coordinate in source_domain to  
     *                                          be mapped into target_domain
     * @param   Rectangle   source_domain       [a,b] containing u
     * @param   Rectangle   target_domain       [c,d] to map u) into
     * @return  double                          u translated into [c,d] by a
     *                                          linear map
     */
    double affine_coordinate_transform(double u, 
            Interval source_domain, 
            Interval target_domain);
    /**
     * The same as above, but performed on both components of uv. Could be
     * collapse into remap_samples(), but left around for the convenience of not
     * dump a Point2 into a NumVector. TODO kill this function and overload
     * remap_samples() instead.
     * @param   Point2      u                  coordinate in source_domain to  
     *                                          be mapped into target_domain
     * @param   Rectangle   source_domain       [a,b]x[c,d] containing uv
     * @param   Rectangle   target_domain       [e,f]x[g,h] to map uv into
     * @return  Point2                          u translated into [e,f]x[g,h] by a
     *                                          linear map
     */
    Point2 affine_coordinate_transform(Point2 uv, 
            Rectangle source_domain, 
            Rectangle target_domain);

    /**
     * Generate a set of num_samples^2 2D tensor-product equispaced samples.
     * @param    int    num_samples             number of sample points along
     *                                          each dimension (endpoints
     *                                          included).
     * @param    Rectangle domain               domain boundaries within which
     *                                          we'd like to sample
     * @return   NumMatrix                      2 x num_samples^2 list of sample
     *                                          points equispaced on [0,1]^2
     */
    //NumMatrix equispaced(int num_samples, Rectangle domain);

    /**
     * Applies affine_coordinate_transform element-wise to samples.
     */
    void remap_samples(Interval source_domain, Interval target_domain, 
            NumVector& samples);
    
    /**
     * Generate a set of num_samples samples on volume domain.
     * @template NumVector (*sampling_func)(int, Interval)
     *                                          function pointer to a 1d sampling 
     *                                          pattern defined on Interval.
     *                                          such as equispaced, chebyshev,
     *                                          etc. 
     * @param    int    num_samples             number of sample points along
     *                                          each dimension (endpoints
     *                                          included).
     * @param    Cube    domain              domain boundaries within which
     *                                          we'd like to sample
     * @return   NumVector                      num_samples list of sample
     *                                          points equispaced on domain 
     */
    template <NumVector (*sampling_func)(int, Interval)>
    NumMatrix sample_3d(int num_samples, Cube domain){
        int parameter_dim = 3;
        NumMatrix samples(parameter_dim, num_samples*num_samples*num_samples);

        NumVector sample_points_x = sampling_func(num_samples, get<0>(domain));
        NumVector sample_points_y = sampling_func(num_samples, get<1>(domain));
        NumVector sample_points_z = sampling_func(num_samples, get<2>(domain));

        for(int k = 0; k < num_samples; k++){
            for(int j = 0; j < num_samples; j++){
                for(int i = 0; i < num_samples; i++){
                    int index = k*num_samples*num_samples + j*num_samples + i;
                    Point3 sample_point(sample_points_x(i), sample_points_y(j), sample_points_z(k));

                    for(int d =0; d < parameter_dim; d++)
                        samples(d,index) = sample_point(d);

                }
            }
        }
        return samples;
    }
    
    /**
     * Generate a set of num_samples samples on domain.
     * @template NumVector (*sampling_func)(int, Interval)
     *                                          function pointer to a 1d sampling 
     *                                          pattern defined on Interval.
     *                                          such as equispaced, chebyshev,
     *                                          etc. 
     * @param    int    num_samples             number of sample points along
     *                                          each dimension (endpoints
     *                                          included).
     * @param    Interval   domain              domain boundaries within which
     *                                          we'd like to sample
     * @return   NumVector                      num_samples list of sample
     *                                          points equispaced on domain 
     */
    template <NumVector (*sampling_func)(int, Interval)>
    NumMatrix sample_2d(int num_samples, Rectangle domain){
        NumMatrix samples(parameter_dim, num_samples*num_samples);

        NumVector sample_points_x = sampling_func(num_samples, domain.first);
        NumVector sample_points_y = sampling_func(num_samples, domain.second);

        for(int j = 0; j < num_samples; j++){
            for(int i = 0; i < num_samples; i++){
                int index = j*num_samples + i;
                Point2 sample_point(sample_points_x(i), sample_points_y(j));

                for(int d =0; d < parameter_dim; d++)
                    samples(d,index) = sample_point(d);

            }
        }
        return samples;
    }

    template <NumVector (*sampling_func)(int, Interval)>
    NumVector sample_1d(int num_samples, Rectangle domain){
        return sampling_func(num_samples, domain.first);
    }



    /**
     * Generate a set of num_samples equispaced points on domain [a,b].
     * endpoints are included.
     * @param    int    num_samples             number of sample points along
     *                                          each dimension (endpoints
     *                                          included).
     * @param    Interval   domain              domain boundaries within which
     *                                          we'd like to sample
     * @return   NumVector                      num_samples list of sample
     *                                          points equispaced on domain 
     */
    NumVector equispaced(int num_samples, Interval domain);
    
    /**
     * Generate a set of num_samples Chebyshev points on domain [a,b]. endpoints
     * are included.
     * @param    int    num_samples             number of sample points along
     *                                          each dimension (endpoints
     *                                          included).
     * @param    Interval   domain              domain boundaries within which
     *                                          we'd like to sample
     * @return   NumVector                      num_samples list of sample
     *                                          points equispaced on domain 
     */
    NumVector chebyshev1(int num_samples, Interval domain);
    NumVector chebyshev2(int num_samples, Interval domain);
    
};
#endif
