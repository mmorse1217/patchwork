import sympy as sp
import re
from sympy.utilities.codegen import codegen
from sympy.abc import x
from sympy.functions.special.bsplines import bspline_basis 
from sympy.polys.polyfuncs import horner
from sympy.functions.elementary.piecewise import Piecewise
from sympy import plot


def plot_bsp(n,k=0):
    m = 10
    knots = range(-n, m+n+1)
    bsp = bspline_basis(n, knots, 0, x)
    ddx_bsp = sp.diff(bsp, x, k)
    bounds = [(i,i+1) for i in xrange(-n,2)]
    p = plot(0, (x,-1,1))
    for i in range(len(ddx_bsp.args)):
        print i
        print ddx_bsp.args[i]
        lower_bound, upper_bound = bounds[i]
        p.append(sp.plot(ddx_bsp.args[i][0], (x, lower_bound, upper_bound))[0])
    p.show()
#plot_bsp(5,1)
    
def generate_bspline_code(n, k=0):
    m = 10
    knots = range(-n, m+n+1)
    bsp = bspline_basis(n, knots, 0, x)
    ddx_bsp = sp.diff(bsp, x, k)
    new_expr_cond_pairs = []
    
    # express each piecewise part of B-spline in horner form
    for i in range(len(ddx_bsp.args)):
        expr_cond_pair = ( horner(ddx_bsp.args[i][0]), ddx_bsp.args[i][1])
        new_expr_cond_pairs.append(expr_cond_pair)

    bsp_horner = Piecewise(*new_expr_cond_pairs)

    [(c_name, c_code), (h_name, c_header)] = codegen(
            ("bspline_deg"+str(n)+'_deriv'+str(k), bsp_horner),
            "C", 
            "bspline_codegen", 
            header=False, 
            empty=False)
    return c_name, c_code, h_name, c_header

def post_process(c_code, c_header):
    pattern = '\#[A-z \"\_\.\<\>]+'
    #kill off header includes
    c_code = re.sub(pattern, '', c_code)
    #kill off header guard 
    c_header = re.sub(pattern, '', c_header)
    return c_code, c_header

def dump_code():
    cpp_name = 'bspline_codegen.cpp'
    hpp_name = 'bspline_codegen.hpp'
    cpp_code = '#include "bspline_codegen.hpp"\n'
    hpp_code = '#ifndef __BSPLINE_CODEGEN_HPP__\n #define __BSPLINE_CODEGEN_HPP__\n#include <math.h>\n'
    n = 5
    for i in xrange(3, 6):
        for k in xrange(3):
            _,ith_c_code, __, ith_c_header = generate_bspline_code(i,k)
            ith_c_code, ith_c_header = post_process(ith_c_code, ith_c_header)
            cpp_code += ith_c_code
            hpp_code += ith_c_header

    hpp_code += '#endif\n'

    f = open(cpp_name, 'w')
    f.write(cpp_code)
    f.close()

    f = open(hpp_name, 'w')
    f.write(hpp_code)
    f.close()

    print cpp_code
    print hpp_code
