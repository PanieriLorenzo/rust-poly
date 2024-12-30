"""Generate quadratic equation formulas using sympy"""

from sympy import symbols, I, sqrt, re, simplify, im, expand, cse, trigsimp, expand_trig
import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)
from sympy_to_fpcore import codegen, parse


ra, rb, rc, ia, ib, ic = symbols("ra rb rc ia ib ic", rational=True)
re_sqrt_delta, im_sqrt_delta = symbols("re_sqrt_delta im_sqrt_delta", rational=True)
a = ra + ia * I
b = rb + ib * I
c = rc + ic * I
delta = sqrt(b**2 - 4 * a * c)
re_sqrt_delta = re(delta)
im_sqrt_delta = im(delta)
print(re_sqrt_delta)
print(im_sqrt_delta)

sqrt_delta = re_sqrt_delta + im_sqrt_delta * I
quad_plus = (-b + sqrt_delta) / (2 * a)
quad_minus = (-b - sqrt_delta) / (2 * a)
re_quad_plus = simplify(re(quad_plus))
# im_quad_plus = simplify(im(quad_plus))
# re_quad_minus = simplify(re(quad_minus))
# im_quad_minus = simplify(im(quad_minus))

# optimizations = {""}


# print(codegen(parse(str(re_quad_plus)), "re_quad_plus"))
# print()
# print(codegen(parse(str(im_quad_plus)), "im_quad_plus"))
# print()
# print(codegen(parse(str(re_quad_minus)), "re_quad_minus"))
# print()
# print(codegen(parse(str(im_quad_minus)), "im_quad_minus"))
# print()
