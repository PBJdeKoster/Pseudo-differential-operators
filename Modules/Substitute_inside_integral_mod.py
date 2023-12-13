# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 10:34:37 2020

@author: Pascal de Koster
"""

from sympy import symbols, Function, diff, I, simplify, \
Matrix, expand, integrate, pprint, Integral, preorder_traversal

''' Substitute_inside_integral
Input arguments
fun:        (Matrix) function with a symbolic expression, possibly containing 
            an integral
term_old:   The term to be replaced, but which .subs() does not cover. 
term_new:   The new term

Ouput arguments
fun_new:    (Matrix) function with term_old replaced with term_new, BUT ONLY
            INSIDE THE INTEGRALS!

This function only replaces terms within integrals, which are otherwise skipped
by .subs(). Refrain from using this function without first using .subs(). This
is not a standalone replacement for .subs(), but rather an addition.
'''

def Substitute_inside_integral(fun, term_old, term_new ):
    fun_new = fun.expand()
    #Check if fun has multiple entries, such as in a matrix or list. Make sure
    #to return the same type for the output as the input
    if hasattr(fun, '__len__'):
        N = len(fun)
    else:
        fun_new = Matrix([fun_new])
        N = 1
    
    #Creat a list of terms present in the expression, e.g.
    #f = a_x*int(b*a_x, x)  + 1
    #=> [add:       a_x*int(b*a_x) dx+1,
    #       mul:       a_x*int(b*a_x,x),
    #           diff:      diff(a,x)    ,
    #               function:   a,
    #               tuple:      (x,1),
    #                  tuple:     x,
    #                  number:    1,
    #           integral:  integral(b*a_x,x),
    #               mul: b*a_x,
    #                   function: b,
    #                   diff:     (a,x),
    #                       function:  a,
    #                       tuple:     (x,1),
    #                           symbol:   x,
    #                           number:   x,
    #               symbol: x,
    #               number: 1,
    #       number:    1]
    #The list consists only of the elements on the right.
    #The types in front of each element are not included; I added them for 
    #clarification. The indentation is also my own doing, to show the hierachy;
    #this is also not included in the list of elements.
    
    for i_N in range(0,N):
        fun_tmp = fun_new[i_N]
        List_of_elements = [el for el in preorder_traversal(fun_tmp)]
        for i_el in range(0,len(List_of_elements)):
            el_i = List_of_elements[i_el]
            if isinstance(el_i, Integral):
                el_i_args = el_i.args
                int_f_old = el_i_args[0]
                int_dx = el_i_args[1]
                int_f_new = int_f_old.subs(term_old,term_new)
                if not (int_f_new == int_f_old):
                    fun_new[i_N] = fun_new[i_N].subs(\
                           Integral(int_f_old,int_dx),\
                           Integral(int_f_new,int_dx)\
                           )    
        fun_new[i_N] = fun_new[i_N].expand()#.doit(manual=True)
    
    if not hasattr(fun, '__len__'):
        fun_new = fun_new[0]
    
    return fun_new
