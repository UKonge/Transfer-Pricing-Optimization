# -*- coding: utf-8 -*-
"""
Created on Sun Apr 25 09:51:40 2021

@author: Utkarsh

Here, the MILP model for transfer pricing optimization is developed. 
The goal is to decide the production quantities and transfer prices among various subsidies to maximize 
the profit after tax on income in various countries.

"""

from pyomo.environ import *

comp_centers = 4
fac_centers = 2
dist_centers = 3


m = ConcreteModel()

# Variables
m.btp_comp = Var(RangeSet(comp_centers),domain=NonNegativeReals)
m.btp_fac = Var(RangeSet(fac_centers),domain=NonNegativeReals)
m.btp_dist = Var(RangeSet(dist_centers),domain=NonNegativeReals)
m.btl_comp = Var(RangeSet(comp_centers),domain=NonNegativeReals)
m.btl_fac = Var(RangeSet(fac_centers),domain=NonNegativeReals)
m.btl_dist = Var(RangeSet(dist_centers),domain=NonNegativeReals)
m.g = Var(RangeSet(comp_centers),domain=NonNegativeReals)
m.h = Var(RangeSet(fac_centers),domain=NonNegativeReals)
m.x = Var(RangeSet(comp_centers),RangeSet(fac_centers),domain=NonNegativeReals)
m.y = Var(RangeSet(fac_centers),RangeSet(dist_centers),domain=NonNegativeReals)
m.tp_fc = Var(RangeSet(comp_centers),RangeSet(fac_centers),domain=NonNegativeReals)
m.tp_df = Var(RangeSet(fac_centers),RangeSet(dist_centers),domain=NonNegativeReals)


