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
m.TP_fc = Var(RangeSet(comp_centers),RangeSet(fac_centers),domain=NonNegativeReals)
m.TP_df = Var(RangeSet(fac_centers),RangeSet(dist_centers),domain=NonNegativeReals)

# Constraints
# Capacity Constraints
m.cap_cons_comp = ConstraintList()
for i in range(comp_centers):
    m.cap_cons_comp.add(expr=m.g[i+1]<=QC[i])

m.cap_cons_fac = ConstraintList()
for i in range(comp_centers):
    m.cap_cons_fac.add(expr=m.h[i+1]<=QF[i])
    
# Balance constraints
m.bal_comp = ConstraintList()
for i in range(comp_centers):
    m.bal_comp.add(expr=sum(m.x[i+1,j+1] for j in range(fac_centers))==m.g[i+1])
    
m.bal_fac = ConstraintList()
for j in range(fac_centers):
    m.bal_comp.add(expr=sum(m.y[j+1,k+1] for k in range(dist_centers))==m.h[j+1])

# Intermediate products and final products balance constraints
m.int_final_bal = ConstraintList()
for j in range(fac_centers):
    m.int_final_bal.add(expr=m.x[1,j+1]+m.x[2,j+1] == m.y[j+1,k+1] for k in range(dist_centers))) # Component 1
for j in range(fac_centers):
    m.int_final_bal.add(expr=m.x[3,j+1]+m.x[4,j+1] == m.y[j+1,k+1] for k in range(dist_centers))) # Component 2
    
# Demand satisfaction
m.dem_cons = ConstraintList()
for k in range(dist_centers):
    m.dem_cons.add(expr=sum(m.y[j+1,k] for j in range(fac_centers))==D[k])
    
# Transfer prices bounds
m.tp_bnds_comp = ConstraintList()
for i in range(comp_centers):
    for j in range(fac_centers):
        m.tp_bnds_comp.add(expr=(lbC[i][j]*m.x[i+1,j+1],m.TP_fc[i+1,j+1],ubC[i][j]*m.x[i+1,j+1]))
    
m.tp_bnds_fac = ConstraintList()
for j in range(fac_centers):
    for k in range(dist_centers):
        m.tp_bnds_fac.add(expr=(lbF[i][j]*m.y[j+1,k+1],m.TP_df[j+1,k+1],ubF[i][j]*m.y[j+1,k+1))

# Profit before tax computation constraints
m.comp_btp_cons = ConstraintList() # component production centers profit before tax
for i in range(comp_centers):
    m.comp_btp_cons.add(expr=m.btp_comp[i+1]-m.btl_comp[i+1] == sum((m.TP_fc[i+1,j+1]-f_fc*a[i][j]*m.x[i+1,j+1])for j in range(comp_centers))-m.g[i]*(pC[i]+c[i]))
                        
m.fac_btp_cons = ConstraintList() # Final product production facilities
for j in range(fac_centers):
    m.fac_btp_cons.add(expr=m.btp_fac[j+1]-m.btl_fac[j+1] == sum((m.TP_df[j+1,k+1]-f_df*b[j][k]*m.y[j+1,k+1]) for k in range(dist_centers)) - e[j]*m.h[j+1] - sum((m.TP_fc[i+1,j+1]+(1-f_fc)*a[i][j]*m.x[i+1,j+1]) for i in range(comp_centers)))
    
m.dist_btp_cons = ConstraintList()
for k in range(dist_centers):
    m.dist_btp_cons.add(expr= m.btp_dist[k+1]-m.btl_dist[k+1] == sum((s[k]*m.y[j+1,k+1]-m.TP_df[j+1,k+1]-(1-f_df)*b[j][k]*m.y[j+1,k+1])for j in range(fac_centers)))
    
















