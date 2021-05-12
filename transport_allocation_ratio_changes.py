# -*- coding: utf-8 -*-
"""
Created on Fri May  7 21:49:27 2021

@author: 91824
"""
from pyomo.environ import *
import pandas as pd

comp_centers = 4
fac_centers = 2
dist_centers = 3

D = [1200, 1000, 900]
QF = [2000, 2000]
QC = [2000, 2000, 2000, 2000]
a = [[30, 100],[70, 90],[70, 90],[110, 40]]
b = [[40, 150, 160],[210, 200, 180]]
s = [800, 830, 880]
c = [80, 70, 80, 50]
pC = [50, 30, 30, 20]
e = [140, 120]
n = [0.65, 0.8, 0.8, 0.8]
m1 = [0.65, 0.8]
r = [0.65, 0.7, 0.75]
lbC = [[160, 230],[170, 190],[180, 200],[180, 110]]
ubC = [[340, 530],[350, 490],[360, 500],[360, 410]]
lbF = [[520, 630, 640],[630, 620, 600]]
ubF = [[800, 830, 880],[800, 830, 880]]
f_fc_list = [1.0]
f_df_list = [1.0]
profits = {}

for f_fc in f_fc_list:
    for f_df in f_df_list:

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
        m.x = Var(RangeSet(comp_centers),RangeSet(fac_centers),domain=NonNegativeIntegers)
        m.y = Var(RangeSet(fac_centers),RangeSet(dist_centers),domain=NonNegativeIntegers)
        m.TP_fc = Var(RangeSet(comp_centers),RangeSet(fac_centers),domain=NonNegativeReals)
        m.TP_df = Var(RangeSet(fac_centers),RangeSet(dist_centers),domain=NonNegativeReals)
        
        # Constraints
        # Capacity Constraints
        m.cap_cons_comp = ConstraintList()
        for i in range(comp_centers):
            m.cap_cons_comp.add(expr=m.g[i+1]<=QC[i])
        
        m.cap_cons_fac = ConstraintList()
        for j in range(fac_centers):
            m.cap_cons_fac.add(expr=m.h[j+1]<=QF[j])
            
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
            m.int_final_bal.add(expr=m.x[1,j+1]+m.x[2,j+1] == sum(m.y[j+1,k+1] for k in range(dist_centers))) # Component 1
        for j in range(fac_centers):
            m.int_final_bal.add(expr=m.x[3,j+1]+m.x[4,j+1] == sum(m.y[j+1,k+1] for k in range(dist_centers))) # Component 2
            
        # Demand satisfaction
        m.dem_cons = ConstraintList()
        for k in range(dist_centers):
            m.dem_cons.add(expr=sum(m.y[j+1,k+1] for j in range(fac_centers))==D[k])
            
        # Transfer prices bounds
        m.tp_bnds_comp = ConstraintList()
        for i in range(comp_centers):
            for j in range(fac_centers):
                m.tp_bnds_comp.add(expr=lbC[i][j]*m.x[i+1,j+1] <= m.TP_fc[i+1,j+1])
            
        for i in range(comp_centers):
            for j in range(fac_centers):
                m.tp_bnds_comp.add(expr=m.TP_fc[i+1,j+1] <= ubC[i][j]*m.x[i+1,j+1])
        
        m.tp_bnds_fac = ConstraintList()
        for j in range(fac_centers):
            for k in range(dist_centers):
                m.tp_bnds_fac.add(expr=lbF[j][k]*m.y[j+1,k+1] <= m.TP_df[j+1,k+1])
        
        for j in range(fac_centers):
            for k in range(dist_centers):
                m.tp_bnds_fac.add(expr=m.TP_df[j+1,k+1] <= ubF[j][k]*m.y[j+1,k+1])
        
        # Profit before tax computation constraints
        m.comp_btp_cons = ConstraintList() # component production centers profit before tax
        for i in range(comp_centers):
            m.comp_btp_cons.add(expr=m.btp_comp[i+1]-m.btl_comp[i+1] == sum((m.TP_fc[i+1,j+1]-f_fc*a[i][j]*m.x[i+1,j+1])for j in range(fac_centers))-m.g[i+1]*(pC[i]+c[i]))
                                
        m.fac_btp_cons = ConstraintList() # Final product production facilities
        for j in range(fac_centers):
            m.fac_btp_cons.add(expr=m.btp_fac[j+1]-m.btl_fac[j+1] == sum((m.TP_df[j+1,k+1]-f_df*b[j][k]*m.y[j+1,k+1]) for k in range(dist_centers)) - e[j]*m.h[j+1] - sum((m.TP_fc[i+1,j+1]+(1-f_fc)*a[i][j]*m.x[i+1,j+1]) for i in range(comp_centers)))
            
        m.dist_btp_cons = ConstraintList()
        for k in range(dist_centers):
            m.dist_btp_cons.add(expr= m.btp_dist[k+1]-m.btl_dist[k+1] == sum((s[k]*m.y[j+1,k+1]-m.TP_df[j+1,k+1]-(1-f_df)*b[j][k]*m.y[j+1,k+1])for j in range(fac_centers)))
            
        # Objective
        #m.obj = Objective(expr = sum((1-n[i])*m.btp_comp[i+1]-m.btl_comp[i+1] for i in range(comp_centers))+sum((1-m1[j])*m.btp_fac[j+1]-m.btl_fac[j+1] for j in range(fac_centers))+sum((1-r[k])*m.btp_dist[k+1]-m.btl_dist[k+1] for k in range(dist_centers)),sense=maximize)
        
        m.obj = Objective(expr = sum((n[i])*m.btp_comp[i+1]-m.btl_comp[i+1] for i in range(comp_centers))+sum((m1[j])*m.btp_fac[j+1]-m.btl_fac[j+1] for j in range(fac_centers))+sum((r[k])*m.btp_dist[k+1]-m.btl_dist[k+1] for k in range(dist_centers)),sense=maximize)
        
        # Solving the model
        opt = SolverFactory('cplex')
        results = opt.solve(m,tee=True)
        print("Solver status:", results.solver.status)
        print("Solver termination condition:", results.solver.termination_condition)
        
        
        cols = ['Variable','Variable Name','Number','TP_variable','TP_fixed_mean','TP_fixed_min','TP_fixed_max']
        sol = pd.DataFrame(data=None,columns=cols,index=range(35))
        
        ind = 0
        
        for i in range(comp_centers):
            sol.at[ind,'Variable'] = 'Comp Prod'
            sol.at[ind,'Variable Name'] = 'g'
            sol.at[ind,'Number'] = i+1
            sol.at[ind,'TP_variable'] = m.g[i+1].value
            ind += 1
            
        for i in range(fac_centers):
            sol.at[ind,'Variable'] = 'Final Prod'
            sol.at[ind,'Variable Name'] = 'h'
            sol.at[ind,'Number'] = i+1
            sol.at[ind,'TP_variable'] = m.h[i+1].value
            ind += 1
        
        for i in range(comp_centers):
            for j in range(fac_centers):
                sol.at[ind,'Variable'] = 'Comp ship'
                sol.at[ind,'Number'] = (i+1,j+1)
                sol.at[ind,'Variable Name'] = 'x'
                sol.at[ind,'TP_variable'] = m.x[i+1,j+1].value
                ind += 1
        
        for i in range(fac_centers):
            for j in range(dist_centers):
                sol.at[ind,'Variable'] = 'Prod ship'
                sol.at[ind,'Number'] = (i+1,j+1)
                sol.at[ind,'Variable Name'] = 'y'
                sol.at[ind,'TP_variable'] = m.y[i+1,j+1].value
                ind += 1
        
        for i in range(comp_centers):
            for j in range(fac_centers):
                sol.at[ind,'Variable'] = 'TP c,f'
                sol.at[ind,'Number'] = (i+1,j+1)
                sol.at[ind,'Variable Name'] = 'tp'
                if m.TP_fc[i+1,j+1].value == 0 and m.x[i+1,j+1].value == 0:
                    sol.at[ind,'TP_variable'] = 0
                else:
                    sol.at[ind,'TP_variable'] = m.TP_fc[i+1,j+1].value/m.x[i+1,j+1].value
                ind += 1
            
        for i in range(fac_centers):
            for j in range(dist_centers):
                sol.at[ind,'Variable'] = 'TP f,d'
                sol.at[ind,'Number'] = (i+1,j+1)
                sol.at[ind,'Variable Name'] = 'tp'
                if m.TP_df[i+1,j+1].value == 0 and m.y[i+1,j+1].value == 0:
                    sol.at[ind,'TP_variable'] = 0
                else:
                    sol.at[ind,'TP_variable'] = m.TP_df[i+1,j+1].value/m.y[i+1,j+1].value
                ind += 1
                
        sol.at[ind,'Variable'] = 'Profit'
        sol.at[ind,'TP_variable'] = m.obj() 
        
        
        print("Profits after tax")
        print("  Component")
        for i in range(comp_centers):
            print("    Comp manuf:",i+1,"Profit:",n[i]*m.btp_comp[i+1].value,"Loss:",m.btl_comp[i+1].value)
        print("  Facilities")
        for i in range(fac_centers):
            print("    Fac manuf:",i+1,"Profit:",m1[i]*m.btp_fac[i+1].value,"Loss:",m.btl_fac[i+1].value)
        print("  Distribution centers")
        for i in range(dist_centers):
            print("    Dist center:",i+1,"Profit:",r[i]*m.btp_dist[i+1].value,"Loss:",m.btl_dist[i+1].value)
        
        profits[f_fc,f_df] = m.obj() 

'''
##########################################################
##########################################################
##########################################################

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
m.x = Var(RangeSet(comp_centers),RangeSet(fac_centers),domain=NonNegativeIntegers)
m.y = Var(RangeSet(fac_centers),RangeSet(dist_centers),domain=NonNegativeIntegers)
m.TP_fc = Var(RangeSet(comp_centers),RangeSet(fac_centers),domain=NonNegativeReals)
m.TP_df = Var(RangeSet(fac_centers),RangeSet(dist_centers),domain=NonNegativeReals)

# Constraints
# Capacity Constraints
m.cap_cons_comp = ConstraintList()
for i in range(comp_centers):
    m.cap_cons_comp.add(expr=m.g[i+1]<=QC[i])

m.cap_cons_fac = ConstraintList()
for j in range(fac_centers):
    m.cap_cons_fac.add(expr=m.h[j+1]<=QF[j])
    
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
    m.int_final_bal.add(expr=m.x[1,j+1]+m.x[2,j+1] == sum(m.y[j+1,k+1] for k in range(dist_centers))) # Component 1
for j in range(fac_centers):
    m.int_final_bal.add(expr=m.x[3,j+1]+m.x[4,j+1] == sum(m.y[j+1,k+1] for k in range(dist_centers))) # Component 2
    
# Demand satisfaction
m.dem_cons = ConstraintList()
for k in range(dist_centers):
    m.dem_cons.add(expr=sum(m.y[j+1,k+1] for j in range(fac_centers))==D[k])
    
# Transfer prices bounds
m.tp_bnds_comp = ConstraintList()
for i in range(comp_centers):
    for j in range(fac_centers):
        m.tp_bnds_comp.add(expr=(lbC[i][j]+ubC[i][j])/2*m.x[i+1,j+1] <= m.TP_fc[i+1,j+1])
    
#for i in range(comp_centers):
#    for j in range(fac_centers):
#        m.tp_bnds_comp.add(expr=m.TP_fc[i+1,j+1] <= ubC[i][j]*m.x[i+1,j+1])

m.tp_bnds_fac = ConstraintList()
for j in range(fac_centers):
    for k in range(dist_centers):
        m.tp_bnds_fac.add(expr=(lbF[j][k]+ubF[j][k])/2*m.y[j+1,k+1] == m.TP_df[j+1,k+1])

#for j in range(fac_centers):
#    for k in range(dist_centers):
#        m.tp_bnds_fac.add(expr=m.TP_df[j+1,k+1] <= ubF[j][k]*m.y[j+1,k+1])

# Profit before tax computation constraints
m.comp_btp_cons = ConstraintList() # component production centers profit before tax
for i in range(comp_centers):
    m.comp_btp_cons.add(expr=m.btp_comp[i+1]-m.btl_comp[i+1] == sum((m.TP_fc[i+1,j+1]-f_fc*a[i][j]*m.x[i+1,j+1])for j in range(fac_centers))-m.g[i+1]*(pC[i]+c[i]))
                        
m.fac_btp_cons = ConstraintList() # Final product production facilities
for j in range(fac_centers):
    m.fac_btp_cons.add(expr=m.btp_fac[j+1]-m.btl_fac[j+1] == sum((m.TP_df[j+1,k+1]-f_df*b[j][k]*m.y[j+1,k+1]) for k in range(dist_centers)) - e[j]*m.h[j+1] - sum((m.TP_fc[i+1,j+1]+(1-f_fc)*a[i][j]*m.x[i+1,j+1]) for i in range(comp_centers)))
    
m.dist_btp_cons = ConstraintList()
for k in range(dist_centers):
    m.dist_btp_cons.add(expr= m.btp_dist[k+1]-m.btl_dist[k+1] == sum((s[k]*m.y[j+1,k+1]-m.TP_df[j+1,k+1]-(1-f_df)*b[j][k]*m.y[j+1,k+1])for j in range(fac_centers)))
    
# Objective
#m.obj = Objective(expr = sum((1-n[i])*m.btp_comp[i+1]-m.btl_comp[i+1] for i in range(comp_centers))+sum((1-m1[j])*m.btp_fac[j+1]-m.btl_fac[j+1] for j in range(fac_centers))+sum((1-r[k])*m.btp_dist[k+1]-m.btl_dist[k+1] for k in range(dist_centers)),sense=maximize)

m.obj = Objective(expr = sum((n[i])*m.btp_comp[i+1]-m.btl_comp[i+1] for i in range(comp_centers))+sum((m1[j])*m.btp_fac[j+1]-m.btl_fac[j+1] for j in range(fac_centers))+sum((r[k])*m.btp_dist[k+1]-m.btl_dist[k+1] for k in range(dist_centers)),sense=maximize)

# Solving the model
opt = SolverFactory('cplex')
results = opt.solve(m,tee=True)
print("Solver status:", results.solver.status)
print("Solver termination condition:", results.solver.termination_condition)


ind = 0

for i in range(comp_centers):
    #sol.at[ind,'Variable'] = 'Comp Prod'
    #sol.at[ind,'Variable Name'] = 'g'
    #sol.at[ind,'Number'] = i+1
    sol.at[ind,'TP_fixed_mean'] = m.g[i+1].value
    ind += 1
    
for i in range(fac_centers):
    #sol.at[ind,'Variable'] = 'Final Prod'
    #sol.at[ind,'Variable Name'] = 'h'
    #sol.at[ind,'Number'] = i+1
    sol.at[ind,'TP_fixed_mean'] = m.h[i+1].value
    ind += 1

for i in range(comp_centers):
    for j in range(fac_centers):
        #sol.at[ind,'Variable'] = 'Comp ship'
        #sol.at[ind,'Number'] = (i+1,j+1)
        #sol.at[ind,'Variable Name'] = 'x'
        sol.at[ind,'TP_fixed_mean'] = m.x[i+1,j+1].value
        ind += 1

for i in range(fac_centers):
    for j in range(dist_centers):
        #sol.at[ind,'Variable'] = 'Prod ship'
        #sol.at[ind,'Number'] = (i+1,j+1)
        #sol.at[ind,'Variable Name'] = 'y'
        sol.at[ind,'TP_fixed_mean'] = m.y[i+1,j+1].value
        ind += 1

for i in range(comp_centers):
    for j in range(fac_centers):
        #sol.at[ind,'Variable'] = 'TP c,f'
        #sol.at[ind,'Number'] = (i+1,j+1)
        #sol.at[ind,'Variable Name'] = 'tp'
        if m.TP_fc[i+1,j+1].value == 0 and m.x[i+1,j+1].value == 0:
            sol.at[ind,'TP_fixed_mean'] = 0
        else:
            sol.at[ind,'TP_fixed_mean'] = m.TP_fc[i+1,j+1].value/m.x[i+1,j+1].value
        ind += 1
    
for i in range(fac_centers):
    for j in range(dist_centers):
        #sol.at[ind,'Variable'] = 'TP f,d'
        #sol.at[ind,'Number'] = (i+1,j+1)
        #sol.at[ind,'Variable Name'] = 'tp'
        if m.TP_df[i+1,j+1].value == 0 and m.y[i+1,j+1].value == 0:
            sol.at[ind,'TP_fixed_mean'] = 0
        else:
            sol.at[ind,'TP_fixed_mean'] = m.TP_df[i+1,j+1].value/m.y[i+1,j+1].value
        ind += 1
        
#sol.at[ind,'Variable'] = 'Profit'
sol.at[ind,'TP_fixed_mean'] = m.obj() 

########################################################
########################################################
########################################################
########################################################

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
m.x = Var(RangeSet(comp_centers),RangeSet(fac_centers),domain=NonNegativeIntegers)
m.y = Var(RangeSet(fac_centers),RangeSet(dist_centers),domain=NonNegativeIntegers)
m.TP_fc = Var(RangeSet(comp_centers),RangeSet(fac_centers),domain=NonNegativeReals)
m.TP_df = Var(RangeSet(fac_centers),RangeSet(dist_centers),domain=NonNegativeReals)

# Constraints
# Capacity Constraints
m.cap_cons_comp = ConstraintList()
for i in range(comp_centers):
    m.cap_cons_comp.add(expr=m.g[i+1]<=QC[i])

m.cap_cons_fac = ConstraintList()
for j in range(fac_centers):
    m.cap_cons_fac.add(expr=m.h[j+1]<=QF[j])
    
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
    m.int_final_bal.add(expr=m.x[1,j+1]+m.x[2,j+1] == sum(m.y[j+1,k+1] for k in range(dist_centers))) # Component 1
for j in range(fac_centers):
    m.int_final_bal.add(expr=m.x[3,j+1]+m.x[4,j+1] == sum(m.y[j+1,k+1] for k in range(dist_centers))) # Component 2
    
# Demand satisfaction
m.dem_cons = ConstraintList()
for k in range(dist_centers):
    m.dem_cons.add(expr=sum(m.y[j+1,k+1] for j in range(fac_centers))==D[k])
    
# Transfer prices bounds
m.tp_bnds_comp = ConstraintList()
for i in range(comp_centers):
    for j in range(fac_centers):
        m.tp_bnds_comp.add(expr=(lbC[i][j])*m.x[i+1,j+1] <= m.TP_fc[i+1,j+1])
    
#for i in range(comp_centers):
#    for j in range(fac_centers):
#        m.tp_bnds_comp.add(expr=m.TP_fc[i+1,j+1] <= ubC[i][j]*m.x[i+1,j+1])

m.tp_bnds_fac = ConstraintList()
for j in range(fac_centers):
    for k in range(dist_centers):
        m.tp_bnds_fac.add(expr=(lbF[j][k])*m.y[j+1,k+1] == m.TP_df[j+1,k+1])

#for j in range(fac_centers):
#    for k in range(dist_centers):
#        m.tp_bnds_fac.add(expr=m.TP_df[j+1,k+1] <= ubF[j][k]*m.y[j+1,k+1])

# Profit before tax computation constraints
m.comp_btp_cons = ConstraintList() # component production centers profit before tax
for i in range(comp_centers):
    m.comp_btp_cons.add(expr=m.btp_comp[i+1]-m.btl_comp[i+1] == sum((m.TP_fc[i+1,j+1]-f_fc*a[i][j]*m.x[i+1,j+1])for j in range(fac_centers))-m.g[i+1]*(pC[i]+c[i]))
                        
m.fac_btp_cons = ConstraintList() # Final product production facilities
for j in range(fac_centers):
    m.fac_btp_cons.add(expr=m.btp_fac[j+1]-m.btl_fac[j+1] == sum((m.TP_df[j+1,k+1]-f_df*b[j][k]*m.y[j+1,k+1]) for k in range(dist_centers)) - e[j]*m.h[j+1] - sum((m.TP_fc[i+1,j+1]+(1-f_fc)*a[i][j]*m.x[i+1,j+1]) for i in range(comp_centers)))
    
m.dist_btp_cons = ConstraintList()
for k in range(dist_centers):
    m.dist_btp_cons.add(expr= m.btp_dist[k+1]-m.btl_dist[k+1] == sum((s[k]*m.y[j+1,k+1]-m.TP_df[j+1,k+1]-(1-f_df)*b[j][k]*m.y[j+1,k+1])for j in range(fac_centers)))
    
# Objective
#m.obj = Objective(expr = sum((1-n[i])*m.btp_comp[i+1]-m.btl_comp[i+1] for i in range(comp_centers))+sum((1-m1[j])*m.btp_fac[j+1]-m.btl_fac[j+1] for j in range(fac_centers))+sum((1-r[k])*m.btp_dist[k+1]-m.btl_dist[k+1] for k in range(dist_centers)),sense=maximize)

m.obj = Objective(expr = sum((n[i])*m.btp_comp[i+1]-m.btl_comp[i+1] for i in range(comp_centers))+sum((m1[j])*m.btp_fac[j+1]-m.btl_fac[j+1] for j in range(fac_centers))+sum((r[k])*m.btp_dist[k+1]-m.btl_dist[k+1] for k in range(dist_centers)),sense=maximize)

# Solving the model
opt = SolverFactory('cplex')
results = opt.solve(m,tee=True)
print("Solver status:", results.solver.status)
print("Solver termination condition:", results.solver.termination_condition)


ind = 0

for i in range(comp_centers):
    #sol.at[ind,'Variable'] = 'Comp Prod'
    #sol.at[ind,'Variable Name'] = 'g'
    #sol.at[ind,'Number'] = i+1
    sol.at[ind,'TP_fixed_min'] = m.g[i+1].value
    ind += 1
    
for i in range(fac_centers):
    #sol.at[ind,'Variable'] = 'Final Prod'
    #sol.at[ind,'Variable Name'] = 'h'
    #sol.at[ind,'Number'] = i+1
    sol.at[ind,'TP_fixed_min'] = m.h[i+1].value
    ind += 1

for i in range(comp_centers):
    for j in range(fac_centers):
        #sol.at[ind,'Variable'] = 'Comp ship'
        #sol.at[ind,'Number'] = (i+1,j+1)
        #sol.at[ind,'Variable Name'] = 'x'
        sol.at[ind,'TP_fixed_min'] = m.x[i+1,j+1].value
        ind += 1

for i in range(fac_centers):
    for j in range(dist_centers):
        #sol.at[ind,'Variable'] = 'Prod ship'
        #sol.at[ind,'Number'] = (i+1,j+1)
        #sol.at[ind,'Variable Name'] = 'y'
        sol.at[ind,'TP_fixed_min'] = m.y[i+1,j+1].value
        ind += 1

for i in range(comp_centers):
    for j in range(fac_centers):
        #sol.at[ind,'Variable'] = 'TP c,f'
        #sol.at[ind,'Number'] = (i+1,j+1)
        #sol.at[ind,'Variable Name'] = 'tp'
        if m.TP_fc[i+1,j+1].value == 0 and m.x[i+1,j+1].value == 0:
            sol.at[ind,'TP_fixed_min'] = 0
        else:
            sol.at[ind,'TP_fixed_min'] = m.TP_fc[i+1,j+1].value/m.x[i+1,j+1].value
        ind += 1
    
for i in range(fac_centers):
    for j in range(dist_centers):
        #sol.at[ind,'Variable'] = 'TP f,d'
        #sol.at[ind,'Number'] = (i+1,j+1)
        #sol.at[ind,'Variable Name'] = 'tp'
        if m.TP_df[i+1,j+1].value == 0 and m.y[i+1,j+1].value == 0:
            sol.at[ind,'TP_fixed_min'] = 0
        else:
            sol.at[ind,'TP_fixed_min'] = m.TP_df[i+1,j+1].value/m.y[i+1,j+1].value
        ind += 1
        
#sol.at[ind,'Variable'] = 'Profit'
sol.at[ind,'TP_fixed_min'] = m.obj() 


########################################################
########################################################
########################################################
########################################################

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
m.x = Var(RangeSet(comp_centers),RangeSet(fac_centers),domain=NonNegativeIntegers)
m.y = Var(RangeSet(fac_centers),RangeSet(dist_centers),domain=NonNegativeIntegers)
m.TP_fc = Var(RangeSet(comp_centers),RangeSet(fac_centers),domain=NonNegativeReals)
m.TP_df = Var(RangeSet(fac_centers),RangeSet(dist_centers),domain=NonNegativeReals)

# Constraints
# Capacity Constraints
m.cap_cons_comp = ConstraintList()
for i in range(comp_centers):
    m.cap_cons_comp.add(expr=m.g[i+1]<=QC[i])

m.cap_cons_fac = ConstraintList()
for j in range(fac_centers):
    m.cap_cons_fac.add(expr=m.h[j+1]<=QF[j])
    
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
    m.int_final_bal.add(expr=m.x[1,j+1]+m.x[2,j+1] == sum(m.y[j+1,k+1] for k in range(dist_centers))) # Component 1
for j in range(fac_centers):
    m.int_final_bal.add(expr=m.x[3,j+1]+m.x[4,j+1] == sum(m.y[j+1,k+1] for k in range(dist_centers))) # Component 2
    
# Demand satisfaction
m.dem_cons = ConstraintList()
for k in range(dist_centers):
    m.dem_cons.add(expr=sum(m.y[j+1,k+1] for j in range(fac_centers))==D[k])
    
# Transfer prices bounds
m.tp_bnds_comp = ConstraintList()
for i in range(comp_centers):
    for j in range(fac_centers):
        m.tp_bnds_comp.add(expr=(ubC[i][j])*m.x[i+1,j+1] <= m.TP_fc[i+1,j+1])
    
#for i in range(comp_centers):
#    for j in range(fac_centers):
#        m.tp_bnds_comp.add(expr=m.TP_fc[i+1,j+1] <= ubC[i][j]*m.x[i+1,j+1])

m.tp_bnds_fac = ConstraintList()
for j in range(fac_centers):
    for k in range(dist_centers):
        m.tp_bnds_fac.add(expr=(ubF[j][k])*m.y[j+1,k+1] == m.TP_df[j+1,k+1])

#for j in range(fac_centers):
#    for k in range(dist_centers):
#        m.tp_bnds_fac.add(expr=m.TP_df[j+1,k+1] <= ubF[j][k]*m.y[j+1,k+1])

# Profit before tax computation constraints
m.comp_btp_cons = ConstraintList() # component production centers profit before tax
for i in range(comp_centers):
    m.comp_btp_cons.add(expr=m.btp_comp[i+1]-m.btl_comp[i+1] == sum((m.TP_fc[i+1,j+1]-f_fc*a[i][j]*m.x[i+1,j+1])for j in range(fac_centers))-m.g[i+1]*(pC[i]+c[i]))
                        
m.fac_btp_cons = ConstraintList() # Final product production facilities
for j in range(fac_centers):
    m.fac_btp_cons.add(expr=m.btp_fac[j+1]-m.btl_fac[j+1] == sum((m.TP_df[j+1,k+1]-f_df*b[j][k]*m.y[j+1,k+1]) for k in range(dist_centers)) - e[j]*m.h[j+1] - sum((m.TP_fc[i+1,j+1]+(1-f_fc)*a[i][j]*m.x[i+1,j+1]) for i in range(comp_centers)))
    
m.dist_btp_cons = ConstraintList()
for k in range(dist_centers):
    m.dist_btp_cons.add(expr= m.btp_dist[k+1]-m.btl_dist[k+1] == sum((s[k]*m.y[j+1,k+1]-m.TP_df[j+1,k+1]-(1-f_df)*b[j][k]*m.y[j+1,k+1])for j in range(fac_centers)))
    
# Objective
#m.obj = Objective(expr = sum((1-n[i])*m.btp_comp[i+1]-m.btl_comp[i+1] for i in range(comp_centers))+sum((1-m1[j])*m.btp_fac[j+1]-m.btl_fac[j+1] for j in range(fac_centers))+sum((1-r[k])*m.btp_dist[k+1]-m.btl_dist[k+1] for k in range(dist_centers)),sense=maximize)

m.obj = Objective(expr = sum((n[i])*m.btp_comp[i+1]-m.btl_comp[i+1] for i in range(comp_centers))+sum((m1[j])*m.btp_fac[j+1]-m.btl_fac[j+1] for j in range(fac_centers))+sum((r[k])*m.btp_dist[k+1]-m.btl_dist[k+1] for k in range(dist_centers)),sense=maximize)

# Solving the model
opt = SolverFactory('cplex')
results = opt.solve(m,tee=True)
print("Solver status:", results.solver.status)
print("Solver termination condition:", results.solver.termination_condition)


ind = 0

for i in range(comp_centers):
    #sol.at[ind,'Variable'] = 'Comp Prod'
    #sol.at[ind,'Variable Name'] = 'g'
    #sol.at[ind,'Number'] = i+1
    sol.at[ind,'TP_fixed_max'] = m.g[i+1].value
    ind += 1
    
for i in range(fac_centers):
    #sol.at[ind,'Variable'] = 'Final Prod'
    #sol.at[ind,'Variable Name'] = 'h'
    #sol.at[ind,'Number'] = i+1
    sol.at[ind,'TP_fixed_max'] = m.h[i+1].value
    ind += 1

for i in range(comp_centers):
    for j in range(fac_centers):
        #sol.at[ind,'Variable'] = 'Comp ship'
        #sol.at[ind,'Number'] = (i+1,j+1)
        #sol.at[ind,'Variable Name'] = 'x'
        sol.at[ind,'TP_fixed_max'] = m.x[i+1,j+1].value
        ind += 1

for i in range(fac_centers):
    for j in range(dist_centers):
        #sol.at[ind,'Variable'] = 'Prod ship'
        #sol.at[ind,'Number'] = (i+1,j+1)
        #sol.at[ind,'Variable Name'] = 'y'
        sol.at[ind,'TP_fixed_max'] = m.y[i+1,j+1].value
        ind += 1

for i in range(comp_centers):
    for j in range(fac_centers):
        #sol.at[ind,'Variable'] = 'TP c,f'
        #sol.at[ind,'Number'] = (i+1,j+1)
        #sol.at[ind,'Variable Name'] = 'tp'
        if m.TP_fc[i+1,j+1].value == 0 and m.x[i+1,j+1].value == 0:
            sol.at[ind,'TP_fixed_max'] = 0
        else:
            sol.at[ind,'TP_fixed_max'] = m.TP_fc[i+1,j+1].value/m.x[i+1,j+1].value
        ind += 1
    
for i in range(fac_centers):
    for j in range(dist_centers):
        #sol.at[ind,'Variable'] = 'TP f,d'
        #sol.at[ind,'Number'] = (i+1,j+1)
        #sol.at[ind,'Variable Name'] = 'tp'
        if m.TP_df[i+1,j+1].value == 0 and m.y[i+1,j+1].value == 0:
            sol.at[ind,'TP_fixed_max'] = 0
        else:
            sol.at[ind,'TP_fixed_max'] = m.TP_df[i+1,j+1].value/m.y[i+1,j+1].value
        ind += 1
        
#sol.at[ind,'Variable'] = 'Profit'
sol.at[ind,'TP_fixed_max'] = m.obj() 


'''