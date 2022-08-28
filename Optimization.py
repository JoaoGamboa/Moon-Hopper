# -*- coding: utf-8 -*-
"""
Created on Fri May  6 01:58:36 2022

@author: Joao Gamboa
"""

from gekko import GEKKO
import numpy as np
import matplotlib.pyplot as plt
import Initial_Guess


twr=1.3 #Thrust-to-Weight ratio
Isp=243 #Isp of the engine 

##Nt_values for different ranges, for a max final horizontal speed of 1 m/s
#3.2 -> 59
#2.501 -> 33,49
#2 -> 32
#1.499 -> 28
#1 -> 61
x_final=3.2 #Change here the values of the range of the hopper

m_total=100 #This value is kept as 100, in order to have the final mass in percentage, the algorithm does not depend on the initial mass of the hopper

#Initial constants 
g0=9.81 #Earth's Gravity
gm=1.62 #Moon's Gravity

deg2rad=np.pi/180
gamma0=(89)*deg2rad

n=0
def optimize_trajectory(nt,t_final,x_final):
    """
    This function optimizes the trajectory
    
    Inputs:
        t_final: Initial guess for the time of the trajectory (Output of the Initial Guess)
        x_final: Distance of the trajectory (Output of the initial guess)
        nt: number of steps
        
    Outputs:
        Plot of all the state variables as a function of time       
        
    """    
    
    global n
    m= GEKKO() #Initialize Gekko
    
    m.time = np.linspace(0,t_final,nt) #Time array
    p = np.zeros(nt) # mark final time point
    p[-1] = 1.0 #Array with 0's everywhere except 1 in the last entry
    final = m.Param(value=p)
    
    # Initialize all the Variables  
    x=m.Var(value=0,lb=-1e3,ub=1e3)
    y=m.Var(value=0,lb=0,ub=1e3)
    v = m.Var(value=0,lb=0,ub=1e2)
    gamma = m.Var(value=gamma0,lb=-360*deg2rad,ub=360*deg2rad)
    alpha=m.Var(value=0,lb=-360*deg2rad,ub=360*deg2rad)
    alpha_dot=m.Var(value=0,lb=-1e2,ub=1e2)
    theta=m.Var(value=gamma0,lb=-360*deg2rad,ub=360*deg2rad)
    theta_dot=m.Var(value=0,lb=-360*deg2rad,ub=360*deg2rad)
    mass=m.Var(value=m_total,lb=0,ub=m_total)

    m_dot_max=twr*m_total*gm/(Isp*g0)  	    

    #Initialize the control variables
    m_dot = m.MV(value=m_dot_max,lb=0,ub=m_dot_max)
    m_dot.STATUS=1
    
    M_rw_J = m.MV(value=0,lb=-1,ub=1) #This is the value of M_rw/J_z
    M_rw_J.STATUS=1
    
    
    # Equations  
    m.Equation(x.dt()==v*m.cos(gamma))
    m.Equation(y.dt()==v*m.sin(gamma))
    m.Equation(mass*v.dt()==(Isp*m_dot*g0)*m.cos(alpha)-mass*gm*m.sin(gamma))
    m.Equation(mass*v*gamma.dt()==(Isp*m_dot*g0)*m.sin(alpha)-mass*gm*m.cos(gamma))    
    m.Equation(alpha+gamma==theta)
    m.Equation(theta.dt()==theta_dot.dt())
    m.Equation(theta_dot.dt()==M_rw_J)
    m.Equation(mass.dt()==-m_dot)
    
    #Constraints
    
    m.Minimize(final*(x-x_final)**2)
    m.Minimize(final*(alpha+gamma-90*deg2rad)**2)    
    m.Equation(m.abs(v*m.sin(gamma))*final<=3)
    # m.Equation(m.abs(v*m.cos(gamma))*final<=0)  #Final Horizontal Velocity less than 0 m/s
    m.Equation(y*final<=0)
    m.Obj(-mass*final) # Objective function

    #Options for the solver    
    m.options.MAX_ITER = 1e3     # adjust maximum iterations
    m.options.IMODE = 6 # optimal control mode
    m.options.NODES = 3 #Ipopt solver
    
    m.solve(disp=False,debug=False) # solve with IPOPT
        
    print("Max of M_Rw/J_z is {} for nt={}".format(max([abs(i) for i in M_rw_J.value]),nt))
        
    #Get the node of the t_bal in the time array 
    t_bal_node=np.where(np.array(m_dot.value)==min(m_dot.value))[0][0] 
    # t_bal_node=len(m.time)-1
    print(t_bal_node)
    print(m.time[t_bal_node])
    
    #Plot all the state variables
    
    plt.figure(n) # plot results
    plt.plot(x.value[0:t_bal_node+1],y.value[0:t_bal_node+1],'g-',label=r'Ascent Phase {}'.format(nt))
    plt.plot(x.value[t_bal_node::],y.value[t_bal_node::],'y-',label=r'Ballistic Phase'.format(nt))
    plt.legend(loc='best')
    plt.xlabel('x (m)',fontsize=15)
    plt.ylabel('y (m)',fontsize=15)
    plt.xlim([0, x[-1]])
    plt.ylim([0, max(y)])
    plt.title("Trajectory for Range " +"of {Range} meters".format(Range=round(x_final,1)))
    plt.grid(True)

    plt.figure(n+1) # plot results
    plt.plot(m.time[0:t_bal_node+1],v.value[0:t_bal_node+1],'m-',label=r"Ascent Phase")
    plt.plot(m.time[t_bal_node::],v.value[t_bal_node::],color="yellowgreen",label=r"Ballistic Phase")
    plt.title("Velocity for Range " +"of {Range} meters".format(Range=round(x_final,1)))
    plt.legend(loc='best')
    plt.xlabel('Time(s)',fontsize=15)
    plt.ylabel('V(m/s)',fontsize=15)
    plt.grid(True)
    
    plt.figure(n+2) # plot results
    plt.plot(m.time[0:t_bal_node+1],np.multiply([x for x in v.value],[np.sin(x) for x in gamma.value])[0:t_bal_node+1],'m-',label=r"Ascent Phase")
    plt.plot(m.time[t_bal_node::],np.multiply([x for x in v.value],[np.sin(x) for x in gamma.value])[t_bal_node::],color="yellowgreen",label=r"Ballistic Phase")
    plt.legend(loc='best')
    plt.title("Vertical Velocity for Range " +"of {Range} meters".format(Range=round(x_final,1)))
    plt.xlabel('Time(s)',fontsize=15)
    plt.ylabel(r'$V_y$ (m/s)',fontsize=15)
    plt.grid(True)

    
    plt.figure(n+3) # plot results
    plt.plot(m.time[0:t_bal_node+1],[x * 180/np.pi for x in gamma.value][0:t_bal_node+1],'b-',label=r"Ascent Phase")
    plt.plot(m.time[t_bal_node::],[x * 180/np.pi for x in gamma.value][t_bal_node::],color="orange",label=r"Ballistic Phase")
    plt.title("Flight Path Angle for Range " +"of {Range} meters".format(Range=round(x_final,1)))
    plt.legend(loc='best')
    plt.xlabel('Time(s)',fontsize=15)
    plt.ylabel(r'$\gamma (\degree)$',fontsize=15)    
    plt.grid(True)
    
    plt.figure(n+4)
    plt.plot(m.time[0:t_bal_node+1],[100*x/m_dot_max for x in m_dot.value][0:t_bal_node+1],'r--',label=r"Ascent Phase")
    plt.plot(m.time[t_bal_node::],[100*x/m_dot_max for x in m_dot.value][t_bal_node::],"g--",label=r"Ballistic Phase")
    plt.title("Percentage of Mass Flow rate for Range " +"of {Range} meters".format(Range=round(x_final,1)))
    plt.legend(loc='best')
    plt.xlabel('Time(s)',fontsize=15)
    plt.ylabel(r'$m_{dot}/m_{dot_{max}} (\%)$',fontsize=15)
    plt.grid(True)

    
    plt.figure(n+5)
    # plt.plot(m.time,mass.value,'y',label="Mass")
    plt.plot(m.time[0:t_bal_node+1],mass.value[0:t_bal_node+1],'b-',label=r"Ascent Phase")
    plt.plot(m.time[t_bal_node::],mass.value[t_bal_node::],color="violet",label=r"Ballistic Phase")
    plt.title("Mass for Range " +"of {Range} meters".format(Range=round(x_final,1)))
    plt.legend(loc='best')
    plt.xlabel('Time(s)',fontsize=15)
    plt.ylabel('m (kg)',fontsize=15)
    plt.grid(True)

    
    plt.figure(n+6)
    plt.plot(m.time[0:t_bal_node+1],[x * 180/np.pi for x in alpha.value][0:t_bal_node+1],'b-',label=r"Ascent Phase")
    plt.plot(m.time[t_bal_node::],[x * 180/np.pi for x in alpha.value][t_bal_node::],color="orange",label=r"Ballistic Phase")
    plt.title("Angle of attack for Range " +"of {Range} meters".format(Range=round(x_final,1)))    
    plt.legend(loc='best')
    plt.xlabel('Time(s)',fontsize=15)
    plt.ylabel(r'$\alpha (\degree)$ ',fontsize=15)
    plt.grid(True)

    
    plt.figure(n+7)
    # plt.plot(m.time,np.add([x * 180/np.pi for x in alpha.value],[x * 180/np.pi for x in gamma.value]),'b-',label=r"Pitch Angle")
    plt.plot(m.time[0:t_bal_node+1],[x * 180/np.pi for x in theta.value][0:t_bal_node+1],'b-',label=r"Ascent Phase")
    plt.plot(m.time[t_bal_node::],[x * 180/np.pi for x in theta.value][t_bal_node::],color="orange",label=r"Ballistic Phase")
    plt.title("Pitch Angle for Range " +"of {Range} meters".format(Range=round(x_final,1)))
    plt.legend(loc='best')
    plt.xlabel('Time(s)',fontsize=15)
    plt.ylabel(r'$\theta (\degree)$',fontsize=15)
    plt.grid(True)

    
    plt.figure(n+8)
    # plt.plot(m.time,[x /Jz for x in M_rw.value],'r--',label=r'$M_{rw}$ / $J_Z$')
    plt.plot(m.time[0:t_bal_node+1],[x  for x in M_rw_J.value][0:t_bal_node+1],'r--',label=r"Ascent Phase")
    plt.plot(m.time[t_bal_node::],[x  for x in M_rw_J.value][t_bal_node::],"g--",label=r"Ballistic Phase")
    plt.title("Moment from reaction wheels as a fraction of $J_Z$ for Range " +"of {Range} meters".format(Range=round(x_final,1)))
    plt.legend(loc='best')
    plt.xlabel('Time(s)',fontsize=15)
    plt.ylabel(r'$M_{rw}/J_Z$',fontsize=15)
    plt.grid(True)

    
    plt.figure(n+9) # plot results
    # plt.plot(m.time,np.multiply([x for x in v.value],[np.cos(x) for x in gamma.value]),'m-',label=r'Horizontal Velocity')
    plt.plot(m.time[0:t_bal_node+1],np.multiply([x for x in v.value],[np.cos(x) for x in gamma.value])[0:t_bal_node+1],'m-',label=r"Ascent Phase")
    plt.plot(m.time[t_bal_node::],np.multiply([x for x in v.value],[np.cos(x) for x in gamma.value])[t_bal_node::],color="yellowgreen",label=r"Ballistic Phase")
    plt.title("Horizontal Velocity for Range " +"of {Range} meters".format(Range=round(x_final,1)))    
    plt.legend(loc='best')
    plt.xlabel('Time(s)',fontsize=15)
    plt.ylabel(r'$V_{x} (m/s)$',fontsize=15)
    plt.grid(True)

    print(mass.value[-1])
    n+=10
    # print(M_rw.VALUE)
    
def main():	
    step=0.1 #Guess for the step
    t_final=Initial_Guess.Initial_guess_Hop(x_final,twr)[2]
    # t_final=5.5
    print(t_final)
    interval=10
    nt_values= np.linspace(int(t_final/step)-interval,int(t_final/step)+interval,interval*2+1) #Try for a range of nt to see which works
    # nt_values=[40,41,42,43,44,45,46,47,48,49,50] #Insert here the values of nt that work 
    nt_values=[45] #Insert here the values of nt that work 
    print(nt_values)
    for i in nt_values:
        optimize_trajectory(int(i),t_final,x_final)

if __name__ == "__main__":
    main()

#Do not forget to change the condition of the horizontal velocity

##Nt_values for different ranges, for a max final horizontal speed of 1 m/s
#4 -> 36
#3.2 -> 59
#2.501 -> 33,49
#2 -> 32
#1.499 -> 28
#1 -> 61

##Nt_values for different ranges, for a max final horizontal speed of 0 m/s
#3.2 -> 51,53 (t=5.5)
#3.199 -> 54 (Best one, t=5.5) 