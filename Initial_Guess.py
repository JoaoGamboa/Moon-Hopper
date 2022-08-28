import numpy as np
from scipy.optimize import fsolve
from scipy import integrate
import matplotlib.pyplot as plt

#Initial constants 
g0=9.81 #Earth's Gravity
gm=1.62 #Moon's Gravity

deg2rad=np.pi/180
def Initial_guess_Hop(x_final,twr):
    """
    This function computes the initial guess for the hoppping trajectory
    
    Inputs:
        x_final: the range of the hop
	twr: Thrust-to-Weight ration
    
    Outputs:
        t_final: the time the trajectory takes     
	v_ballistic: the velocity in the end of the ascent phase
	y_ballistic: the y coordinate at the end of the ascent phase 
	t_final: The time the trajectory takes to complete
	t_ballistic: The time of the ascent phase
	gamma_ballistic/deg2rad: The flight path angle at the end of the ascent phase, in degrees
	x_coord: The array containing all the x coordinate points of the trajectory
	y_coord: The array containing all the x coordinate points of the trajectory   
    """
    A=twr*gm
    V_y_final=2.5 #Final vertical velocity in absolute module
       
    gamma0=90*deg2rad #Initial gamma
    #The five equations that define the system 
    gamma=lambda t,t_ballistic,gamma_ballistic: gamma_ballistic+((gamma0-gamma_ballistic)/(t_ballistic**2))*(t-t_ballistic)**2
    eq1 =lambda v_ballistic,t_ballistic,gamma_ballistic: integrate.quad(lambda t: A-gm*np.sin(gamma(t,t_ballistic,gamma_ballistic)),0,t_ballistic)[0]-v_ballistic
    eq2 =lambda y_ballistic,t_ballistic,gamma_ballistic: integrate.quad(lambda t: integrate.quad(lambda t: A-gm*np.sin(gamma(t,t_ballistic,gamma_ballistic)),0,t)[0]*np.sin(gamma(t,t_ballistic,gamma_ballistic)),0,t_ballistic)[0]-y_ballistic
    eq3=lambda v_ballistic,t_final,t_ballistic,gamma_ballistic: v_ballistic*np.sin(gamma_ballistic)-gm*(t_final-t_ballistic)+V_y_final
    eq4=lambda t_ballistic,y_ballistic,v_ballistic,t_final,gamma_ballistic: (v_ballistic*np.sin(gamma_ballistic))/gm * (1+np.sqrt(1+2*(y_ballistic)*gm/(v_ballistic*np.sin(gamma_ballistic))**2))+t_ballistic-t_final
    eq5= lambda v_ballistic,gamma_ballistic,t_final,t_ballistic: v_ballistic*np.cos(gamma_ballistic)*(t_final-t_ballistic)+integrate.quad(lambda t: integrate.quad(lambda t: A-gm*np.sin(gamma(t,t_ballistic,gamma_ballistic)),0,t)[0]*np.cos(gamma(t,t_ballistic,gamma_ballistic)),0,t_ballistic)[0]-x_final

    #Solving the equations
    equations=lambda v_ballistic,y_ballistic,t_final,t_ballistic,gamma_ballistic: (eq1(v_ballistic,t_ballistic,gamma_ballistic),eq2(y_ballistic,t_ballistic,gamma_ballistic),eq3(v_ballistic,t_final,t_ballistic,gamma_ballistic),eq4(t_ballistic,y_ballistic,v_ballistic,t_final,gamma_ballistic),eq5(v_ballistic,gamma_ballistic,t_final,t_ballistic))
    equations_2=lambda x :equations(x[0],x[1],x[2],x[3],x[4])
    solution=fsolve(equations_2,(2,0.1,5,0.4,75*deg2rad))
    v_ballistic,y_ballistic,t_final,t_ballistic,gamma_ballistic=solution[0],solution[1],solution[2],solution[3],solution[4]
    
    
    #Calculating the final horizontal velocity
    # print(v_ballistic*np.cos(gamma_ballistic))
    
    #Rounding the values 
    x_final=round(x_final,2)
    t_final=round(t_final,2)
    
    step=0.01

    #Trajectory calculated from the initial guess subscript 1 denotes the ascent phase and 2 the ballistic Phase
    x_coord_1=[integrate.quad(lambda t: integrate.quad(lambda t: A-gm*np.sin(gamma(t,t_ballistic,gamma_ballistic)),0,t)[0]*np.cos(gamma(t,t_ballistic,gamma_ballistic)),0,t)[0] for t in np.arange(0,t_ballistic,step)]
    x_coord_2=[x_coord_1[-1]+v_ballistic*np.cos(gamma_ballistic)*(t-t_ballistic) for t in np.arange(t_ballistic,t_final,step)]

    y_coord_1=[integrate.quad(lambda t: integrate.quad(lambda t: A-gm*np.sin(gamma(t,t_ballistic,gamma_ballistic)),0,t)[0]*np.sin(gamma(t,t_ballistic,gamma_ballistic)),0,t)[0] for t in np.arange(0,t_ballistic,step)]
    y_coord_2=[y_coord_1[-1]+v_ballistic*np.sin(gamma_ballistic)*(t-t_ballistic)-0.5*gm*(t-t_ballistic)**2 for t in np.arange(t_ballistic,t_final,step)]

    x_coord=x_coord_1+x_coord_2
    y_coord=y_coord_1+y_coord_2
    # print(len(x_coord))
    return v_ballistic,y_ballistic,t_final,t_ballistic,gamma_ballistic/deg2rad,x_coord,y_coord

def main():
    twr=1.3
    x_final=3.2
    Initial_guess_Hop(x_final,twr)	