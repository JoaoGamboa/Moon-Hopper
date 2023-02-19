import numpy as np
import matplotlib.pyplot as plt  
import os
import pandas as pd

#Initial constants 
g0=9.81 #Earth's Gravity
gm=1.62 #Moon's Gravity

path=os.getcwd()+"\\Plots\\"+"Data for Range of 3.2.csv"
name="df32"
command=pd.read_csv("{}".format(path))
exec(name+"="+"command")   #read the csv file (put 'r' before the path string to address any special characters in the path, such as '\'). Don't forget to put the file name at the end of the path + ".csv"


#Parameters to be changed
m_prop_m_total=(100-df32.iloc[-1,4])/100 #Percentage of the mass needed for each hop (0.0022 for normal method, 0.0033 for alternative)
D_covered=200
R_hop=3.2

def subsystem_size(pl_percentage):
    """
    This function calculates the mass of the subsystems of the hopp
    
    Inputs:
        pl_percentage: percentage of the total mass that is payload (8% to 13%)
        
    Outputs:
        m_total_2-m_total_1:difference between initial mass and the mass after the iteration
        m_total_1:initial guess of the total mass
        m_total_2:Mass after the iteration        
    """    
    
    m_pl=3.3
    m_total_1=m_pl/pl_percentage #The pl_percentage change from 8% to 12%
    
    m_propellant=m_total_1-m_total_1*((1-m_prop_m_total)**(np.ceil(D_covered/R_hop))) #Each hop uses 0.3% in the worst case, the descent into the pit uses 1.33% of the total mass
    print(m_propellant)
    m_dry=m_total_1-m_propellant-m_pl
    m_inert= m_dry+m_pl

    #Estimation of Thrust Subsystem -> Taken from MIT Talaris Propulsion sizing
    m_e= 0.1*(1.3*gm*m_total_1)**(2/3) #Mass of the engine
    m_tank=2/3*m_propellant**(2/3)
    m_prop_system=(m_e+m_propellant+m_tank) * 1.1
    m_structural=3/20*m_total_1 * 1.1
    #1.1 is a safe margin
    
    m_power= ((23/0.22) * 2 )/100 * 1.1 #power from payload/ 0.22 (smad) /density + safe margin table 14.20
    
    m_mobility_system=5.3/30.9*m_dry *1.1 #Amalia rover -> change to percentage + 10%
     
    #SMAD Appendix A 
    m_cdh=0.04*m_dry *1.1
    m_ttc=0.07*m_dry *1.1
    m_thermal=0.06*m_dry *1.1
    
    # m_comms_avionics_power_thermal=32.03/91.91*m_total #Percentage of the MIT hopper
    
    m_other= 0.04*m_dry *1.1 #Paper -> A new sizing methodology

    m_total_2=m_prop_system+m_structural+m_power+m_mobility_system+m_pl+m_cdh+m_ttc+m_thermal+m_other
    # print(m_total_1)
    # print(m_total_2)
    return m_total_2-m_total_1, m_total_1,m_total_2

print(subsystem_size(.1))
diff=[]
pl_percentage=np.arange(0.10,0.17,0.01)
for i in pl_percentage:
    i=round(i,2)
    diff.append(subsystem_size(i)[0])
    print("Difference between the masses of {} kg for a payload percentage of {}% and an initial mass of {} kg".format(round(subsystem_size(i)[0],2),int(100*i),round(subsystem_size(i)[1],1)))
    

plt.plot([100*x for x in pl_percentage],[abs(subsystem_size(x)[0]) for x in pl_percentage],color="brown")
plt.title("Convergence of system subsizing")
plt.grid(True)
plt.xlabel('Payload percentage (%)',fontsize=15)
plt.ylabel(r'Mass Difference',fontsize=15)    
plt.savefig(os.getcwd()+"\\plots\\" + "system_subsizing_convergence.png")


# plt.plot([x for x in np.arange(0.08,0.13,0.001)],[abs(subsystem_size(x)[0]) for x in np.arange(0.08,0.13,0.001)],color="brown")
# plt.title("Convergence of system subsizing")
# plt.grid(True)
# plt.xlabel('Payload percentage (decimal)',fontsize=15)
# plt.ylabel(r'Mass Difference',fontsize=15)    

