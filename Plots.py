# -*- coding: utf-8 -*-
"""
Created on Sat Jun 18 14:17:17 2022

@author: Joao Gamboa
"""

import pandas as pd
import os
import matplotlib.pyplot as plt  
import numpy as np
import Initial_Guess

twr=1.3
m_total=100

for i in [1,1.5,2,2.5,3.2]:
    path=os.getcwd()+"\\Plots\\"+"Data for Range of {}.csv".format(i)
    name="df{}".format(int(i*10))
    command=pd.read_csv("{}".format(path))
    exec(name+"="+"command")   #read the csv file (put 'r' before the path string to address any special characters in the path, such as '\'). Don't forget to put the file name at the end of the path + ".csv"


#%%

plt.plot(df10.iloc[:,1],df10.iloc[:,2],'g-',label=r'$Range=1 m$')
plt.plot(df15.iloc[:,1],df15.iloc[:,2],'y-',label=r'$Range=1.5 m$')
plt.plot(df20.iloc[:,1],df20.iloc[:,2],'r-',label=r'$Range=2 m$')
plt.plot(df25.iloc[:,1],df25.iloc[:,2],'b-',label=r'$Range=2.5 m$')
plt.plot(df32.iloc[:,1],df32.iloc[:,2],'m-',label=r'$Range=3.2 m$')
plt.legend(loc='best')
plt.grid(True)
plt.xlabel('x(m)')
plt.ylabel('y(m)')
plt.savefig(os.getcwd()+"\\plots\\" + "All_trajectories.png")


#%%
plt.figure(2)
plt.plot([1,1.5,2,2.5,3.2],[100*(m_total-df10.iloc[-1,4])/m_total,100*(m_total-df15.iloc[-1,4])/m_total,100*(m_total-df20.iloc[-1,4])/m_total,100*(m_total-df25.iloc[-1,4])/m_total,100*(m_total-df32.iloc[-1,4])/m_total],color="orange",label="Optimized")
plt.plot([1,1.5,2,2.5,3.2],[100*Initial_Guess.Initial_guess_Hop(x_final, twr)[7] for x_final in [1,1.5,2,2.5,3.2]],'b-',label="Initial Guess")
plt.legend(loc='best')
plt.grid(True)
plt.ylabel(r'$\frac{|m(t_f)-m_0|}{m_0} (\%)$',fontsize=15)
plt.xlabel('$x_{final} (m)$',fontsize=15)
plt.savefig(os.getcwd()+"\\plots\\" + "Mass_all_traj.png")
print("Mass of propellant over the total mass for a range of 3.2m is {}",100*(m_total-df32.iloc[-1,4])/m_total)


#%%

deg2rad= np.pi/180 
plt.plot(df10.iloc[:,0],df10.iloc[:,3]/deg2rad,'g-',label=r'$Range=1 m$')
plt.plot(df15.iloc[:,0],df15.iloc[:,3]/deg2rad,'y-',label=r'$Range=1.5 m$')
plt.plot(df20.iloc[:,0],df20.iloc[:,3]/deg2rad,'r-',label=r'$Range=2 m$')
plt.plot(df25.iloc[:,0],df25.iloc[:,3]/deg2rad,'b-',label=r'$Range=2.5 m$')
plt.plot(df32.iloc[:,0],df32.iloc[:,3]/deg2rad,'m-',label=r'$Range=3.2 m$')
plt.title("Pitch angle for all the trajectories")
plt.grid(True)
plt.legend(loc='best')
plt.ylabel(r'$\theta (\degree)$')
plt.xlabel('Time (s)')
plt.savefig(os.getcwd()+"\\plots\\" + "Pitch Angle all traj.png")

#%%

plt.plot(df10.iloc[:,0],df10.iloc[:,5],'g-',label=r'$Range=1 m$')
plt.plot(df15.iloc[:,0],df15.iloc[:,5],'y-',label=r'$Range=1.5 m$')
plt.plot(df20.iloc[:,0],df20.iloc[:,5],'r-',label=r'$Range=2 m$')
plt.plot(df25.iloc[:,0],df25.iloc[:,5],'b-',label=r'$Range=2.5 m$')
plt.plot(df32.iloc[:,0],df32.iloc[:,5],'m-',label=r'$Range=3.2 m$')
plt.title("Profile of the Reaction Wheels for all the trajectories")
plt.grid(True)
plt.legend(loc='best')
plt.ylabel(r'$M_{rw} \hspace{0.5} /\hspace{0.5}Jz$')
plt.xlabel('Time (s)')
plt.savefig(os.getcwd()+"\\plots\\" + "M_rw all traj.png")

#%%

plt.plot(df10.iloc[:,0],df10.iloc[:,6],'g-',label=r'$Range=1 m$')
plt.plot(df15.iloc[:,0],df15.iloc[:,6],'y-',label=r'$Range=1.5 m$')
plt.plot(df20.iloc[:,0],df20.iloc[:,6],'r-',label=r'$Range=2 m$')
plt.plot(df25.iloc[:,0],df25.iloc[:,6],'b-',label=r'$Range=2.5 m$')
plt.plot(df32.iloc[:,0],df32.iloc[:,6],'m-',label=r'$Range=3.2 m$')
plt.title("Mass flow rate profile for all the trajectories")
plt.grid(True)
plt.legend(loc='best')
plt.ylabel(r'$m_{dot}$ (kg/s)')
plt.xlabel('Time (s)')
plt.savefig(os.getcwd()+"\\plots\\" + "m_dot all traj.png")

#%%
#Mass of propellant consumed as a function of the twr 
twr_val=np.arange(1,5,0.1)
#The value of $\frac{|m(t_f)-m_0|}{m_0}(\%)$ is equal to 100* twr*g/(Isp*g0)* t_bal
plt.figure(1)
plt.plot(twr_val,[100*Initial_Guess.Initial_guess_Hop(3.2,i)[7] for i in twr_val],'r-')
plt.grid(True)
plt.ylabel(r'$\frac{|m(t_f)-m_0|}{m_0}(\%)$',fontsize=15)
plt.xlabel('$twr$',fontsize=15)
plt.title("Mass of propellant used for a hop of 3.2 m")

#%%
#Final horizontal velocity as a function of the range of the hop
plt.figure(2)
plt.plot(twr_val,[Initial_Guess.Initial_guess_Hop(3.2,i)[0]*np.cos(Initial_Guess.Initial_guess_Hop(3.2,i)[4]*deg2rad) for i in twr_val],'r-')
plt.grid(True)
plt.title("Initial Guess for the final horizontal velocity")
plt.ylabel(r'$V_{x_f} (m/s)$',fontsize=15)
plt.xlabel('$twr(m)$',fontsize=15)


#%%
#Final horizontal velocity as a function of the range of the hop
plt.figure(2)
distance=np.arange(1,3.3,0.1)
twr=1.3
plt.plot(distance,[Initial_Guess.Initial_guess_Hop(i,twr)[0]*np.cos(Initial_Guess.Initial_guess_Hop(i,twr)[4]*deg2rad) for i in distance],'r-')
plt.grid(True)
plt.title("Initial Guess for the final horizontal velocity")
plt.ylabel(r'$V_{x_f} (m/s)$',fontsize=15)
plt.xlabel('$x_{final} (m)$',fontsize=15)

#%%
plt.figure(3)
t_bal_node_IG=np.where(np.array([round(i,2) for i in Initial_Guess.Initial_guess_Hop(3.2,twr)[6]])==round(Initial_Guess.Initial_guess_Hop(3.2,twr)[1],2))[0][0]#Point where y=y_ball
t_bal_node_opt=28 #Node taken from the Opt_Hop.py file
plt.plot(Initial_Guess.Initial_guess_Hop(3.2,twr)[5][0:t_bal_node_IG+1],Initial_Guess.Initial_guess_Hop(3.2,twr)[6][0:t_bal_node_IG+1],'b-',label=r'Initial Guess Ascent Phase')
plt.plot(Initial_Guess.Initial_guess_Hop(3.2,twr)[5][t_bal_node_IG::],Initial_Guess.Initial_guess_Hop(3.2,twr)[6][t_bal_node_IG::],color="orange",label=r'Initial Guess Ballistic Phase')
plt.plot(df32.iloc[:,1][0:t_bal_node_opt+1],df32.iloc[:,2][0:t_bal_node_opt+1],'g-',label=r'Optimized Ascent Phase')
plt.plot(df32.iloc[:,1][t_bal_node_opt::],df32.iloc[:,2][t_bal_node_opt::],'y-',label=r'Optimized Ballistic Phase')
plt.legend(loc='best')
plt.grid(True)
plt.title("Initial Guess for trajectory for Range of " +"{Range} meters".format(Range=round(x_final,1)))
plt.ylabel('$y(m)$',fontsize=15)
plt.xlabel('$x (m)$',fontsize=15)
plt.savefig(os.getcwd()+"\\plots\\" + "IG_opt_traj.png")

