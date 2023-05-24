#%%
#Libraries
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams.update({'font.size': 13, 'lines.linewidth': 2.5})
from matplotlib.widgets import Slider, Button

#%%
#Explicit equations
UA=16500;
FA0=0.0376;
FI0=0;
CpA=163;
CpB=83;
CT0=18.8;
mc=0.111;
Cpcoolant=34.5;
CpC=71;
To=1035;
CpI=32.68;
Ta0=1250
Ea=284.52
A=8.198*10**14
def ODEfun(Yfuncvec, V, UA,FA0,FI0,CpA,CpB,CT0,mc,Cpcoolant,CpC,To,CpI,Ta0,Ea,A):
    Ta=Yfuncvec[0]
    T= Yfuncvec[1]
    X= Yfuncvec[2]
        #Explicit Equation Inline
    FT0=FI0+FA0;
    vo=(FT0/CT0)*1000;
    Ca0=(FA0/vo)*1000;
    delCp=CpB+CpC-CpA;
    thetaI=FI0/FA0;
    sumcp=CpA+thetaI*CpI;
    deltaH = 80770 + delCp * (T - 298); 
    k = A*np.exp(-Ea/(0.008314*T)); 
    ra = 0 - k*(Ca0 * (1 - X) * (To / T)) / (1 + X);
    Qg=ra * deltaH;
    Qr=UA*(T-Ta);
    # Differential equations
    dTadV = 0*UA * (T - Ta) / mc / Cpcoolant; 
    dTdV = (Qg - Qr) / (FA0 * (sumcp + X * delCp)); 
    dXdV = 0 - (ra / FA0); 
    return np.array([dTadV,dTdV,dXdV])

Vspan = np.linspace(0, 0.001, 1000) # Range for the independent variable
y0 = np.array([Ta0,To,0]) # Initial values for the dependent variables

#%%
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
fig.suptitle("""LEP-12-2b:Production of Acetic Anhydride (Constant Ta)""", fontweight='bold', x = 0.2,y=0.97)
plt.subplots_adjust(left  = 0.48)
fig.subplots_adjust(wspace=0.23,hspace=0.3)
sol = odeint(ODEfun, y0, Vspan, (UA,FA0,FI0,CpA,CpB,CT0,mc,Cpcoolant,CpC,To,CpI,Ta0,Ea,A))
Ta =sol[:, 0]
T = sol[:, 1]
X = sol[:, 2]
FT0=FI0+FA0;
vo=(FT0/CT0)*1000;
Ca0=(FA0/vo)*1000;
delCp=CpB+CpC-CpA;
thetaI=FI0/FA0;
sumcp=CpA+thetaI*CpI;
deltaH = 80770 + delCp * (T - 298); 
k = A*np.exp(-Ea/(0.008314*T)); 
ra = 0 - k*(Ca0 * (1 - X) * (To / T)) / (1 + X);
rate = 0 - ra;
Qg=ra * deltaH;
Qr=UA*(T-Ta);
p1,p2= ax1.plot(Vspan,T,Vspan,Ta)
ax1.legend(['T','$T_a$'], loc='upper right')
ax1.set_xlabel(r'$Volume  {(m^3)}$', fontsize='medium')
ax1.set_ylabel('Temperature (K)', fontsize='medium')
ax1.ticklabel_format(style='sci',scilimits=(0,0),axis='x')
ax1.set_ylim(900,1500)
ax1.set_xlim(0,0.001)
ax1.grid()

p3 = ax2.plot(Vspan,X)[0]
ax2.legend(['X'], loc='upper left')
ax2.set_ylim(0,1)
ax2.set_xlim(0,0.001)
ax2.grid()
ax2.set_xlabel(r'$Volume  {(m^3)}$', fontsize='medium')
ax2.set_ylabel('Conversion', fontsize='medium')
ax2.ticklabel_format(style='sci',scilimits=(0,0),axis='x')

p5,p6 = ax3.plot(Vspan,Qg,Vspan,Qr)
ax3.legend(['$Q_g$','$Q_r$'], loc='lower right')
#ax3.set_ylim(0,100)
ax3.set_xlim(0,0.001)
ax3.grid()
ax3.set_xlabel(r'$Volume  {(m^3)}$', fontsize='medium')
ax3.set_ylabel(r'$Q {(J/m^3.s)}$', fontsize='medium')
ax3.ticklabel_format(style='sci',scilimits=(3,4),axis='y')
ax3.ticklabel_format(style='sci',scilimits=(0,0),axis='x')

p4 = ax4.plot(Vspan, rate)[0]
ax4.legend(['$-r_A$'], loc='upper right')
ax4.set_ylim(0,100)
ax4.set_xlim(0,0.001)
ax4.grid()
ax4.set_xlabel(r'$Volume  {(m^3)}$', fontsize='medium')
ax4.set_ylabel('$Rate {(mol/m^3.s)}$', fontsize='medium')
ax4.ticklabel_format(style='sci',scilimits=(0,0),axis='x')
#ax1.axis('off')
ax1.text(-0.0025,100,'Note: While we used the expression k=$k_1$*exp(E/R*(1/$T_1$ - 1/$T_2$)) \n         in the textbook, in python we have to use k=A*exp(-E/RT) \n          in order to explore all the variables.',wrap = True, fontsize=13,
        bbox=dict(facecolor='none', edgecolor='red', pad=10))
ax1.text(-0.0024,320,'Differential Equations'
         '\n\n'
         r'$\dfrac{dT_a}{dV} = \dfrac{Ua*(T-T_a)}{m_c*C_{P_{coolant}}}*0 $'
         '\n'
         r'$\dfrac{dT}{dV} = \dfrac{r_{A}\Delta H_{Rx}-Ua*(T-T_a)}{F_{A0}.\left(\sum_{i}\theta_iC_{P_i}+X \Delta C_P\right)}$'
         '\n'
         r'$\dfrac{dX}{dV} = \dfrac{-r_{A}}{F_{A0}}$'
                  '\n \n'
                  
         'Explicit Equations'
                  '\n\n'  
         r'$A = 8.198*10^{14}\thinspace s^{-1}$'
              '\n'
         r'$F_{A0}=0.0376$'
         '\n'
         r'$F_{T0}=F_{A0}+F_{I0}$' 
         '\n'
          r'$v_{0}=\dfrac{F_{T0}}{C_{T0}}$' 
         '\n'
         r'$C_{A0}=\dfrac{F_{A0}}{v_{0}}$' 
         '\n'
         r'$\theta_I = \dfrac{F_{I0}}{F_{A0}}}$'
                  '\n'
         r'$\sum_{i}\theta_iC_{pi} = C_{P_A} +\theta_IC_{P_I} $'
         '\n'
         r'$\Delta H_{Rx}=80770+ \Delta C_P* (T-298) $'
         '\n'
         r'$\Delta C_P= C_{P_B}+C_{P_C}-C_{P_A}$'
         '\n'
          r'$k = A*exp\left(\dfrac{-E}{0.00831*T}\right)$'
         '\n'
         r'$r_A = -\dfrac{kC_{A0}(1-X)}{1+X} \dfrac{T_0}{T}$'
                  '\n'
         r'$rate = -r_A$'
                 '\n'
         r'$Q_g = r_{A}.\Delta H_{Rx}$'
                  '\n'
         r'$Q_r = Ua(T-T_a)$', ha='left', wrap = True, fontsize=12,
        bbox=dict(facecolor='none', edgecolor='black', pad=10.0), fontweight='bold')

axcolor = 'black'
ax_FA0 = plt.axes([0.24, 0.80, 0.15, 0.015], facecolor=axcolor)
ax_FI0 = plt.axes([0.24, 0.76, 0.15, 0.015], facecolor=axcolor)
ax_CpA = plt.axes([0.24, 0.72, 0.15, 0.015], facecolor=axcolor)
ax_CT0 = plt.axes([0.24, 0.68, 0.15, 0.015], facecolor=axcolor)
ax_To = plt.axes([0.24, 0.64, 0.15, 0.015], facecolor=axcolor)
ax_UA = plt.axes([0.24, 0.60, 0.15, 0.015], facecolor=axcolor)
ax_Ta0 = plt.axes([0.24, 0.56, 0.15, 0.015], facecolor=axcolor)
ax_Ea = plt.axes([0.24, 0.52, 0.15, 0.015], facecolor=axcolor)

sFA0 = Slider(ax_FA0, r'F$_{A0}$($\frac{mol}{s}$)', 0.01, 1, valinit=0.0376,valfmt='%1.3f')
sFI0 = Slider(ax_FI0, r'F$_{I0}$ ($\frac{mol}{s}$)', 0, 10, valinit=0,valfmt='%1.1f')
sCpA= Slider(ax_CpA, r'C$_{P_A}$ ($\frac{J}{mol.k}$)', 50, 500, valinit=163,valfmt='%1.1f')
sCT0 = Slider(ax_CT0,r'C$_{T0}$ ($\frac{mol}{m^3}$)', 10, 40, valinit= 18.8,valfmt='%1.1f')
sTo = Slider(ax_To,  r'T$_0$ ($K$)', 950, 1500, valinit=1035,valfmt='%1.0f')
sUA = Slider(ax_UA, r'Ua ($\frac{J}{mol^3.s.K}$)', 1000, 30000, valinit=16500,valfmt='%1.0f')
sTa0 = Slider(ax_Ta0, r'T$_a$ ($K$)', 900, 1400, valinit=1250,valfmt='%1.0f')
sEa= Slider(ax_Ea, r'$E$ ($\frac{kJ}{mol}$)', 250, 350, valinit=284.5,valfmt='%1.1f')

def update_plot2(val):
    FA0 = sFA0.val
    FI0 =sFI0.val
    CpA =sCpA.val
    CT0 =sCT0.val
    To = sTo.val
    UA = sUA.val
    Ta0 = sTa0.val
    y0 = np.array([Ta0,To,0])
    Ea= sEa.val
    sol = odeint(ODEfun, y0, Vspan, (UA,FA0,FI0,CpA,CpB,CT0,mc,Cpcoolant,CpC,To,CpI,Ta0,Ea,A))
    Ta = sol[:, 0]
    T = sol[:, 1]
    X = sol[:, 2]
    FT0=FI0+FA0;
    vo=(FT0/CT0)*1000;
    Ca0=(FA0/vo)*1000;
    delCp=CpB+CpC-CpA;
    deltaH = 80770 + delCp * (T - 298); 
    k = A*np.exp(-Ea/(0.008314*T)); 
    ra = 0 - k*(Ca0 * (1 - X) * (To / T)) / (1 + X);
    rate = 0 - ra;   
    Qg=ra * deltaH;
    Qr=UA*(T-Ta);
    p1.set_ydata(T)
    p2.set_ydata(Ta)
    p3.set_ydata(X)
    p4.set_ydata(rate)
    p5.set_ydata(Qg)
    p6.set_ydata(Qr)
    fig.canvas.draw_idle()

sFA0.on_changed(update_plot2)
sFI0.on_changed(update_plot2)
sCpA.on_changed(update_plot2)
sCT0.on_changed(update_plot2)
sTo.on_changed(update_plot2)
sUA.on_changed(update_plot2)
sTa0.on_changed(update_plot2)
sEa.on_changed(update_plot2)

resetax = plt.axes([0.26, 0.84, 0.09, 0.04])
button = Button(resetax, 'Reset variables', color='cornflowerblue', hovercolor='0.975')

def reset(event):
    sFA0.reset()
    sFI0.reset()
    sCpA.reset()
    sCT0.reset()
    sTo.reset()
    sUA.reset()
    sTa0.reset()
    sEa.reset()
button.on_clicked(reset)
    
