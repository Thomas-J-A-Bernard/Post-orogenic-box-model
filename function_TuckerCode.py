"""
Created on Wed Mar 27 11:12:05 2019
@author: MarK Naylor (Python 3.6)
Description:
Original code of Tucker and van der Beek, (2013)
"""

import numpy as np
import math as mt

#-----------------------------------------------------------------------------#
def Tucker_Code(dt, nt, tstart, trt, tbt, mut, Tc, rhos, rhoc, rhoi, Lr, Lb, f, Vt, epeirt, epeibt, Hr, Hb, etar, etab, wr, Psir, wb, Psib, Fe, Fb, Fr):  
  for i in range(1, int(nt)):
    ## Increment current time
    #t = tstart + i*dt
    
    ## Switch off tectonics after a while
    #if t>tectonic_duration:
    #    mu = 0
    
    ## Time variation option: update values for new time step
    #if (opt_time_variation):
    mu = mut[i-1]
    tr = trt[i-1]
    tb = tbt[i-1]
    Kra = 2*(mu/Tc)*(rhos/rhoc)*f     # range accretion coefficient
    Kbs = (rhoc/rhos)*(Lr/Lb)*(1/tr)  # sed flux to basin coefficient
    Kba = 2*(mu/Tc)*(Lr/Lb)*f         # basin accretion loss coefficient
    
    epeir = epeirt[i-1]
    epeib = epeibt[i-1]
    
    ## Epeirogeny: update cumulative "deep" uplift
    cumepeir = 0.0
    cumepeib = 0.0
    cumepeir = cumepeir + epeir*dt
    cumepeib = cumepeib + epeib*dt    
    
    ## Compute time derivatives
    #detardt = mu + Kra*etab[i-1] - (1/tr)*((1-Psir)*etar[i-1]-(etab[i-1]-Psib*etar[i-1]))
    detardt = mu + Kra*etab[i-1] - (1/tr)*((Hr[i-1] - Hb[i-1]))
    #detabdt = Kbs*((1-Psir)*etar[i-1]-(etab[i-1]-Psib*etar[i-1])) - Kba*etab[i-1] - (1/tb)*(etab[i-1]-Psib*etar[i-1])
    detabdt = Kbs*((Hr[i-1]-Hb[i-1])) - Kba*etab[i-1] - (1/tb)*(Hb[i-1])
    
    ## Update new values
    etar[i] = etar[i-1] + detardt*dt              # Update range thickness
    etab[i] = etab[i-1] + detabdt*dt              # Update basin thickness
    wr[i] = Psir*etar[i]                          # Update Flexure beneath range
    wb[i] = Psib*etar[i]                          # Update Flexure beneath basin
    Hr[i] = etar[i] - wr[i] #+ cumepeir           # Update range height
    Hb[i] = etab[i] - wb[i] #+ cumepeib           # Update basin height
    Fe[i] = ((Hr[i]-Hb[i])*rhoc*Lr)/trt[i]*dt     # Update sedimentary flux from the range to the basin
    Fb[i] = (Hb[i]*rhoi*Lb)/tbt[i]*dt             # Update sediemntary flux from the basin to the base level
    Fr[i] = Vt[i]*Tc*rhoc*dt                      # Update tectonic flux to the range
    
    #wxr = (rhoc*etar[i]/(2*(rhom-rhoi)))*(2-np.exp(-Lambda*b_range)*np.cos(Lambda*b_range)-np.exp(-Lambda*a_range)*np.cos(Lambda*a_range))
    #wxb = -(rhoc*etar[i]/(2*(rhom-rhoi)))*(np.exp(-Lambda*a_basin)*np.cos(Lambda*a_basin)-np.exp(-Lambda*b_basin)*np.cos(Lambda*b_basin))
    
  tt = tstart + dt*np.arange(0,nt);
  tt = tt/1.e6;
  
  return Hr, Hb, etar, etab, wr, wb, Fe, Fb, Fr, tt

#-----------------------------------------------------------------------------#
def Update_tectonic_parameter(nt, T1_time_steps, T2_time_steps, T3_time_steps, T4_time_steps, T5_time_steps, T6_time_steps, T7_time_steps, Tc, Lt, V1, V2, V3, V4, V5):
    
    Vt = np.zeros(nt)
    Vt[T1_time_steps:T2_time_steps] = V1
    Vt[T2_time_steps:T3_time_steps] = V2
    Vt[T3_time_steps:T4_time_steps] = V3
    Vt[T4_time_steps:T5_time_steps] = V4
    Vt[T5_time_steps:T7_time_steps] = V5
  
    V_T5 = V4
    for v in range(T5_time_steps, T6_time_steps):
        Vt[v] = V_T5 - (V4-V5)/(T6_time_steps-T5_time_steps)
        V_T5 = Vt[v]
       
    mut = Vt*Tc/Lt
    
    return Vt, mut

#-----------------------------------------------------------------------------#
def Update_flexural_parameter(Te, E, nu, g, Lr, rhom, rhoi, rhoc):
    
    D = E*Te**3/(12*(1-nu**2))                  # Flexural rigidity, Nm
    Lambda = ((rhom-rhoi)*g/(4*D))**(1./4.)     # inverse flexural parameter, L-1
    tlaml = 2*Lambda*Lr                         # density ratio
    drat = rhoc/(rhom-rhoi)                     # 2xlambdaxLr
    
    Psir = drat*np.exp(-tlaml)*(np.exp(tlaml)*(-1+2*tlaml)+np.cos(tlaml)-np.sin(tlaml))/(2*tlaml)                                                                                               # average deflexion beneath the range
    
    Lb = 3*np.pi/(4*Lambda)-Lr                  # basin length, m
    Llm = Lambda*Lr
    pi34 = 3*np.pi/4              
     
    Psib = (rhoc/(2*(rhom-rhoi)))*(1/(-4*Llm+3*np.pi))*2*np.exp(-2*Llm-pi34)*(np.sqrt(2)*np.exp(Llm)*(-1+np.exp(2*Llm))*np.cos(Llm)+np.exp(pi34)*(np.exp(2*Llm)-np.cos(2*Llm)+np.sin(2*Llm)))   # average deflexion neneath the basin
    
    return Psir, Psib, Lb, Lambda, pi34

#-----------------------------------------------------------------------------#
def Update_response_time_parameter(nt, Lr, Lb, kappar, kappab, T1_time_steps, T7_time_steps):
    
    tr = (Lr**2)/kappar; tb = (Lb**2)/kappab        # range and basin response time                          
    trt = np.zeros(nt); tbt = np.zeros(nt)
    trt[T1_time_steps:T7_time_steps] = tr
    tbt[T1_time_steps:T7_time_steps] = tb
    
    return tr, tb, trt, tbt

#-----------------------------------------------------------------------------#
def Update_epeiorogenic_parameter(nt, pi34, Lr, Lambda, epei_crest, epeistart, epeiend):
    
    Lrp = Lr*Lambda                                                                             # dimensionless range width (scaled by flex param)
    epeir = epei_crest*(1/Lrp)*(1-np.exp(-Lrp)*np.cos(Lrp))
    epeib = epei_crest*(1/(pi34-Lrp))*(np.exp(-Lrp)*np.cos(Lrp)-np.exp(-pi34)*np.cos(pi34))

    epeirt = np.zeros(nt)
    epeibt = np.zeros(nt)
    epeirt[epeistart:epeiend] = epeir
    epeibt[epeistart:epeiend] = epeib
    
    return epeirt, epeibt