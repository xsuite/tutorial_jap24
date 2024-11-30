from cpymad.madx import Madx
import matplotlib.pyplot as plt
import numpy as np

plt.ion()
plt.close('all')

madx = Madx()
madx.call("sps/sps.seq")
madx.call("sps/strengths/lhc_q20.str")
madx.call("sps/toolkit/macro.madx")
madx.beam(particle='proton',pc=26)
madx.use(sequence="sps")

qx0=20.13
qy0=20.18
madx.input(f"exec, sps_match_tunes({qx0},{qy0});")

# madx.input("select, flag=error, clear; select, flag=error, pattern=qd.11710; ealign, dx=-0.0048;")
# madx.input("select, flag=error, clear; select, flag=error, pattern=qf.11810; ealign, dx=-0.00297;")
# madx.input("select, flag=error, clear; select, flag=error, pattern=qda.11910; ealign, dx=-0.0048;")

madx.input("savebeta,sequence=sps,label=betastart,place=#S;")
tw = madx.twiss(pt=0.5e-3)
dx = tw['dx']

tw = madx.twiss()

f,ax = plt.subplots(1)
# ax.plot(tw0['s'],tw0['dx'])
# ax.plot(tw0['s'],tw0['dx']/np.sqrt(tw0['betx']))
ax.plot(tw['s'],tw['dx']-dx)
# ax.plot(tw['s'],tw['x'])
# ax.plot(tw['s'],tw['y'])

dx_norm = dx/np.sqrt(tw['betx'])
# ax.plot(tw['s'],tw['x'])

""
madx.input("select,FLAG=ERROR,RANGE=QF.11810/QF.31810,CLASS=RBEND,PATTERN=MB*;")
# madx.input("EFCOMP, DKN:={-10e-6};")
""
# tw = madx.twiss()

# ax.plot(tw['s'],tw['dx'])
# ax.plot(tw['s'],tw['dx']/np.sqrt(tw['betx']))
# ax.plot(tw['s'],tw['dx']/np.sqrt(tw['betx'])-dx_norm) 
# ax.plot(tw['s'],tw['dx']-dx) 


madx.input('actcse.31632,volt=4.5,lag=0.75,freq=400;')


tw1 = madx.twiss(beta0="betastart") #,px=9.05684792e-05)
f,ax = plt.subplots(1)
# ax.plot(tw1['s'],tw1['x']*1e3)
i = [i for i,n in enumerate(tw.name) if 'bph' in n]
ax.bar(tw1['s'][i],tw1['x'][i]*1e3,60,color='g')
ax.set_ylim(-7,7)

# f,ax = plt.subplots(1)
# ax.plot(tw1['s'],tw1['pt'])

# f,ax = plt.subplots(1)
# ax.plot(tw1['s'],tw1['t'])
