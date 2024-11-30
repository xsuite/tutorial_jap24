import xtrack as xt
from cpymad.madx import Madx
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

plt.ion()
plt.close('all')


def parseYASPfile(file):
    with open(file,"r") as fid:
        filelines = fid.readlines()
    bpm_names = []
    bpm_pos = []
    bpm_stat = []
    bpm_plane = []
    header_ended = False
    for i,line in enumerate(filelines):
        s = list(filter(None, line.replace('\n','').split(' ')))
        if s[0] == '*':
            s.pop(s.index('*'))
            col_name = s.index('NAME')
            col_pos = s.index('POS')
            col_stat = s.index('STATUS-TAG')
            col_plane = s.index('PLANE')
            header_ended = True
            # print(col_name,col_plane,col_pos,col_stat)
        elif header_ended:
            if s[0] == '#':
                break
            # print(s)
            bpm_names.append(s[col_name].lower())
            bpm_pos.append(float(s[col_pos])/1.e6)
            bpm_plane.append(s[col_plane])
            bpm_stat.append(s[col_stat])
            # print(s)
    return bpm_names, bpm_plane, bpm_pos, bpm_stat

def get_difference_trajectory(file_orbit,file_trajectory,plane,first_bpm=None,bpm_blacklist=[]):
    orbit_bpm_names, orbit_bpm_plane, orbit_bpm_pos, orbit_bpm_stat = parseYASPfile("YASP_injectionOrbitReference.txt")
    traj_bpm_names, traj_bpm_plane, traj_bpm_pos, traj_bpm_stat = parseYASPfile("YASP_injectionTrajectory.txt")

    if first_bpm:
        before_first_bpm = True
    else:
        before_first_bpm = False
    bpm_names_out = []
    bpm_pos_out = []
    for p, o_n, o_p, o_s, t_n, t_p, t_s in zip(orbit_bpm_plane,orbit_bpm_names, orbit_bpm_pos, orbit_bpm_stat,traj_bpm_names, traj_bpm_pos, traj_bpm_stat):
        if p==plane:
            assert o_n == t_n
            if o_n == first_bpm:
                before_first_bpm = False
            if before_first_bpm==False:
                if o_s == t_s == 'OK':
                    if o_n in bpm_blacklist:
                        continue
                    bpm_names_out.append(o_n)
                    bpm_pos_out.append(t_p-o_p)
    return bpm_names_out, bpm_pos_out


bpm_names, bpm_pos = get_difference_trajectory("YASP_injectionOrbitReference.txt","YASP_injectionTrajectory.txt",'H',"bph.12008",["bph.42008"])


# MAD-X 

madx = Madx()
madx.call("sps/sps.seq")
madx.call("sps/strengths/lhc_q20.str")
madx.call("sps/toolkit/macro.madx")
madx.beam(particle='proton',pc=26)
madx.use(sequence="sps")
madx.input("seqedit,sequence=sps; cycle,start=BPH.12008; flatten; endedit;")
madx.use(sequence="sps")

qx0=20.13
qy0=20.18
madx.input(f"exec, sps_match_tunes({qx0},{qy0});")
madx.input("actcse.31632,volt:=actcse31632_volt,lag:=actcse31632_lag,freq=200;")
madx.input("acl.31735,volt:=acl31735_volt,lag:=acl31735_lag,freq=800;")


# madx.input("select, flag=error, clear; select, flag=error, pattern=qd.11710; ealign, dx=-0.0048;")
# madx.input("select, flag=error, clear; select, flag=error, pattern=qf.11810; ealign, dx=-0.00297;")
# madx.input("select, flag=error, clear; select, flag=error, pattern=qda.11910; ealign, dx=-0.0048;")

'''
madx.input("savebeta,sequence=sps,label=betastart,place=#S;")
tw = madx.twiss()
madx.input("select,flag=mytwiss,clear;")
for n in bpm_names:
    madx.input("select,flag=mytwiss, pattern=%s;"%n)
tw = madx.twiss(table="mytwiss",beta0="betastart")
tw.x[tw.selected_rows()]
'''

line = xt.Line.from_madx_sequence(madx.sequence.sps,deferred_expressions=True)
line.particle_ref = xt.Particles(q0=1, mass0=xt.PROTON_MASS_EV,p0c=26e9)

line.vars['actcse31632_volt'] = 4.5 #
line.vars['actcse31632_lag'] = 0.5
line.vars['acl31735_volt'] = 0.45
line.vars['acl31735_lag'] = 0.5*0

# tw_p = line.twiss(method='4d')
tw_p = line.twiss(method='4d')
indx = [i for i,n in enumerate(tw_p.name) if n in bpm_names]
tw_init = tw_p.get_twiss_init(tw_p.name[0])


''
tw_p = line.twiss(method='4d')
bpm_indx = [i for i,n in enumerate(tw_p.name) if n in bpm_names]
tw_init = tw_p.get_twiss_init(tw_p.name[0])
def trajectory_xtrack(line,bpm_indx,tw_init,x,px,delta,zeta):
    tw_init.x = x
    tw_init.px = px
    tw_init.delta = delta
    tw_init.zeta = zeta
    tw = line.twiss(start=tw_p.name[0],end=tw_p.name[-1],init=tw_init)
    return tw.x[indx], tw.delta[indx]

class ActionTrajectory(xt.Action):
    
    def __init__(self, line, bpm_indx, bpm_names, tw_init, target_trajectory):
        self.line = line
        self.bpm_indx = bpm_indx
        self.bpm_names = bpm_names
        self.tw_init = tw_init
        self.target_trajectory = target_trajectory
        line.vars['x_start'] = 0
        line.vars['px_start'] = 0
        line.vars['delta_start'] = 0
        line.vars['zeta_start'] = 0

    def run(self):
        x = self.line.vars['x_start']._value
        px = self.line.vars['px_start']._value
        delta = self.line.vars['delta_start']._value
        zeta = self.line.vars['zeta_start']._value

        trajectory_x,trajectory_delta = trajectory_xtrack(self.line,self.bpm_indx,self.tw_init,x,px,delta,zeta)
        difference = trajectory_x-self.target_trajectory
        result = {}
        for i in range(len(bpm_indx)):
            result[bpm_names[i]] = trajectory_x[i]
        out_dict = {'std': np.std(difference), 'mean': np.mean(difference), 'trajectory_x': trajectory_x,'trajectory_delta': trajectory_delta}
        out_dict.update(result)
        return out_dict


    
# Build action object
action = ActionTrajectory(line,bpm_indx, bpm_names, tw_init, bpm_pos)
print(action.run())

targets = []
for i in range(len(bpm_indx)):
    targets.append(action.target(bpm_names[i],value=bpm_pos[i],tol=1e-4))


opt = line.match(solve=True,n_steps_max=10,verbose=True,assert_within_tol=False,
        vary=[xt.VaryList(['x_start','px_start'],step=1e-4,tag='t'),
              xt.VaryList(['delta_start'],step=1e-4,tag='ld'),
              xt.VaryList(['zeta_start'],step=1e-1,tag='lz',weight=None, limits=(-0.8, 0.8))],
        targets=targets) #[action.target('std',0,tol=1e-4,tag='std'),action.target('mean',0,tol=1e-4,tag='mean')])
plt.figure();plt.plot(action.run()['trajectory_x']);plt.plot(bpm_pos)

'''

def trajectory_(line,bpm_names):
    tw_p = line.twiss(method='4d')
    indx = [i for i,n in enumerate(tw_p.name) if n in bpm_names]
    tw_init = tw_p.get_twiss_init(tw_p.name[0])
    def twiss(args):
        x=args[0]
        px=args[1]
        delta=args[2]
        zeta=args[3]
        tw_init.x = x
        tw_init.px = px
        tw_init.delta = delta
        tw_init.zeta = zeta
        tw = line.twiss(start=tw_p.name[0],end=tw_p.name[-1],init=tw_init)
        return tw.x[indx]
    return twiss
trajectory = trajectory_(line,bpm_names)

def fnc2optimize_(target_trajectory):
    def difference(args):
        return np.std(trajectory(args)-target_trajectory)
    return difference

fnc2optimize = fnc2optimize_(bpm_pos)

res = minimize(fnc2optimize,x0=[-5.900e-04,-1.025e-05,5.663e-04,-1.860e-01],method='Nelder-Mead',bounds=((-1e-2,1e-2),(-1e-4,1e-4),(-1.5e-3,1.5e-3),(-0.75,0.75)))


plt.figure()
plt.plot(trajectory(res.x))
plt.plot(bpm_pos)
'''



