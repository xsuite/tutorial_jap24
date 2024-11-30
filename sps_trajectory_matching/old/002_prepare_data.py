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

dct = {
    'bpm_name': np.array(bpm_names),
    'x_position': np.array(bpm_pos),
}

tt_data = xt.Table(data=dct, index='bpm_name')

xt.json.dump(dct, 'bpm_data.json')