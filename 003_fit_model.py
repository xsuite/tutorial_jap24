import xtrack as xt
import numpy as np
import matplotlib.pyplot as plt

line = xt.Line.from_json('sps.json')
line.cycle('bph.12008')

tw = line.twiss4d()

line.vars.default_to_zero = True
line.set('actcse.31632', voltage='actcse31632_volt', lag='actcse31632_lag', frequency=200e6)
line.set('acl.31735', voltage='acl31735_volt', lag='acl31735_lag', frequency=800e6)
line.vars.default_to_zero = False




# Switch on RF
line['actcse31632_volt'] = 4.5e6
line['actcse31632_lag'] = 180.
line['acl31735_volt'] = 0.45e6
line['acl31735_lag'] = 0.

# Load data for injection trajectory
bpm_data = xt.json.load('bpm_data.json')
bpm_meas = xt.Table(data={
    'name': np.array(bpm_data['bpm_name']),
    'x_meas': np.array(bpm_data['x_position']),
})

line['x_inj'] = 0.
line['px_inj'] = 0.
line['delta_inj'] = 0.
line['zeta_inj'] = 0.

# line.vars.update({
#  'x_inj': -0.0005635430558554506,
#  'px_inj': -9.739671862890309e-06,
#  'delta_inj': 0.0005548018233742175,
#  'zeta_inj': -0.34736514218425374})

tw_inj = line.twiss(
    betx=1, bety=1,
    x=line['x_inj'], px=line['px_inj'],
    delta=line['delta_inj'], zeta=line['zeta_inj'])

class ActionTrajectory(xt.Action):

    def run(self):
        tw = line.twiss(
            betx=1, bety=1,
            x=line['x_inj'], px=line['px_inj'],
            delta=line['delta_inj'], zeta=line['zeta_inj'])
        return tw

action = ActionTrajectory()

tw_inj_bpms = tw_inj.rows[bpm_meas['name']]

targets = []
for nn in bpm_meas['name']:
    targets.append(action.target('x', bpm_meas['x_meas', nn], at=nn))

opt = line.match(
    solve=False,
    vary=[
        xt.Vary('x_inj', step=1e-4),
        xt.Vary('px_inj', step=1e-6),
        xt.Vary('delta_inj', step=1e-5),
        xt.Vary('zeta_inj', step=1e-3),
    ],
    targets=targets
)

opt.step(10)

tw_inj_opt = line.twiss(
    betx=1, bety=1,
    x=line['x_inj'], px=line['px_inj'],
    delta=line['delta_inj'], zeta=line['zeta_inj'])

plt.close('all')
plt.figure(1)
plt.plot(tw_inj_opt.s, tw_inj_opt.x)
plt.plot(tw_inj_bpms.s, tw_inj_opt.rows[bpm_meas['name']].x, 'o')
plt.plot(tw_inj_bpms.s, bpm_meas['x_meas'], 'o')

plt.figure(2)
plt.plot(tw_inj_bpms.s, tw_inj_opt.rows[bpm_meas['name']].x)
plt.plot(tw_inj_bpms.s, bpm_meas['x_meas'])

plt.show()