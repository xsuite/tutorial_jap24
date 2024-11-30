import xtrack as xt

env = xt.load_madx_lattice('../../acc-models-lhc/lhc.seq',
                           reverse_lines='lhcb2')
lhc = xt.Multiline(lines={'b1': env['lhcb1'].copy(), 'b2': env['lhcb2'].copy()})

particle_ref = xt.Particles(energy0=6.8e12, mass0=xt.PROTON_MASS_EV)
lhc.b1.particle_ref = particle_ref
lhc.b2.particle_ref = particle_ref

lhc.b1.cycle('ip7')
lhc.b2.cycle('ip7_reversed') # !!!!!

lhc.b2.twiss_default['reverse'] = True

# Load optics
lhc.vars.load_madx(
    "../../acc-models-lhc/strengths/ATS_Nominal/2024/ats_30cm.madx")

lhc['on_disp'] = 1.
lhc['on_xx1_v'] = 160.0
lhc['on_x1_v'] = 160.0
lhc['on_sep1_h'] = 0.
lhc['on_sep2h'] = -1.0
lhc['on_x2v'] = 200.0
lhc['on_xx5_h'] = 160.0
lhc['on_x5_h'] = 160.0
lhc['on_sep5_v'] = 0.
lhc['on_sep8h'] = -1.0
lhc['on_x8v'] = 200.0

twb1 = lhc.b1.twiss4d()
twb1.plot('x y')

# opt = lhc.b1.match_knob(
#             run=False,
#             knob_name='shift_h_ip5.b1',
#             ele_start='e.ds.l1.b1',
#             ele_stop='s.ds.r1.b1',
#             betx=1, bety=1, x=0, y=0, px=0, py=0, # <- initial conditions
#             vary=[
#                 xt.VaryList(['acbyvs4.l5b1', 'acbcv6.r5b1', 'acbyv4.r5b1', 'acbcv5.l5b1'],
#                             step=1e-6),
#             ],
#             targets=[
#                 xt.TargetSet(x=1e-3, y=0, px=0, py=0, at='ip5'),
#                 xt.TargetSet(x=0, px=0, y=0, py=0, at=xt.END),
#             ])

