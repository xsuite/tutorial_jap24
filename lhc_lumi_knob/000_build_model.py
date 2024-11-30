import xtrack as xt

env = xt.load_madx_lattice('../../acc-models-lhc/lhc.seq',
                           reverse_lines='lhcb2')
lhc = xt.Multiline(lines={'b1': env['lhcb1'].copy(), 'b2': env['lhcb2'].copy()})

particle_ref = xt.Particles(energy0=6.8e12, mass0=xt.PROTON_MASS_EV)
lhc.b1.particle_ref = particle_ref
lhc.b2.particle_ref = particle_ref

lhc.b2.twiss_default['reverse'] = True

# Load optics
lhc.vars.load_madx(
    "../../acc-models-lhc/strengths/ATS_Nominal/2024/ats_30cm.madx")

# Wipe all orbit knobs
for kk in lhc.vars.get_table().rows['on.*'].name:
    lhc.vars[kk] = 0

lhc['on_disp'] = 1.
lhc['on_x1_v'] = 160.
lhc['on_xx1_v'] = 160.

twb1 = lhc.b1.twiss4d()
twb1.plot('y')
