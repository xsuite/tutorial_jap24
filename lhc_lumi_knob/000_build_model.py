import xtrack as xt
import matplotlib.pyplot as plt

lhc = xt.Multiline.from_json('../../acc-models-lhc/xsuite/lhc.json')

# env = xt.load_madx_lattice('../../acc-models-lhc/lhc.seq',
#                            reverse_lines='lhcb2')
# lhc = xt.Multiline(lines={'b1': env['lhcb1'].copy(), 'b2': env['lhcb2'].copy()})

particle_ref = xt.Particles(energy0=6.8e12, mass0=xt.PROTON_MASS_EV)
lhc.b1.particle_ref = particle_ref
lhc.b2.particle_ref = particle_ref

lhc.b1.cycle('ip7')
lhc.b2.cycle('ip7')

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
# twb1.plot('x y')

opt_b1 = lhc.b1.match_knob(
            run=False,
            knob_name='shift_h_ip1.b1',
            start='e.ds.l1.b1', end='s.ds.r1.b1',
            betx=1, bety=1, x=0, y=0, px=0, py=0, # <- initial conditions
            vary=xt.VaryList(['acbch5.r1b1', 'acbyhs4.r1b1',
                              'acbch6.l1b1', 'acbyh4.l1b1'],
                              step=1e-6),
            targets=[
                xt.TargetSet(x=1e-3, px=0, at='ip1'),
                xt.TargetSet(x=0, px=0, at=xt.END),
            ])
opt_b1.solve()
opt_b1.generate_knob()

opt_b2 = lhc.b2.match_knob(
            run=False,
            knob_name='shift_h_ip1.b2',
            start='e.ds.l1.b2', end='s.ds.r1.b2',
            betx=1, bety=1, x=0, y=0, px=0, py=0, # <- initial conditions
            vary=xt.VaryList(['acbyh4.r1b2', 'acbch6.r1b2',
                              'acbyhs4.l1b2', 'acbch5.l1b2'],
                              step=1e-6),
            targets=[
                xt.TargetSet(x=1e-3, px=0, at='ip1'),
                xt.TargetSet(x=0, px=0, at=xt.END),
            ])
opt_b2.solve()
opt_b2.generate_knob()

lhc['shift_h_ip1.b1'] = 2
lhc['shift_h_ip1.b2'] = -1
tw_b1 = lhc.b1.twiss4d(zero_at='ip1', strengths=True)
tw_b2 = lhc.b2.twiss4d(zero_at='ip1', strengths=True)

plt.close('all')
plt.figure(1)
plt.plot(tw_b1.s, tw_b1.x*1e3, color='b', label='b1')
plt.plot(tw_b2.s, tw_b2.x*1e3, color='r', label='b2')
plt.legend()
plt.xlim(-100, 100)
plt.ylim(-5, 5)
plt.xlabel('s [m]')
plt.ylabel('x [mm]')
plt.grid(True, alpha=0.5)


lhc['shift_h_ip1.b1'] = 2
lhc['shift_h_ip1.b2'] = -1
lhc.install_beambeam_interactions(
    clockwise_line="b1",
    anticlockwise_line="b2",
    ip_names=["ip1", "ip2", "ip5", "ip8"],
    delay_at_ips_slots=[0, 891, 0, 2670],
    num_long_range_encounters_per_side={"ip1": 25, "ip2": 20, "ip5": 25, "ip8": 20},
    num_slices_head_on=11,
    harmonic_number=35640,
    bunch_spacing_buckets=10,
    sigmaz=0.076,
)
lhc.build_trackers()

lhc.b2.twiss_default['reverse'] = False # !!!!!!!!!!!!!!!!!!!!!!!!!!!!

lhc.configure_beambeam_interactions(
    num_particles=1.8e11, nemitt_x=2.e-6, nemitt_y=2.e-6, crab_strong_beam=False
)

fp_sep_on = lhc.b1.get_footprint(
    freeze_longitudinal=True,
    nemitt_x=2.5e-6,
    nemitt_y=2.5e-6,
    r_range=(0.1, 5),
    linear_rescale_on_knobs=[
        xt.LinearRescale(knob_name="beambeam_scale", v0=0.0, dv=0.1)
    ],
)

lhc['shift_h_ip1.b1'] = 0
lhc['shift_h_ip1.b2'] = 0

fp_sep_off = lhc.b1.get_footprint(
    freeze_longitudinal=True,
    nemitt_x=2.5e-6,
    nemitt_y=2.5e-6,
    r_range=(0.1, 5),
    linear_rescale_on_knobs=[
        xt.LinearRescale(knob_name="beambeam_scale", v0=0.0, dv=0.1)
    ],
)

plt.figure()
fp_sep_off.plot(color='C2', alpha=0.5, label='separation off')
fp_sep_on.plot(color='C3', alpha=1, label='separation on')

plt.show()

lhc.to_json('lhc_bb.json')
