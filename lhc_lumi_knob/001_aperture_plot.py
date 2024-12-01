import xtrack as xt
import numpy as np
import matplotlib.pyplot as plt

from build_lines_with_apertures import build_lines_with_apertures

lhc = xt.Multiline.from_json('lhc_bb.json')

lhc['shift_h_ip1.b1'] = 2
lhc['shift_h_ip1.b2'] = -1


lhc['beambeam_scale'] = 0.

build_lines_with_apertures(lhc, 'APERTURE_EYETS 2023-2024.seq')

lhc.b1.cycle('ip5')
lhc.b2.cycle('ip5')

arc_vars = [f'ab.a{ab}' for ab in ['12', '23', '34', '45', '56', '67', '78', '81']]
old_vars = {}
for var in arc_vars:
    old_vars[var] = lhc.vars[var]
    lhc.vars[var] = 0

print('Start twiss')
tw1_thin = lhc.b1.twiss4d()
tw2_thin = lhc.b2.twiss4d(reverse=True)
print('Done twiss')

print('Start survey')
sv1 = lhc.b1.survey()
sv2 = lhc.b2.survey().reverse()
print('Done survey')




# Compute offsets
# ===============

def offset_elements(line, survey):
    tt = line.get_table()
    apertypes = ['LimitEllipse', 'LimitRect', 'LimitRectEllipse', 'LimitRacetrack']
    # aper_idx = np.isin(tt.element_type, apertypes)
    aper_idx = np.where([tt['element_type', nn] in apertypes for nn in survey.name])[0]
    mech_sep_arr = survey.s * 0 + np.nan
    for ii in aper_idx:
        nn = survey.name[ii]
        el = line[nn]
        mech_sep = el.extra['mech_sep'] # * dir
        x = survey['X', nn]
        el.shift_x = mech_sep / 2 - x

        mech_sep_arr[ii] = mech_sep
    survey['mech_sep'] = mech_sep_arr

    return aper_idx

# Convenience function to compute aperture size and beam sizes
# ============================================================

def get_aperture_size(el):
    if hasattr(el, 'min_x'):
        return el.min_x, el.max_x
    if hasattr(el, 'max_x'):
        return -el.max_x, el.max_x
    return -el.a, el.a


def compute_beam_size(survey, twiss):
    sx = survey.X
    s = twiss.s
    x = twiss.x
    bx = twiss.betx
    dx = twiss.dx
    nemitt_x = 2.5e-6
    gamma0 = twiss.gamma0
    n_sigmas = 13.
    sigma_delta = 8e-4
    # sigx = 13 * np.sqrt(2.5e-6 / 450 * 0.938 * bx) + abs(dx) * 8e-4
    sigx = n_sigmas * np.sqrt(nemitt_x / gamma0 * bx) + abs(dx) * sigma_delta

    return s, sx, x, sigx


# Make plots
# ==========

def plot_apertures(line, twiss, survey):
    tt = line.get_table()
    apertypes = ['LimitEllipse', 'LimitRect', 'LimitRectEllipse', 'LimitRacetrack']
    # aper_idx = np.isin(tt.element_type, apertypes)
    aper_idx = np.where([tt['element_type', nn] in apertypes for nn in survey.name])[0]

    tw_ap = twiss.rows[aper_idx]
    sv_ap = survey.rows[aper_idx]
    ap_extent = np.array([get_aperture_size(line[nn]) for nn in tw_ap.name])
    ap_offset = np.array([line[nn].shift_x for nn in tw_ap.name])

    upper = ap_offset + ap_extent[:, 0] + sv_ap.X
    lower = ap_offset + ap_extent[:, 1] + sv_ap.X

    plt.fill_between(tw_ap.s, upper, lower, alpha=1., color='lightgrey')
    plt.plot(tw_ap.s, upper, color="k")
    plt.plot(tw_ap.s, lower, color="k")


def plot_beam_size(twiss, survey, color):
    s, sx, x, sigx = compute_beam_size(survey, twiss)
    plt.fill_between(s, x - sigx + sx, x + sigx + sx, alpha=0.5, color=color)

lhc.build_trackers()

plt.close('all')

offset_elements(lhc.b1, sv1)
offset_elements(lhc.b2, sv2)

plot_apertures(lhc.b1, tw1_thin, sv1)
plot_apertures(lhc.b2, tw2_thin, sv2)

plot_beam_size(tw1_thin, sv1, color='b')
plot_beam_size(tw2_thin, sv2, color='r')

plt.show()