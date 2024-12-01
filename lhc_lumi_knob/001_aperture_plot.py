import xtrack as xt
import numpy as np
import matplotlib.pyplot as plt

from build_lines_with_apertures import ApertureModel

lhc = xt.Multiline.from_json('lhc_bb.json')

lhc['shift_h_ip1.b1'] = 2
lhc['shift_h_ip1.b2'] = -1


lhc['beambeam_scale'] = 0.

aper_model = ApertureModel(lhc=lhc, aperture_file_path='APERTURE_EYETS 2023-2024.seq')

print('Start twiss')
tw1_aper = aper_model.lhcb1_aper.twiss4d()
tw2_aper = aper_model.lhcb2_aper.twiss4d(reverse=True)
print('Done twiss')






# Compute offsets
# ===============



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

plot_apertures(aper_model.lhcb1_aper, tw1_aper, aper_model.survey_b1)
plot_apertures(aper_model.lhcb2_aper, tw2_aper, aper_model.survey_b2)

plot_beam_size(tw1_aper, aper_model.survey_b1, color='b')
plot_beam_size(tw2_aper, aper_model.survey_b2, color='r')

plt.show()