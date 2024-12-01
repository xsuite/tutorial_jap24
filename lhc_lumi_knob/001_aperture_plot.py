import xtrack as xt
import numpy as np
import matplotlib.pyplot as plt

lhc = xt.Multiline.from_json('lhc_bb.json')

lhc['beambeam_scale'] = 0.

lhc_sequences = '''
    lhcb1: sequence;
    endsequence;
    lhcb2: sequence;
    endsequence;
'''

with open('APERTURE_EYETS 2023-2024.seq', 'r') as f:
    aperture_file = f.read()

input_string = lhc_sequences + aperture_file

aper_env = xt.Environment()
aper_loader = xt.mad_parser.loader.MadxLoader(env=aper_env)
builders = aper_loader.load_string(input_string, build=False)
builder_ap1, builder_ap2 = builders



def aperture_line(builder, table_s_ref):
    # The aperture file expects there to already be some markers in the line,
    # so we add them here.
    # print(f'Building aperture-only version of {line.name}')

    relational_markers = set(p.from_ for p in builder.components)

    for ip in relational_markers:
        builder.new(ip, 'marker', at=table_s_ref['s', ip])

    return builder.build()

lhc.b1.cycle('ip1')
lhc.b2.cycle('ip1')
table_s_b1 = lhc.b1.twiss4d(betx=1, bety=1, x=0, y=0, px=0, py=0).rows['ip.*'].cols['s']
table_s_b2 = lhc.b2.twiss4d(betx=1, bety=1, x=0, y=0, px=0, py=0).rows['ip.*'].cols['s']

circum_b2 = lhc.b2.get_length()
for nn in table_s_b2.name:
    table_s_b2['s', nn] = circum_b2 - table_s_b2['s', nn]

ind_patch = table_s_b2.rows.indices['ip1.l1'][0]
table_s_b2.s[ind_patch] = circum_b2

line_aper1 = aperture_line(builder_ap1, table_s_b1)
line_aper2 = aperture_line(builder_ap2, table_s_b2)

def insert_apertures(line, line_aper, reverse=False):
    print(f'Inserting apertures into {line.name}')
    tt = line_aper.get_table()
    tt_apertures = tt.rows[tt.element_type != 'Marker']
    tt_apertures = tt_apertures.rows[tt_apertures.element_type != 'Drift']

    if reverse:
        tt_apertures.s = line.get_length() - tt_apertures.s

    line._insert_thin_elements_at_s(elements_to_insert=[
        (tt_apertures['s', nn], [(nn, line_aper[nn])]) for nn in tt_apertures.name
        if not nn == '_end_point'
    ])

lhc.discard_trackers()
insert_apertures(lhc.b1, line_aper1)
insert_apertures(lhc.b2, line_aper2, reverse=True)

# lhc.b1.cycle('ip5')
# lhc.b2.cycle('ip5')


arc_vars = [f'ab.a{ab}' for ab in ['12', '23', '34', '45', '56', '67', '78', '81']]
old_vars = {}
for var in arc_vars:
    old_vars[var] = lhc.vars[var]
    lhc.vars[var] = 0

print('Start twiss')
tw1_thin = lhc.b1.twiss4d()
tw2_thin = lhc.b2.twiss4d()
print('Done twiss')

print('Start survey')
sv1 = lhc.b1.survey()
sv2 = lhc.b2.survey()
print('Done survey')


# Compute offsets
# ===============

def offset_elements(line, survey, dir=1):
    tt = line.get_table()
    apertypes = ['LimitEllipse', 'LimitRect', 'LimitRectEllipse', 'LimitRacetrack']
    aper_idx = np.isin(tt.element_type, apertypes)
    for nn in tt.rows[aper_idx].name:
        el = line[nn]
        mech_sep = el.extra['mech_sep'] * dir
        x = survey['X', nn]
        el.shift_x = mech_sep / 2 - x

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
    sigx = 13 * np.sqrt(2.5e-6 / 450 * 0.938 * bx) + abs(dx) * 8e-4

    return s, sx, x, sigx


# Make plots
# ==========

def plot_apertures(line, twiss, survey):
    aper_idx = offset_elements(line, survey, dir=1 if line.name == 'lhcb1' else -1)

    tw_ap = twiss.rows[aper_idx]
    sv_ap = survey.rows[aper_idx]
    ap_extent = np.array([get_aperture_size(line[nn]) for nn in tw_ap.name])
    ap_offset = np.array([line[nn].shift_x for nn in tw_ap.name])

    upper = ap_offset + ap_extent[:, 0] + sv_ap.X
    lower = ap_offset + ap_extent[:, 1] + sv_ap.X

    plt.fill_between(tw_ap.s, upper, lower, alpha=0.5, color="k")
    plt.plot(tw_ap.s, upper, color="k")
    plt.plot(tw_ap.s, lower, color="k")


def plot_beam_size(twiss, survey, color):
    s, sx, x, sigx = compute_beam_size(survey, twiss)
    plt.fill_between(s, x - sigx + sx, x + sigx + sx, alpha=0.5, color=color)

plt.close('all')

plot_apertures(lhc.b1, tw1_thin, sv1)
plot_apertures(lhc.b2, tw2_thin, sv2)

plot_beam_size(tw1_thin, sv1, color='b')
plot_beam_size(tw2_thin, sv2, color='r')

plt.show()