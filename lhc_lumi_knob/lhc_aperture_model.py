import xtrack as xt
import numpy as np
import matplotlib.pyplot as plt


class LhcApertureModel:

    def __init__(self, lhc, aperture_file_path,
                 cycle_to='ip5'):
        self.lhc = lhc
        self.aperture_file_path = aperture_file_path

        assert 'b1' in self.lhc.lines
        assert 'b2' in self.lhc.lines

        # Shallow copies
        self.lhcb1_aper = self.lhc.b1.select()
        self.lhcb2_aper = self.lhc.b2.select()
        self.lhcb1_aper.particle_ref = self.lhc.b1.particle_ref
        self.lhcb2_aper.particle_ref = self.lhc.b2.particle_ref

        build_lines_with_apertures(self.lhcb1_aper, self.lhcb2_aper,
                                   self.aperture_file_path)

        self.lhcb1_aper.cycle(cycle_to)
        self.lhcb2_aper.cycle(cycle_to)

        # Go to mid beam reference frame (straight)
        arc_vars = [f'ab.a{ab}' for ab in ['12', '23', '34', '45', '56', '67', '78', '81']]
        old_vars = {}
        for var in arc_vars:
            old_vars[var] = lhc.vars[var]
            lhc.vars[var] = 0

        # Compute surveys
        self.survey_b1 = self.lhcb1_aper.survey()
        self.survey_b2 = self.lhcb2_aper.survey().reverse()

        # Restore
        for var in arc_vars:
            lhc.vars[var] = old_vars[var]

        # Offset apertures
        offset_elements(self.lhcb1_aper, self.survey_b1)
        offset_elements(self.lhcb2_aper, self.survey_b2)

    def plot_horizontal_aperture_and_beam_envelopes(self, zero_at=None):

        if zero_at is not None:
            s_zero = self.survey_b1['s', zero_at]
        else:
            s_zero = 0

        print('Start twiss')
        tw1_aper = self.lhcb1_aper.twiss4d()
        tw2_aper = self.lhcb2_aper.twiss4d(reverse=True)
        print('Done twiss')
        plot_horizontal_apertures(self.lhcb1_aper, tw1_aper, self.survey_b1, s_zero=s_zero)
        plot_horizontal_apertures(self.lhcb2_aper, tw2_aper, self.survey_b2, s_zero=s_zero)

        plot_horizontal_beam_size(tw1_aper, self.survey_b1, color='b', s_zero=s_zero)
        plot_horizontal_beam_size(tw2_aper, self.survey_b2, color='r', s_zero=s_zero)

        plt.xlabel('s [m]')
        plt.ylabel('x [m]')



def build_lines_with_apertures(lhcb1, lhcb2, aper_file_path):
    with open(aper_file_path, 'r') as f:
        aperture_file = f.read()

    lhc_sequences = '''
        lhcb1: sequence;
        endsequence;
        lhcb2: sequence;
        endsequence;
    '''

    input_string = lhc_sequences + aperture_file

    aper_env = xt.Environment()
    aper_loader = xt.mad_parser.loader.MadxLoader(env=aper_env)
    builders = aper_loader.load_string(input_string, build=False)
    builder_ap1, builder_ap2 = builders

    lhcb1.cycle('ip1')
    lhcb2.cycle('ip1')
    table_s_b1 = lhcb1.get_table().rows['ip.*'].cols['s']
    table_s_b2 = lhcb2.get_table().rows['ip.*'].cols['s']

    circum_b2 = lhcb2.get_length()
    table_s_b2['s'] = circum_b2 - table_s_b2['s']

    # PATCH !!!!!!!!!!!!!!!!!!!
    ind_patch = table_s_b2.rows.indices['ip1.l1'][0]
    table_s_b2.s[ind_patch] = circum_b2
    ind_patch = table_s_b2.rows.indices['ip1'][0]
    table_s_b2.s[ind_patch] = 0

    line_aper1 = aperture_line(builder_ap1, table_s_b1)
    line_aper2 = aperture_line(builder_ap2, table_s_b2)

    lhcb1.build_tracker() # To resolve parents
    lhcb2.build_tracker()
    lhcb1.discard_tracker()
    lhcb2.discard_tracker()
    insert_apertures(lhcb1, line_aper1)
    insert_apertures(lhcb2, line_aper2, reverse=True)

def aperture_line(builder, table_s_ref):
    # The aperture file expects there to already be some markers in the line,
    # so we add them here.
    # print(f'Building aperture-only version of {line.name}')

    relational_markers = set(p.from_ for p in builder.components)

    for ip in relational_markers:
        builder.new(ip, 'marker', at=table_s_ref['s', ip])

    return builder.build()

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

def offset_elements(line, survey):
    tt = line.get_table()
    apertypes = ['LimitEllipse', 'LimitRect', 'LimitRectEllipse', 'LimitRacetrack']
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

def get_horizontal_aperture_size(el):
    if hasattr(el, 'min_x'):
        return el.min_x, el.max_x
    if hasattr(el, 'max_x'):
        return -el.max_x, el.max_x
    return -el.a, el.a


def compute_horizontal_beam_size(survey, twiss, nemitt_x=2.5e-6, sigma_delta=8e-4):
    sx = survey.X
    s = survey.s
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

def plot_horizontal_apertures(line, twiss, survey, s_zero=0):
    tt = line.get_table()
    apertypes = ['LimitEllipse', 'LimitRect', 'LimitRectEllipse', 'LimitRacetrack']
    aper_idx = np.where([tt['element_type', nn] in apertypes for nn in survey.name])[0]

    tw_ap = twiss.rows[aper_idx]
    sv_ap = survey.rows[aper_idx]
    ap_extent = np.array([get_horizontal_aperture_size(line[nn]) for nn in tw_ap.name])
    ap_offset = np.array([line[nn].shift_x for nn in tw_ap.name])

    upper = ap_offset + ap_extent[:, 0] + sv_ap.X
    lower = ap_offset + ap_extent[:, 1] + sv_ap.X

    plt.fill_between(tw_ap.s - s_zero, upper, lower, alpha=1., color='lightgrey')
    plt.plot(sv_ap.s - s_zero, upper, color="k")
    plt.plot(sv_ap.s - s_zero, lower, color="k")

def plot_horizontal_beam_size(twiss, survey, color, s_zero=0):
    s, sx, x, sigx = compute_horizontal_beam_size(survey, twiss)
    plt.fill_between(s - s_zero, x - sigx + sx, x + sigx + sx, alpha=0.5, color=color)
