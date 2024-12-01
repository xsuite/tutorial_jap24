import xtrack as xt
import numpy as np


def build_lines_with_apertures(lhc, aper_file_path):
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

    lhc.b1.cycle('ip1')
    lhc.b2.cycle('ip1')
    table_s_b1 = lhc.b1.get_table().rows['ip.*'].cols['s']
    table_s_b2 = lhc.b2.get_table().rows['ip.*'].cols['s']

    circum_b2 = lhc.b2.get_length()
    table_s_b2['s'] = circum_b2 - table_s_b2['s']

    # PATCH !!!!!!!!!!!!!!!!!!!
    ind_patch = table_s_b2.rows.indices['ip1.l1'][0]
    table_s_b2.s[ind_patch] = circum_b2
    ind_patch = table_s_b2.rows.indices['ip1'][0]
    table_s_b2.s[ind_patch] = 0

    line_aper1 = aperture_line(builder_ap1, table_s_b1)
    line_aper2 = aperture_line(builder_ap2, table_s_b2)

    lhc.build_trackers() # To resolve parents
    lhc.discard_trackers()
    insert_apertures(lhc.b1, line_aper1)
    insert_apertures(lhc.b2, line_aper2, reverse=True)

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
