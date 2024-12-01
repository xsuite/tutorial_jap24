import xtrack as xt
import numpy as np
import matplotlib.pyplot as plt

from build_lines_with_apertures import ApertureModel

lhc = xt.Multiline.from_json('lhc_bb.json')

lhc['shift_h_ip1.b1'] = 2
lhc['shift_h_ip1.b2'] = -1


lhc['beambeam_scale'] = 0.

aper_model = ApertureModel(lhc=lhc,
                           aperture_file_path='APERTURE_EYETS 2023-2024.seq')

plt.close('all')
plt.figure(figsize=(6.4, 6.4))
aper_model.plot_horizontal_aperture_and_beam_envelopes(zero_at='ip1')
plt.xlim(-300, 300)
plt.ylim(-0.20, 0.20)
plt.show()