import xtrack as xt
from cpymad.madx import Madx

madx = Madx()
madx.call("sps/sps.seq")
madx.call("sps/strengths/lhc_q20.str")
madx.call("sps/toolkit/macro.madx")
madx.beam(particle='proton',pc=26)
madx.use(sequence='sps')

line = xt.Line.from_madx_sequence(madx.sequence.sps,deferred_expressions=True)
line.particle_ref = xt.Particles(q0=1, mass0=xt.PROTON_MASS_EV,p0c=26e9)

line.to_json('sps.json')