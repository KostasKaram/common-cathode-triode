
import numpy
import math
import matplotlib.pyplot as plt
import json
from triode_amp import *


# read parameters file

# parse parameters from file
params = json.load(open('default_params.json'))

# sample rate
fs = params['os_factor'] * 44100.0

# linear gain applied to full-scale sine
gain = math.pow(10.0, params['db_gain'] / 10.0)

# sine signal
f0 = params['sine_freq']
Vi, time = (util.gen_sine(gain, f0, fs, params['sine_duration']))
Vi = gain*Vi

# tube object
tube = triode_model.triode_model(
    G     = params['G'],
    mu    = params['mu'],
    gamma = params['gamma'],
    C     = params['C'],
    Gg    = params['Gg'],
    xi    = params['xi'],
    Cg    = params['Cg'],
    Ig0   = params['Ig0']
)

# circuit object
triode_stage = triode_preamp.triode_preamp(fs)
triode_stage.setup(
    Rg  = params['Rg'],
    Rk  = params['Rk'],
    Rp  = params['Rp'],
    Ri  = params['Ri'],
    Ci  = params['Ci'],
    Cf  = params['Cf'],
    Ck  = params['Ck'],
    Vpp = params['Vpp']
)
triode_stage.discretize()

# pre-allocate output
Vo  = numpy.zeros_like(Vi)

# process signal
for i in range(0, len(Vi)):

    triode_stage.process(Vi[i], tube)
    Vo[i]  = triode_stage.Vout

# plot voltages
plt.figure(1)
plt.plot(time, Vi,  label='Input')
plt.plot(time, Vo,  label='Output')
plt.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
plt.show()