

import pickle
import numpy as np
from dolfyn.adv import api as adv
dat = adv.read('S_B_202.vec')
dat.set_inst2head_rotmat(np.zeros([3,3]))
a = adv.correct_motion(dat, accel_filtfreq=0.1)
f=a