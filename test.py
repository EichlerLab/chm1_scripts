#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib
pylab.ioff()


print matplotlib.is_interactive()

for i in range(3):
    plt.plot(np.random.rand(10))

print "showing?????"
plt.show(block=True)




