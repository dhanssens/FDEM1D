#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    EXAMPLE FILE
    ============
    Example file for FDEM1D.

    Use:
    Calculates the forward response and sensitivity distribution for a given
    layered half-space and loop-loop configuration. Typical characteristics
    of the half-space are stored in the Model object (M) while the sensor
    characteristics are stored in the Sensor object (S). This file is
    used as example.

    Cite:
    Hanssens, D., Delefortrie, S., De Pue, J., Van Meirvenne, M., and De Smedt, P., 2018,
    Frequency domain electromagnetic 1D forward and sensitivity modelling: Overview,
    assessment and practical implementation(including MATLAB, Python and pseudo-code):
    Submitted to IEEE Geoscience and Remote Sensing Magazine.

    Code created by Daan Hanssens and Jan De Pue
    Ghent University, Belgium
    April, 2017
"""

# Import dependencies
import numpy as np
import pylab
import FDEM1D

###############################################################################
# ------------------------------------------------------- User-input -------- #
###############################################################################

#
# Sensor object (S)
#

x = 40.0                                                                        # x-coordinate receiver (m)
y = 0.0                                                                        # y-coordinate receiver (m)
z = -0.16                                                                      # z-coordinate receiver (m) - positive z-axis pointed down
height = 0.16                                                                  # Height of transmitter (m)
freq = 90e3                                                                     # Frequency (Hz)
mom = 1.0                                                                      # Transmitter moment (A.m²)
ori = 'ZZ'                                                                     # Coil orientation (-) - case sensitive

S = FDEM1D.Sensor(x, y, z, height, freq, mom, ori)                             # Sensor object


#
# Model object (M)
#

sus = np.linspace(100.0, 100.0, 100) * 1e-5                                     # Susceptibility of layer(s) (-)
con = np.linspace(.001, .001, 100)                                              # Conductivity of layer(s) (S/m)
perm = np.linspace(1.0, 1.0, 100) * 1e-12                                       # Permittivity of layer(s) (F/m)
thick = np.logspace(np.log10(.1), np.log10(.5), 100)                            # Layer(s) thickness (m)

M = FDEM1D.Model(thick, sus, con, perm)                                         # Model object


###############################################################################
# ------------------------------------------------------- Modelling --------- #
###############################################################################

# Define sensitivity parameter
par_sens = 'sus'                                                                   # Variable to be tested in sensitivity analysis

# Reflection coefficient approach (default)
[FWD_IP_RC, FWD_QP_RC] = FDEM1D.Calculate(S, M).forward()                          # Forward
[SENS_IP_RC, SENS_QP_RC, ERROR_RC] = FDEM1D.Calculate(S, M).sensitivity(par_sens)  # Sensitivity

# Propagation Matrix approach
[FWD_IP_PM, FWD_QP_PM] = FDEM1D.Calculate(S, M, method='PM').forward()                          # Forward
[SENS_IP_PM, SENS_QP_PM, ERROR_PM] = FDEM1D.Calculate(S, M, method='PM').sensitivity(par_sens)  # Sensitivity


###############################################################################
# ------------------------------------------------------- Plot -------------- #
###############################################################################

fig = pylab.figure()
ax = fig.add_subplot(211)
ax.plot(M.depth[:-1], SENS_IP_RC[:-1], '.k')
ax.plot(M.depth[:-1], SENS_IP_PM[:-1], '.r')
ax.legend(['RC', 'PM'])
ax.set_xlabel('Depth (m)')
ax.set_ylabel('IP Sensitivity: %s' % par_sens)
ax = fig.add_subplot(212)
ax.plot(M.depth[:-1], SENS_QP_RC[:-1], '.k')
ax.plot(M.depth[:-1], SENS_QP_PM[:-1], '.r')
ax.legend(['RC', 'PM'])
ax.set_xlabel('Depth (m)')
ax.set_ylabel('QP Sensitivity: %s' % par_sens)
pylab.show()

