#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    FDEM1D
    ======
    Frequency domain electromagnetic (FDEM) 1D forward and sensitivity modelling

    Cite:
    Hanssens, D., Delefortrie, S., De Pue, J., Van Meirvenne, M., and De Smedt, P., 2018,
    Practical aspects of frequency domain electromagnetic forward and sensitivity modelling 
    of a magnetic dipole in a multi-layered half-space:
    Submitted to IEEE Geoscience and Remote Sensing Magazine.

    :AUTHOR: Daan Hanssens, Jan De Pue
    :CONTACT: daan.hanssens@ugent.be

    :REQUIRES: numpy, scipy, copy
"""

# Import
import numpy as np
import copy
from scipy.constants import mu_0, epsilon_0


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                 SENSOR CHARACTERISTICS                                               #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
class Sensor:
    def __init__(self, x, y, z, height, freq, mom, ori):
        """
            Sensor characteristics

            Parameters
            ----------
            x: float
                X-coordinate receiver (m)

            y: float
                Y-coordinate receiver (m)

            z: float
                Z-coordinate receiver (m)

            height: float
                Height of transmitter (m)

            freq: float
                Operating frequency (Hz)

            mom: float
                Transmitter moment (Am**2)

            ori: string
                Coil orientation (Transmitter('X','Y','Z'), Receiver('X','Y','Z'))
        """

        self.xyz = np.array([x, y, z])            # x,y,z-coordinate receiver (m)
        self.height = height                      # Height of transmitter (m)
        self.freq = freq                          # Frequency (Hz)
        self.mom = mom                            # Transmitter moment (A.m**2)
        self.ori = ori                            # Coil orientation (T(X,Y,Z), R(X,Y,Z))

    @property
    def r(self):
        """
            Calculate coil spacing (m)

            Returns
            -------
            coil_spacing: float
                Coil spacing (m)
        """

        return np.sqrt(self.xyz[0] ** 2 + self.xyz[1] ** 2 +
                       (self.xyz[2] + self.height) ** 2)                     # Coil spacing (m)

    @property
    def omega(self):
        """
            Calculate angular frequency (Rad/s)

            Returns
            -------
            angular_frequency: float
                Angular frequency (Rad/s)
        """

        return np.pi * 2 * self.freq                                         # Angular frequency (Rad/s)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                  MODEL CHARACTERISTICS                                               #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
class Model:
    def __init__(self, thick, sus, con, perm):
        """
            Model characteristics (1D)

            Parameters
            ----------
            thick: np.array
                Layer(s) thickness (m)

            sus: np.array
                Magnetic susceptibility of layer(s) (-)

            con: np.array
                Electrical conductivity of layer(s) (F/m)

            perm: np.array
                Dielectric permittivity of layer(s) (S/m)
        """

        self.thick = np.array(thick)                              # Layer(s) thickness (m)
        self.sus = np.array(sus)                                  # Susceptibility of layer(s) (-)
        self.perm = np.array(perm)                                # Permittivity of layer(s) (F/m)
        self.con = np.array(con)                                  # Conductivity of layer(s) (S/m)

    @property
    def depth(self):
        """
            Calculate depth (m) of model

            Returns
            -------
            depth: np.array
                Depth (m)
        """

        return np.cumsum(self.thick)                                         # Depth (m)

    @property
    def nr_layers(self):
        """
            Calculate number of layers

            Returns
            -------
            nr_layers: int
                Number of layers
        """

        return self.depth.size                                              # Number of layer(s) (-)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                       1D  MODELLING                                                  #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
class Calculate:
    def __init__(self, Sensor, Model, method='RC'):
        """
            Set FDEM (1D) characteristics and methods

            Parameters
            ----------
            Sensor: object
                FDEM.Sensor object

            Model: object
                FDEM.Model object

            method: str
                Reflection coefficient ('RC'; default) or Propagation Matrix ('PM')
        """

        self.Sensor = Sensor
        self.Model = Model
        self.Method = method
        if method not in ['RC', 'PM']:
            raise ValueError(
                "Choose an appropriate method: 'PM' for Propagation Matrix or 'RC' for Reflection Coefficient")
        self.OriginalModel = copy.deepcopy(self.Model)

    def original_model(self, par):
        """
            Store original profile with sensitivity parameter

            Parameters
            ----------
            par: str
                Sensitivity of physical property ('con','sus','perm')

            Returns
            -------
            original_model
                Model with set sensitivity parameter
        """

        return getattr(self.OriginalModel, par)

    def reset_model(self):
        """
            Reset model (for sensitivity calculation)
        """

        self.Model = copy.deepcopy(self.OriginalModel)

    def forward(self):
        """
            Calculate the forward response (ppm) of a given layered half-space and loop-loop configuration

            Returns
            -------
            FWD_IP: float
                IP response (ppm)

            FWD_QP: float
                QP response (ppm)
        """

        # Get magnetic fields
        [H, Hn] = self.magnetic_fields()

        # Normalization of magnetic field (ppm)
        H_rsp = 1e6 * H / Hn

        # Get forward response (ppm)
        FWD_IP = np.real(H_rsp)
        FWD_QP = np.imag(H_rsp)

        return FWD_IP, FWD_QP

    def sensitivity(self, par, pertFactor=1e-2):
        """
            Calculate the sensitivity distribution of a given layered half-space and loop-loop configuration

            Parameters
            ----------
            par: str
                Sensitivity of physical property ('con','sus','perm')

            pertFactor: float
                Perturbation factor (default: 1e-2)

            Returns
            -------
            SENS_IP: np.array
                IP sensitivity distribution

            SENS_QP: np.array
                QP sensitivity distribution

            error: float
                Estimated error on sensitivity
        """

        # Initialize
        nlay = self.original_model(par).size
        FWD_IP_p = np.zeros(nlay)
        FWD_QP_p = np.zeros(nlay)
        FWD_IP_n = np.zeros(nlay)
        FWD_QP_n = np.zeros(nlay)

        # Get original response
        [FWD_IP_ori, FWD_QP_ori] = self.forward()

        # Loop over layer(s)
        for i in range(nlay):

            # Copy original model
            M_ori = copy.deepcopy(self.original_model(par))

            # Get altered response (forward)
            pert = M_ori[i] * pertFactor
            M_ori[i] = M_ori[i] + pert
            setattr(self.Model, par, M_ori)  # Set Perturbed profile
            [FWD_IP_p[i], FWD_QP_p[i]] = self.forward()

            # Get altered response (backward)
            M_ori = copy.deepcopy(self.original_model(par))
            M_ori[i] = M_ori[i] - pert
            setattr(self.Model, par, M_ori)  # Set Perturbed profile
            [FWD_IP_n[i], FWD_QP_n[i]] = self.forward()

            # Reset Model
            self.reset_model()

        # First derivative
        SENS_IP = (FWD_IP_p - FWD_IP_ori) / pert
        SENS_QP = (FWD_QP_p - FWD_QP_ori) / pert

        # Second derivative
        SENS_IP_sd = (FWD_IP_p - 2 * FWD_IP_ori + FWD_IP_n) / pert ** 2
        SENS_QP_sd = (FWD_QP_p - 2 * FWD_QP_ori + FWD_QP_n) / pert ** 2

        # Estimate maximum error
        error = [np.max(SENS_IP_sd) * pert / 2, np.max(SENS_QP_sd) * pert / 2]

        return SENS_IP, SENS_QP, error

    def magnetic_fields(self):
        """
            Calculate magnetic fields for an x-directed transmitter and different X,Y,Z-receiver orientations
            ('XX','XY','XZ'), Y-directed transmitter and different X,Y,Z-receiver orientations ('YX','YY','YZ') and
            Z-directed transmitter and different x,y,z-receiver orientations ('ZX','ZY','ZZ').

            Returns
            -------
            H
                Magnetic field (H/m)

            Hn
                Magnetic field used for normalization (H/m)
        """

        if self.Sensor.ori == 'ZZ':

            # Calculate Hzz (HCP)
            H = self.Sensor.mom / (4 * np.pi) * self.digital_filter(self.rTEp2, 0)
            Hn = self.Sensor.mom / (4 * np.pi) * self.digital_filter(self.rTE02, 0)

        elif self.Sensor.ori == 'ZY':

            # Calculate Hzy (NULL)
            H = -self.Sensor.mom / (4 * np.pi) * self.Sensor.xyz[1] / \
                self.Sensor.r * self.digital_filter(self.rTEn2, 1)
            Hn = self.Sensor.mom / (4 * np.pi) * self.digital_filter(self.rTE02, 0)

        elif self.Sensor.ori == 'ZX':

            # Calculate Hzx (PRP)
            H = -self.Sensor.mom / (4 * np.pi) * self.Sensor.xyz[0] / \
                self.Sensor.r * self.digital_filter(self.rTEn2, 1)
            Hn = self.Sensor.mom / (4 * np.pi) * self.digital_filter(self.rTE02, 0)

        elif self.Sensor.ori == 'XX':

            # Calculate Hxx (VCA)
            H = -self.Sensor.mom / (4 * np.pi) * (
            1 / self.Sensor.r - 2 * self.Sensor.xyz[0] ** 2 / self.Sensor.r ** 3) * self.digital_filter(self.rTEn1, 1) -\
                self.Sensor.mom / (4 * np.pi) * self.Sensor.xyz[0] ** 2 / self.Sensor.r ** 2 * self.digital_filter(
                    self.rTEn2, 0)
            Hn = -self.Sensor.mom / (4 * np.pi) * (
            1 / self.Sensor.r - 2 * self.Sensor.xyz[0] ** 2 / self.Sensor.r ** 3) * self.digital_filter(self.rTE01, 1) -\
                 self.Sensor.mom / (4 * np.pi) * self.Sensor.xyz[0] ** 2 / self.Sensor.r ** 2 * self.digital_filter(
                     self.rTE02, 0)

        elif self.Sensor.ori == 'XY':

            # Calculate Hxy (NULL)
            H = self.Sensor.mom / (2 * np.pi) * self.Sensor.xyz[0] * self.Sensor.xyz[
                1] / self.Sensor.r ** 3 * self.digital_filter(self.rTEn1, 1) - \
                self.Sensor.mom / (4 * np.pi) * self.Sensor.xyz[0] * self.Sensor.xyz[
                    1] / self.Sensor.r ** 2 * self.digital_filter(self.rTEn2, 0)
            Hn = -self.Sensor.mom / (4 * np.pi) * (
            1 / self.Sensor.r - 2 * self.Sensor.xyz[0] ** 2 / self.Sensor.r ** 3) * self.digital_filter(self.rTE01, 1) -\
                 self.Sensor.mom / (4 * np.pi) * self.Sensor.xyz[0] ** 2 / self.Sensor.r ** 2 * self.digital_filter(
                     self.rTE02, 0)

        elif self.Sensor.ori == 'XZ':

            # Calculate Hxz (PRP)
            H = self.Sensor.mom / (4 * np.pi) * self.Sensor.xyz[0] / self.Sensor.r * self.digital_filter(self.rTEp2,
                                                                                                         1)
            Hn = -self.Sensor.mom / (4 * np.pi) * (
            1 / self.Sensor.r - 2 * self.Sensor.xyz[0] ** 2 / self.Sensor.r ** 3) * self.digital_filter(self.rTE01, 1) -\
                 self.Sensor.mom / (4 * np.pi) * self.Sensor.xyz[0] ** 2 / self.Sensor.r ** 2 * self.digital_filter(
                     self.rTE02, 0)

        elif self.Sensor.ori == 'YX':

            # Calculate Hyx (NULL)
            H = self.Sensor.mom / (2 * np.pi) * self.Sensor.xyz[0] * self.Sensor.xyz[
                1] / self.Sensor.r ** 3 * self.digital_filter(self.rTEn1, 1) - \
                self.Sensor.mom / (4 * np.pi) * self.Sensor.xyz[0] * self.Sensor.xyz[
                    1] / self.Sensor.r ** 2 * self.digital_filter(self.rTEn2, 0)
            Hn = -self.Sensor.mom / (4 * np.pi) * (
            1 / self.Sensor.r - 2 * self.Sensor.xyz[1] ** 2 / self.Sensor.r ** 3) * self.digital_filter(self.rTE01, 1) -\
                 self.Sensor.mom / (4 * np.pi) * self.Sensor.xyz[1] ** 2 / self.Sensor.r ** 2 * self.digital_filter(
                     self.rTE02, 0)

        elif self.Sensor.ori == 'YY':

            # Calculate Hyy (VCP)
            H = -self.Sensor.mom / (4 * np.pi) * (
            1 / self.Sensor.r - 2 * self.Sensor.xyz[1] ** 2 / self.Sensor.r ** 3) * self.digital_filter(self.rTEn1, 1) -\
                self.Sensor.mom / (4 * np.pi) * self.Sensor.xyz[1] ** 2 / self.Sensor.r ** 2 * self.digital_filter(
                    self.rTEn2, 0)
            Hn = -self.Sensor.mom / (4 * np.pi) * (
            1 / self.Sensor.r - 2 * self.Sensor.xyz[1] ** 2 / self.Sensor.r ** 3) * self.digital_filter(self.rTE01, 1) -\
                 self.Sensor.mom / (4 * np.pi) * self.Sensor.xyz[1] ** 2 / self.Sensor.r ** 2 * self.digital_filter(
                     self.rTE02, 0)

        elif self.Sensor.ori == 'YZ':

            # Calculate Hyz (NULL)
            H = self.Sensor.mom / (4 * np.pi) * self.Sensor.xyz[1] / self.Sensor.r * self.digital_filter(self.rTEp2,
                                                                                                         1)
            Hn = -self.Sensor.mom / (4 * np.pi) * (
            1 / self.Sensor.r - 2 * self.Sensor.xyz[1] ** 2 / self.Sensor.r ** 3) * self.digital_filter(self.rTE01, 1) -\
                 self.Sensor.mom / (4 * np.pi) * self.Sensor.xyz[1] ** 2 / self.Sensor.r ** 2 * self.digital_filter(
                     self.rTE02, 0)

        else:

            raise ValueError('Transmitter-receiver orientation unknown: should be a case sensitive X, Y, Z combination')

        return H, Hn

    def digital_filter(self, funname, order):
        """
            Solves the Hankel Transform of the zeroth or first order by using a Guptasarma and Singh filtering routine.

            Parameters
            ----------
            funname
                Name of function

            order
                Order of filter

            Returns
            -------
            y
                Solved Hankel Transform
        """

        # Load Guptasarma and Singh filter
        if order == 0:

            # Load 120-point filter
            filter_a = -8.3885
            filter_s = 0.090422646867
            filter_w = np.array(
                [9.62801364263000e-07, -5.02069203805000e-06, 1.25268783953000e-05, -1.99324417376000e-05,
                 2.29149033546000e-05, -2.04737583809000e-05, 1.49952002937000e-05, -9.37502840980000e-06,
                 5.20156955323000e-06, -2.62939890538000e-06, 1.26550848081000e-06, -5.73156151923000e-07,
                 2.76281274155000e-07, -1.09963734387000e-07, 7.38038330280000e-08, -9.31614600001000e-09,
                 3.87247135578000e-08, 2.10303178461000e-08, 4.10556513877000e-08, 4.13077946246000e-08,
                 5.68828741789000e-08, 6.59543638130000e-08, 8.40811858728000e-08, 1.01532550003000e-07,
                 1.26437360082000e-07, 1.54733678097000e-07, 1.91218582499000e-07, 2.35008851918000e-07,
                 2.89750329490000e-07, 3.56550504341000e-07, 4.39299297826000e-07, 5.40794544880000e-07,
                 6.66136379541000e-07, 8.20175040653000e-07, 1.01015545059000e-06, 1.24384500153000e-06,
                 1.53187399787000e-06, 1.88633707689000e-06, 2.32307100992000e-06, 2.86067883258000e-06,
                 3.52293208580000e-06, 4.33827546442000e-06, 5.34253613351000e-06, 6.57906223200000e-06,
                 8.10198829111000e-06, 9.97723263578000e-06, 1.22867312381000e-05, 1.51305855976000e-05,
                 1.86329431672000e-05, 2.29456891669000e-05, 2.82570465155000e-05, 3.47973610445000e-05,
                 4.28521099371000e-05, 5.27705217882000e-05, 6.49856943660000e-05, 8.00269662180000e-05,
                 9.85515408752000e-05, 0.000121361571831000, 0.000149454562334000, 0.000184045784500000,
                 0.000226649641428000, 0.000279106748890000, 0.000343716968725000, 0.000423267056591000,
                 0.000521251001943000, 0.000641886194381000, 0.000790483105615000, 0.000973420647376000,
                 0.00119877439042000, 0.00147618560844000, 0.00181794224454000, 0.00223860214971000,
                 0.00275687537633000, 0.00339471308297000, 0.00418062141752000, 0.00514762977308000,
                 0.00633918155348000, 0.00780480111772000, 0.00961064602702000, 0.0118304971234000, 0.0145647517743000,
                 0.0179219149417000, 0.0220527911163000, 0.0271124775541000, 0.0333214363101000, 0.0408864842127000,
                 0.0501074356716000, 0.0612084049407000, 0.0745146949048000, 0.0900780900611000, 0.107940155413000,
                 0.127267746478000, 0.146676027814000, 0.162254276550000, 0.168045766353000, 0.152383204788000,
                 0.101214136498000, -0.00244389126667000, -0.154078468398000, -0.303214415655000, -0.297674373379000,
                 0.00793541259524000, 0.426273267393000, 0.100032384844000, -0.494117404043000, 0.392604878741000,
                 -0.190111691178000, 0.0743654896362000, -0.0278508428343000, 0.0109992061155000, -0.00469798719697000,
                 0.00212587632706000, -0.000981986734159000, 0.000444992546836000, -0.000189983519162000,
                 7.31024164292000e-05, -2.40057837293000e-05, 6.23096824846000e-06, -1.12363896552000e-06,
                 1.04470606055000e-07])

        elif order == 1:

            # Load 140-point filter
            filter_a = -7.91001919
            filter_s = 0.087967143957
            filter_w = np.array(
                [-6.76671159511000e-14, 3.39808396836000e-13, -7.43411889153000e-13, 8.93613024469000e-13,
                 -5.47341591896000e-13, -5.84920181906000e-14, 5.20780672883000e-13, -6.92656254606000e-13,
                 6.88908045074000e-13, -6.39910528298000e-13, 5.82098912530000e-13, -4.84912700478000e-13,
                 3.54684337858000e-13, -2.10855291368000e-13, 1.00452749275000e-13, 5.58449957721000e-15,
                 -5.67206735175000e-14, 1.09107856853000e-13, -6.04067500756000e-14, 8.84512134731000e-14,
                 2.22321981827000e-14, 8.38072239207000e-14, 1.23647835900000e-13, 1.44351787234000e-13,
                 2.94276480713000e-13, 3.39965995918000e-13, 6.17024672340000e-13, 8.25310217692000e-13,
                 1.32560792613000e-12, 1.90949961267000e-12, 2.93458179767000e-12, 4.33454210095000e-12,
                 6.55863288798000e-12, 9.78324910827000e-12, 1.47126365223000e-11, 2.20240108708000e-11,
                 3.30577485691000e-11, 4.95377381480000e-11, 7.43047574433000e-11, 1.11400535181000e-10,
                 1.67052734516000e-10, 2.50470107577000e-10, 3.75597211630000e-10, 5.63165204681000e-10,
                 8.44458166896000e-10, 1.26621795331000e-09, 1.89866561359000e-09, 2.84693620927000e-09,
                 4.26886170263000e-09, 6.40104325574000e-09, 9.59798498616000e-09, 1.43918931885000e-08,
                 2.15798696769000e-08, 3.23584600810000e-08, 4.85195105813000e-08, 7.27538583183000e-08,
                 1.09090191748000e-07, 1.63577866557000e-07, 2.45275193920000e-07, 3.67784458730000e-07,
                 5.51470341585000e-07, 8.26916206192000e-07, 1.23991037294000e-06, 1.85921554669000e-06,
                 2.78777669034000e-06, 4.18019870272000e-06, 6.26794044911000e-06, 9.39858833064000e-06,
                 1.40925408889000e-05, 2.11312291505000e-05, 3.16846342900000e-05, 4.75093313246000e-05,
                 7.12354794719000e-05, 0.000106810848460000, 0.000160146590551000, 0.000240110903628000,
                 0.000359981158972000, 0.000539658308918000, 0.000808925141201000, 0.00121234066243000,
                 0.00181650387595000, 0.00272068483151000, 0.00407274689463000, 0.00609135552241000,
                 0.00909940027636000, 0.0135660714813000, 0.0201692550906000, 0.0298534800308000, 0.0439060697220000,
                 0.0639211368217000, 0.0916763946228000, 0.128368795114000, 0.173241920046000, 0.219830379079000,
                 0.251193131178000, 0.232380049895000, 0.117121080205000, -0.117252913088000, -0.352148528535000,
                 -0.271162871370000, 0.291134747110000, 0.317192840623000, -0.493075681595000, 0.311223091821000,
                 -0.136044122543000, 0.0512141261934000, -0.0190806300761000, 0.00757044398633000, -0.00325432753751000,
                 0.00149774676371000, -0.000724569558272000, 0.000362792644965000, -0.000185907973641000,
                 9.67201396593000e-05, -5.07744171678000e-05, 2.67510121456000e-05, -1.40667136728000e-05,
                 7.33363699547000e-06, -3.75638767050000e-06, 1.86344211280000e-06, -8.71623576811000e-07,
                 3.61028200288000e-07, -1.05847108097000e-07, -1.51569361490000e-08, 6.67633241420000e-08,
                 -8.33741579804000e-08, 8.31065906136000e-08, -7.53457009758000e-08, 6.48057680299000e-08,
                 -5.37558016587000e-08, 4.32436265303000e-08, -3.37262648712000e-08, 2.53558687098000e-08,
                 -1.81287021528000e-08, 1.20228328586000e-08, -7.10898040664000e-09, 3.53667004588000e-09,
                 -1.36030600198000e-09, 3.52544249042000e-10, -4.53719284366000e-11])

        else:

            raise ValueError('Digital filter order should be 0 or 1')

        # Get (complex) lambda values
        nF = filter_w.size
        ind = np.arange(nF)
        l1 = 1.0 / self.Sensor.r
        l2 = 10.0 ** (filter_a + ind * filter_s)
        L = l1 * l2
        L = L.astype('complex128')

        # Evaluate function at lambda
        YF = funname(L)

        # Calculate output, considering weights and r
        y = (np.dot(YF[None, :], filter_w[:, None])) / self.Sensor.r
        y = y[0]  # reduce dimensionality

        return y

    def reflection_coefficient(self, L):
        """
            Calculate the reflection coefficient for a given layered half-space and lambda value

            Parameters
            ----------
            L
                lambda

            Returns
            -------
            rTE
                Reflection Coefficient
        """

        # Calculate mu
        mu = mu_0 * (1.0 + self.Model.sus[None, :])

        # Calculate u0 and u
        u0 = np.sqrt(L ** 2 - self.Sensor.omega ** 2 * mu_0 * epsilon_0)
        u = np.sqrt(L[:, None] ** 2 -
                    self.Sensor.omega ** 2 * mu * self.Model.perm[None, :] +
                    1j * self.Sensor.omega * mu * self.Model.con[None, :])

        # Calculate Y and create Yhat
        Y0 = u0 / (1j * self.Sensor.omega * mu_0)
        Y = u / (1j * self.Sensor.omega * mu)
        Yhat = copy.deepcopy(Y)

        # In case of half-space
        if self.Model.con.size == 1:

            # Calculate rTE
            rTE = (Y0 - Yhat[:, 0]) / (Y0 + Yhat[:, 0])

        # In case of more complex space w0000t
        else:

            # Recursive formula
            for i in range(self.Model.nr_layers - 1)[::-1]:

                # Calculate tanh
                tanh_uh = np.tanh(u[:, i] * self.Model.thick[i])

                # Calculate Yhat
                num = Yhat[:, i + 1] + Y[:, i] * tanh_uh
                den = Y[:, i] + Yhat[:, i + 1] * tanh_uh
                Yhat[:, i] = Y[:, i] * num / den

            # Calculate rTE
            rTE = (Y0 - Yhat[:, 0]) / (Y0 + Yhat[:, 0])

        return rTE

    def propagation_matrix(self, L):
        """
            Calculates the P(2,1)/P(1,1) ratio of the propagation matrix for a given layered half-space and lambda
            value.

            :return: Propagation Matrix
        """

        # Get information
        nF = L.size
        nL = self.Model.sus.size

        # Calculate mu
        mu = mu_0 * (1.0 + self.Model.sus)

        # Calculate u0 and u
        u0 = np.sqrt(L ** 2 - self.Sensor.omega ** 2 * mu_0 * epsilon_0)
        u = np.sqrt(L[:, None] ** 2 -
                    self.Sensor.omega ** 2 * mu[None, :] * self.Model.perm[None, :] +
                    1j * self.Sensor.omega * mu[None, :] * self.Model.con[None, :])

        # Calculate M1 (for first layer)
        M = np.zeros((2, 2, nF, nL), dtype='complex128')
        M[:, :, :, 0] = np.array([[0.5 * (1 + (mu_0 * u[:, 0]) / (mu[0] * u0)),
                                   0.5 * (1 - (mu_0 * u[:, 0]) / (mu[0] * u0))],
                                  [0.5 * (1 - (mu_0 * u[:, 0]) / (mu[0] * u0)),
                                      0.5 * (1 + (mu_0 * u[:, 0]) / (mu[0] * u0))]])

        # Calculate Mn
        M[:, :, :, 1:] = np.array([[0.5 * (1 + (mu[:-1] * u[:, 1:]) / (mu[1:] * u[:, :-1])),
                                    0.5 * (1 - (mu[:-1] * u[:, 1:]) / (mu[1:] * u[:, :-1]))],
                                   [0.5 * (1 - (mu[:-1] * u[:, 1:]) / (mu[1:] * u[:, :-1]))
                                    * np.exp(-2 * u[:, :-1] * self.Model.thick[:-1]),
                                    0.5 * (1 + (mu[:-1] * u[:, 1:]) / (mu[1:] * u[:, :-1]))
                                    * np.exp(-2 * u[:, :-1] * self.Model.thick[:-1])]])

        # Dot product vector
        PM = copy.deepcopy(M[:, :, :, 0])
        for iL in range(1, nL):
            for iF in range(nF):
                PM[:, :, iF] = np.dot(PM[:, :, iF], M[:, :, iF, iL])

        # Get ratio
        PP = PM[1, 0, :] / PM[0, 0, :]

        return PP

    def rTE01(self, L):
        """
            Additional lambda functions for Hankel calculation of Primary and Secondary magnetic fields.

            :return: Output function
        """

        # Define variable u0
        u0 = np.sqrt(L ** 2 - self.Sensor.omega ** 2 * mu_0 * epsilon_0)

        # Output
        y = (np.exp(-u0 * (self.Sensor.xyz[2] + self.Sensor.height))) * L

        return y

    def rTE02(self, L):
        """
            Additional lambda functions for Hankel calculation of Primary and Secondary magnetic fields.

            :return: Output function
        """

        # Define variable u0
        u0 = np.sqrt(L ** 2 - self.Sensor.omega ** 2 * mu_0 * epsilon_0)

        # Output
        y = (np.exp(-u0 * (self.Sensor.xyz[2] + self.Sensor.height))) * L ** 3 / u0

        return y

    def rTEn1(self, L):
        """
            Additional lambda functions for Hankel calculation of Primary and Secondary magnetic fields.

            :return: Output function
        """

        # Calculate rTE
        if self.Method == 'RC':
            rTE = self.reflection_coefficient(L)
        elif self.Method == 'PM':
            rTE = self.propagation_matrix(L)

        # Define variable u0
        u0 = np.sqrt(L ** 2 - self.Sensor.omega ** 2 * mu_0 * epsilon_0)

        # Output
        y = (-rTE * np.exp(u0 * (self.Sensor.xyz[2] - self.Sensor.height))) * L

        return y

    def rTEn2(self, L):
        """
            Additional lambda functions for Hankel calculation of Primary and Secondary magnetic fields.

            :return: Output function
        """

        # Calculate rTE
        if self.Method == 'RC':
            rTE = self.reflection_coefficient(L)
        elif self.Method == 'PM':
            rTE = self.propagation_matrix(L)

        # Define variable u0
        u0 = np.sqrt(L ** 2 - self.Sensor.omega ** 2 * mu_0 * epsilon_0)

        # Output
        y = (-rTE * np.exp(u0 * (self.Sensor.xyz[2] - self.Sensor.height))) * L ** 2

        return y

    def rTEp2(self, L):
        """
            Additional lambda functions for Hankel calculation of Primary and Secondary magnetic fields.

            :return: Output function
        """

        # Calculate rTE
        if self.Method == 'RC':
            rTE = self.reflection_coefficient(L)
        elif self.Method == 'PM':
            rTE = self.propagation_matrix(L)

        # Define variable u0
        u0 = np.sqrt(L ** 2 - self.Sensor.omega ** 2 * mu_0 * epsilon_0)

        # Output
        y = (rTE * np.exp(u0 * (self.Sensor.xyz[2] - self.Sensor.height))) * L ** 3 / u0

        return y
