
import numpy
import math


class triode_preamp(object):
    '''
    Common cathode triode amplifier circuit.
    
    Circuit structure & parameters as presented in the following publications by
    David Te Mao Yeh:
    "AUTOMATED PHYSICAL MODELLING OF NONLINEAR AUDIO CIRCUITS FOR REAL-TIME 
    AUDIO EFFECTS - PART II: BJT AND VACUUM TUBE EXAMPLES"
    "DIGITAL IMPLEMENTATION OF MUSICAL DISTORTION CIRCUITS BY ANALYSIS AND
    SIMULATION"

    Solution by means of Discrete K-method.
    Equations:

    x' = Ax + Bu + Ci
    i = f(v)
    v = Dx + Eu + Fi

    x : state (capacitors)
    u : inputs (input voltage and supply)
    i : nonlinear device currents
    v : nonlinear device voltages
    '''


    def __init__(self, fs, num_inputs=2, num_states=3, num_nl_voltages=2, num_iterations=5):
        '''
        Instantiate circuit object
        '''

        self.fs              = fs
        self.num_inputs      = num_inputs
        self.num_nl_voltages = num_nl_voltages
        self.num_states      = num_states
        self.num_iterations  = num_iterations

        # Circuit matrices
        self.A = numpy.zeros((self.num_states, self.num_states))
        self.B = numpy.zeros((self.num_states, self.num_inputs))
        self.C = numpy.zeros((self.num_states, self.num_nl_voltages))
        self.D = numpy.zeros((self.num_nl_voltages, self.num_states))
        self.E = numpy.zeros((self.num_nl_voltages, self.num_inputs))
        self.F = numpy.zeros((self.num_nl_voltages, self.num_nl_voltages))

        # DK method matrices
        self.H = numpy.zeros((self.num_states, self.num_states))
        self.I = numpy.identity(self.num_states)
        self.K = numpy.zeros((self.num_nl_voltages, self.num_nl_voltages))

        # initialize empty buffers
        self._clear_buffers()


    def setup(self, Rg=7e+4, Rk=15e+2, Rp=1e+5, Ri=1e+6, Ci=0.047e-6, Cf=2.5e-11, Ck=25e-6, Vpp=250):
        '''
        Setup circuit parameters
        '''

        self.Rg = Rg
        self.Rk = Rk
        self.Rp = Rp
        self.Ri = Ri
        self.Ci = Ci
        self.Cf = Cf
        self.Ck = Ck
        self.Vpp = Vpp
        self.Gg = 1.0 / self.Rg
        self.Gk = 1.0 / self.Rk
        self.Gp = 1.0 / self.Rp
        self.Gi = 1.0 / self.Ri

        self.A[0,0] = -1.0*((self.Gi+self.Gg)*self.Gp+self.Gi*self.Gg) / (self.Ci*(self.Gg+self.Gp))
        self.A[0,1] = (self.Gg*self.Gp) / (self.Ci*(self.Gg+self.Gp))
        self.A[1,0] = (self.Gg*self.Gp) / (self.Cf*(self.Gg+self.Gp))
        self.A[1,1] = -1.0*(self.Gg*self.Gp) / (self.Cf*(self.Gg+self.Gp))
        self.A[2,2] = -1.0*self.Gk/self.Ck

        self.B[0,0] = ((self.Gi+self.Gg)*self.Gp+self.Gi*self.Gg) / (self.Ci*(self.Gg+self.Gp))
        self.B[0,1] = -1.0*(self.Gg*self.Gp) / (self.Ci*(self.Gg+self.Gp))
        self.B[1,0] = -1.0*(self.Gg*self.Gp) / (self.Cf*(self.Gg+self.Gp))
        self.B[1,1] = (self.Gg*self.Gp) / (self.Cf*(self.Gg+self.Gp))

        self.C[0,0] = self.Gg / (self.Ci*(self.Gg+self.Gp))
        self.C[0,1] = self.Gg / (self.Ci*(self.Gg+self.Gp))
        self.C[1,0] = self.Gp / (self.Cf*(self.Gg+self.Gp))
        self.C[1,1] = -1.0*self.Gp / (self.Cf*(self.Gg+self.Gp))
        self.C[2,0] = 1.0 / self.Ck
        self.C[2,1] = 1.0 / self.Ck

        self.D[0,0] = -1.0*self.Gg / (self.Gg+self.Gp)
        self.D[0,1] = -1.0*self.Gp / (self.Gg+self.Gp)
        self.D[0,2] = -1.0
        self.D[1,0] = -1.0*self.Gg / (self.Gg+self.Gp)
        self.D[1,1] = self.Gg / (self.Gg+self.Gp)
        self.D[1,2] = -1.0

        self.E[0,0] = self.Gg / (self.Gg+self.Gp)
        self.E[0,1] = self.Gp / (self.Gg+self.Gp)
        self.E[1,0] = self.Gg / (self.Gg+self.Gp)
        self.E[1,1] = self.Gp / (self.Gg+self.Gp)

        self.F[0,0] = -1.0 / (self.Gp+self.Gg)
        self.F[0,1] = -1.0 / (self.Gp+self.Gg)
        self.F[1,0] = -1.0 / (self.Gp+self.Gg)
        self.F[1,1] = -1.0 / (self.Gp+self.Gg)


    def discretize(self):
        '''
        Discretize circuit (trapezoidal rule)
        '''

        # precalculate matrices needed for parameter & state update
        self.alpha = 2.0*self.fs
        self.H     = numpy.linalg.inv(self.alpha*self.I - self.A)
        self.DH    = numpy.dot(self.D, self.H)
        self.K     = numpy.dot(self.DH, self.C) + self.F
        self.HB    = numpy.dot(self.H, self.B)
        self.HC    = numpy.dot(self.H, self.C)
        self.aIA   = self.alpha*self.I+self.A
        self.HaIA  = numpy.dot(self.H, self.aIA)


    def process(self, input_sample, tube):
        '''
        Process single sample
        '''

        # input voltages: input & supply rail
        self.un[0,0] = input_sample
        self.un[1,0] = self.Vpp

        # update state parameter p[n]
        self._update_state_param()

        # Newton solver
        for k in range(0, self.num_iterations):

            # get non-linear currents and derivatives
            tube.calc_currents(self.est[0,0], self.est[1,0])
            self.it[0,0] = tube.Ig
            self.it[1,0] = tube.Ip

            # equation to be solved
            FV = self.pn + numpy.dot(self.K, self.it) - self.est

            # Jacobian
            self.J[0,0] = self.K[0,0]*tube.Ig_Vgk + self.K[0,1]*tube.Ip_Vgk - 1.0
            self.J[0,1] = self.K[0,0]*tube.Ig_Vpk + self.K[0,1]*tube.Ip_Vpk
            self.J[1,0] = self.K[1,0]*tube.Ig_Vgk + self.K[1,1]*tube.Ip_Vgk
            self.J[1,1] = self.K[1,0]*tube.Ig_Vpk + self.K[1,1]*tube.Ip_Vpk - 1.0

            # avoid NaNs
            self.J[self.J==0] = 1e-10

            # solve
            self.est = self.est - numpy.dot(numpy.linalg.inv(self.J),FV)

        # update currents
        tube.calc_currents(self.est[0,0], self.est[1,0])
        self.it[0,0] = tube.Ig
        self.it[1,0] = tube.Ip

        # update state
        self._update_state()

        # write buffers
        self.prev_x = self.xn
        self.prev_i = self.it
        self.prev_u = self.un

        # get output voltages:
        # capacitors
        self.Vci  = self.xn[0,0]
        self.Vcf  = self.xn[1,0]
        self.Vck  = self.xn[2,0]
        # tube
        self.Vgk  = self.est[0,0]
        self.Vpk  = self.est[1,0]
        # output:
        self.Vout = self.Vpk + self.Vck


    def _clear_buffers(self):
        '''
        Reset buffers
        '''

        self.prev_x = numpy.zeros((self.num_states,1))
        self.prev_i = numpy.zeros((self.num_nl_voltages,1))
        self.prev_u = numpy.zeros((self.num_inputs,1))
        self.est    = numpy.zeros((self.num_nl_voltages,1))
        self.un     = numpy.zeros((self.num_inputs,1))
        self.xn     = numpy.zeros((self.num_states,1))
        self.pn     = numpy.zeros((self.num_nl_voltages,1))
        self.it     = numpy.zeros((self.num_nl_voltages,1))
        self.J      = numpy.zeros((self.num_nl_voltages, self.num_nl_voltages))


    def _update_state_param(self):
        '''
        Update state parameter p[n]
        '''

        # trapezoidal rule: 
        # p[n] = DH((aI+A)x[n-1] + B(u[n]) + u[n-1]) + Ci[n-1]) + Eu[n]
        temp = numpy.dot(self.aIA, self.prev_x) + numpy.dot(self.B,( self.un + self.prev_u)) \
                + numpy.dot(self.C, self.prev_i)
        self.pn = numpy.dot(self.DH, temp) + numpy.dot(self.E,  self.un)


    def _update_state(self):
        '''
        Update system state x[n]
        '''

        # trapezoidal rule: 
        # x[n] = H(aI+A)x[n-1] + HB(u[n]) + u[n-1]) + HC(i[n]+i[n-1])
        self.xn = numpy.dot(self.HaIA, self.prev_x) + numpy.dot(self.HB, self.un + self.prev_u) \
                + numpy.dot(self.HC, self.it + self.prev_i)

