
import math


class triode_model(object):
    '''
    12AX7 triode non-linear element model as presented in the paper:
    "A PHYSICALLY-MOTIVATED TRIODE MODEL FOR CIRCUIT SIMULATIONS"
    by Kristjan Dempwolf and Udo Zoelzer, DAFX-11

    Ig     : grid current
    Ip     : plate current
    Ig_Vgk : grid current derivative to grid-cathode voltage
    Ig_Vpk : grid current derivative to plate-cathode voltage
    Ip_Vgk : plate current derivative to grid-cathode voltage
    Ip_Vpk : plate current derivative to plate-cathode voltage
    '''


    def __init__(self, G=1.371e-3, mu=86.9, gamma=1.349, C=4.56, Gg=3.263e-4, xi=1.456, Cg=11.99, Ig0=3.917e-8):
        '''
        Instantiate a triode model.
        Default parameters are for 12AX7 EHX-1.
        '''

        # tube parameters
        self.G     = G
        self.mu    = mu
        self.gamma = gamma
        self.C     = C
        self.Gg    = Gg
        self.xi    = xi
        self.Cg    = Cg
        self.Ig0   = Ig0

        # currents (grid / plate)
        self.Ig     = 1e-12
        self.Ip     = 1e-12
        self.Ig_Vgk = 1e-12
        self.Ig_Vpk = 1e-12
        self.Ip_Vgk = 1e-12
        self.Ip_Vpk = 1e-12


    def calc_currents(self, Vgk, Vpk):
        '''
        Calculate grid/plate currents and derivatives.
        '''

        # Grid current
        exg     = math.exp(self.Cg*Vgk)
        logexg  = math.log(exg + 1.0)
        temp    = (1.0/self.Cg)*math.log(1.0 + exg)
        self.Ig = self.Gg*math.pow(temp, self.xi) + self.Ig0

        # Grid current derivatives
        if (exg + 1.0)*logexg == 0:
            self.Ig_Vgk = 1e-12
        else:
            self.Ig_Vgk = self.Cg*self.Gg*self.xi*(logexg/self.Cg)*exg / ((exg + 1.0)*logexg)
        self.Ig_Vpk = 1e-12

        # Plate current
        ex      = math.exp(self.C*Vgk + self.C*Vpk/self.mu)
        logex   = math.log(ex + 1.0)
        self.Ip = self.G*math.pow((1.0/self.C)*logex, self.gamma)

        # Plate current derivatives
        num = self.C*self.G*self.gamma*math.pow(logex/self.C, self.gamma)*ex
        den = (ex + 1.0)*logex
        if den == 0:
            self.Ip_Vgk = 1e-12
            self.Ip_Vpk = 1e-12
        else:
            self.Ip_Vgk = num /den
            self.Ip_Vpk = num / (self.mu*den)


