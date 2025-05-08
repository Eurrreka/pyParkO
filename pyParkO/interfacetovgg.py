from abc import ABCMeta, abstractmethod
from scipy.fftpack import fft2, ifft2
import numpy as np
from math import factorial


class InterfaceToVgg(metaclass=ABCMeta):
    def __init__(self, delta_bnd, delta_sigma, mu, reference_depth, longrkm, longckm):
        """
        initialize the inputs
        :param delta_bnd: A matrix stores the undulations of the density interface.
        :param delta_sigma: A float for the initial density contrast for the density interface.
        :param mu: A float for the attenuation factor of density contrast for the interface.
        :param reference_depth: A float for the reference depth of the density interface.
        :param longrkm: A numerical value representing the spherical distance (in kilometers) of
                        the ocean area being studied along the meridional direction
        :param longckm: A numerical value representing the spherical distance (in kilometers) of
                        the ocean area being studied along the vertical direction.
        """
        self.delta_bnd = delta_bnd
        self.delta_sigma = delta_sigma
        self.mu = mu
        self.reference_depth = reference_depth
        self.longrkm, self.longckm = longrkm, longckm

    @abstractmethod
    def forward(self, t):
        """
        calculate the vgg from the inputs with iteration t.
        :param t: iteration stands for the order of parker's formula.
        :return: vgg matrix.
        """
        pass


class Interface2Gravity(InterfaceToVgg):
    def __init__(self, delta_bnd, delta_sigma, mu, reference_depth, longrkm, longckm):
        super(Interface2Gravity, self).__init__(delta_bnd, delta_sigma, mu, reference_depth, longrkm, longckm)
        self.frequency = self.__frequency__()
        self.G = 6.67

    def __frequency__(self):
        """
        Inner function to calculate the frequency matrix.
        :return: Frequency matrix.
        """
        nrow, ncol = self.delta_bnd.shape
        frequency = np.zeros((nrow, ncol))
        for i in range(nrow):
            for j in range(ncol):
                ii = i if i <= nrow / 2 else i - nrow
                jj = j if j <= ncol / 2 else j - ncol
                frequency[i, j] = 2 * np.pi * np.sqrt((ii / self.longrkm) ** 2 + (jj / self.longckm) ** 2)
        return frequency

    def forward(self, t):
        # gravity factor
        factor = -2 * np.pi * self.G * self.delta_sigma * np.exp(-self.frequency * self.reference_depth)
        gravity_fft_k = 0
        for it in range(1, t + 1):
            gravity_fft_k += (self.frequency - self.mu) ** (it - 1) * fft2(self.delta_bnd ** it) / factorial(it)
        # summary gravity caused by all matrix
        gravity_fft = factor * gravity_fft_k
        # inverse iff
        gravity_fft[0, 0] = 0
        gravity_vgg = ifft2(gravity_fft).real
        return gravity_vgg


if __name__ == "__main__":
    x, y = np.linspace(0, 4 * np.pi, 20), np.linspace(0, 3 * np.pi, 30)
    xx, yy = np.meshgrid(x, y)
    zz = np.sin(xx + yy)
    miu = 0.01
    parameters = {
        "delta_bnd": zz,
        "delta_sigma": 1.82,
        "mu": miu,
        "reference_depth": 5,
        "longrkm": 500,
        "longckm": 600,
    }
    model = Interface2Gravity(**parameters)
    vgg = model.forward(t=3)
