from scipy.fftpack import fft2, ifft2
import numpy as np
from abc import ABCMeta, abstractmethod
from math import factorial
from scipy.signal.windows import tukey


class VggToInterface(metaclass=ABCMeta):
    def __init__(self, vgg, delta_sigma, mu, reference_depth, longrkm, longckm):
        """
        initialize the inputs and parameters
        :param vgg: A matrix for gravity anomalies.
        :param delta_sigma: A float for the density contrast of the density interface.
        :param mu: A float for the attenuation factor of density contrasts for the interfaces.
        :param reference_depth: A float for the reference depth of the density interfaces.
        :param longrkm: A numerical value representing the spherical distance (in kilometers) of
                        the ocean area being studied along the meridional direction
        :param longckm: A numerical value representing the spherical distance (in kilometers) of
                        the ocean area being studied along the vertical direction.
        """
        self.vgg = vgg
        self.delta_sigma = delta_sigma
        self.mu = mu
        self.reference_depth = reference_depth
        self.longrkm, self.longckm = longrkm, longckm

    @abstractmethod
    def downward(self, t, criteria):
        """
        calculate the density interface of interest with downward iteration steps.
        :param t: iteration for downwards.
        :param criteria: criteria for downwards iteration.
        :return: matrix for undulation of density interface of interest.
        """
        pass

    @classmethod
    def twkey(cls, matrix, edge=0.02):
        """
        tukey the border values for smoothness with outside values 0.
        :param matrix: matrix for tukey.
        :param edge: float for definition of border.
        :return: tukey matrix.
        """
        nrow, ncol = matrix.shape
        tky = np.array([row_tky * col_tky for row_tky in tukey(nrow, edge)
                        for col_tky in tukey(ncol, edge)]).reshape(matrix.shape)
        return tky * matrix


class Gravity2Interface(VggToInterface):
    def __init__(self, vgg, delta_sigma, mu, reference_depth, longrkm, longckm, wh, alpha):
        """
        initialize the inputs and parameters for Gravity2Interface.
        :param vgg: A matrix for gravity anomalies.
        :param delta_sigma: A float for the density contrast of the density interface.
        :param mu: A float for the attenuation factor of density contrasts for the interfaces.
        :param reference_depth: A float for the reference depth of the density interfaces.
        :param longrkm: A numerical value representing the spherical distance (in kilometers) of
                        the ocean area being studied along the meridional direction
        :param longckm: A numerical value representing the spherical distance (in kilometers) of
                        the ocean area being studied along the vertical direction.
        :param wh: A threshold ranging from 0.1 to 0.4 represents the cutoff point utilized
                    for eliminating high-frequency components.
        :param alpha: An integer ranging from 4 to 32 representing the penalty imposed on
                    high-frequency components of Fourier transforms.
        """
        super(Gravity2Interface, self).__init__(vgg, delta_sigma, mu, reference_depth, longrkm, longckm)
        self.G = 6.67
        self.wh = wh
        self.alpha = alpha
        self.frequency = self.__frequency__()
        self.filter = self.__filter__()
        self.delta_bnd = None
        # temp
        self.temp = {}

    def __frequency__(self):
        """
        Inner function for calculating the frequency matrix.
        :return: Frequency matrix
        """
        nrow, ncol = self.vgg.shape
        frequency = np.zeros((nrow, ncol))
        for i in range(nrow):
            for j in range(ncol):
                ii = i if i <= nrow / 2 else i - nrow
                jj = j if j <= ncol / 2 else j - ncol
                frequency[i, j] = 2 * np.pi * np.sqrt((ii / self.longrkm) ** 2 + (jj / self.longckm) ** 2)
        return frequency

    def __filter__(self):
        """
        Inner function for calculating the filter that lowpass the frequency value
        :return: A filter matrix
        """
        nrow, ncol = self.vgg.shape
        filter = np.ones(self.vgg.shape)
        for i in range(nrow):
            for j in range(ncol):
                if self.frequency[i, j] > self.wh:
                    ratio = self.frequency[i, j] / self.wh
                    filter[i, j] = ratio ** (1 - self.alpha) - (1 - self.alpha) * np.log(ratio) * ratio ** (
                                1 - self.alpha)
        return filter

    def bnd_n(self, n):
        bnd_twkey = self.twkey(self.delta_bnd)
        bnd_twkey_n = bnd_twkey ** n
        bnd_fourier = fft2(bnd_twkey_n)
        return bnd_fourier

    def __vgg_fft__(self):
        name = 'vgg-ft'
        if self.temp.get(name) is None:
            vgg_ft = fft2(self.twkey(self.vgg))
            self.temp[name] = vgg_ft
        else:
            vgg_ft = self.temp.get(name)
        return vgg_ft

    def once_downward(self, t):
        fft = self.__vgg_fft__() / 2 / np.pi / self.G / self.delta_sigma / \
              np.exp(-self.frequency * self.reference_depth)
        summary = 0 + 0j
        for n in range(2, t + 1):
            summary += (self.frequency - self.mu) ** (n-1) / factorial(n) * self.bnd_n(n)
        residual_fft = (fft - summary) * self.filter
        residual_fft[0, 0] = 0 + 0j
        self.delta_bnd = ifft2(residual_fft).real

    def downward(self, t=10, criteria=None):
        # calculate the initialized interface
        self.once_downward(t=1)
        temp_interface = self.delta_bnd
        if t == 1:
            return self.delta_bnd
        else:
            for n in range(2, t+1):
                self.once_downward(t=n)
                if criteria is not None:
                    rmse = np.sqrt(np.mean((self.delta_bnd - temp_interface) ** 2))
                    if rmse < criteria:
                        print(
                            f'The error of the first {n} terms and the first {n - 1} terms: {rmse} '
                            f'is less than the threshold {criteria} '
                            f'\n -------------Iteration Completed!!!------------')
                        print(f'Total number of iterations: {n}')
                        break
                    else:
                        print(f'The error of the first {n} terms and the first {n - 1} terms: {rmse}')
                        temp_interface = self.delta_bnd
            return self.delta_bnd - self.reference_depth


if __name__ == "__main__":
    x, y = np.linspace(0, 4 * np.pi, 20), np.linspace(0, 3 * np.pi, 30)
    xx, yy = np.meshgrid(x, y)
    zz = np.sin(xx + yy)
    miu = 0.01
    parameters = {
        "vgg": zz,
        "delta_sigma": 1.82,
        "mu": miu,
        "reference_depth": 5,
        "longrkm": 500,
        "longckm": 600,
        "wh": 0.18,
        "alpha": 8,
    }
    model = Gravity2Interface(**parameters)
    inverse_topography = model.downward(t=8)
