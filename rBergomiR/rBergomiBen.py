#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -----------------------------------------------------------------------------
# IMPORTS

from math import sqrt, log, e
from scipy.special import hyp2f1
import numpy as np

# ------------------------------------------------------------------------------
# CLASS DEFINITIONS


class CholeskyPricer:
    """
    Cholesky Monte Carlo Pricer for the modified rough Bergomi model:
        $$  dS_t/S_t        = \sigma(B_t^H) d(Z_t) \\
            Z 				= \rho B + \sqrt{1-\rho^2} W \\
            \sigma(B_t^H)   = \sigma_0 * \exp{1/2 * \eta * B_t^H} \\
            B_t^H       	= \sqrt(2H) * \int_0^t (t-s)^{H-1/2} dB_t $$
    where $(W,B)$ a $2$-dim Brownian motion and $\rho \in (-1,1)$.
    :param hurst_index: hurst index of the fractional brownian motion
    :param time_steps: number of discretization steps
    :param time_to_maturity: time to expiry of the option
    :param rho: correlation parameter between driving noises
    :param spot_variance: spot variance of the model
    :param spot_price: spot price of the modeled asset
    :param vol_of_vol: volatility of volatility
    :type hurst_index: float
    :type time_steps: int
    :type time_to_maturity: float
    :type rho: float
    :type spot_variance: float
    :type spot_price: float
    :type vol_of_vol: float
    See docstrings of methods for usage.
    """

    def __init__(self, hurst_index, time_steps, time_to_maturity, rho,
                 spot_variance, spot_price, vol_of_vol):

        self.H = hurst_index
        self.N = time_steps
        self.T = time_to_maturity
        self.rho = rho
        self.v0 = spot_variance
        self.S0 = spot_price
        self.eta = vol_of_vol

        # Stores Cholesky decomposition of joint covariance matrix
        self.make_cholesky_matrix()

    def make_cholesky_matrix(self):
        """
        Builds the joint covariance matrix of $Z$ and $\hat{B}$ and then stores
        the Cholesky decomposition as an attribute of the class.
        """

        # Initialisation of variables
        # ---------------------------

        cov = np.zeros((2 * self.N, 2 * self.N))

        t_grid = np.linspace(self.T/self.N, self.T, self.N)

        gamma = 1/2 - self.H

        D_H = sqrt(2 * self.H)/(self.H + 1/2)

        # Definition of auxiliary functions
        # ---------------------------------

        # Translates index to time.
        def itt(x):
            return t_grid[x % self.N]

        # Convenience function: Long expresssion.
        def G(x):
            return 2*self.H*(1/(1-gamma)*x**(-gamma)+gamma/(1-gamma) *
                             x ** (-1-gamma) * hyp2f1(1, 1+gamma, 3-gamma,
                                                      1/x)/(2-gamma))

        # Builds covariance matrix for (Z,\hat{B}) without spot.
        # ------------------------------------------------------

        for i, j in np.ndindex(cov.shape):

            if i < self.N and j < self.N:

                cov[i, j] = min(itt(i), itt(j))

            elif i >= self.N and j >= self.N:

                if i == j:

                    cov[i, j] = itt(i)**(2*self.H)

                else:

                    max_time = max(itt(i), itt(j))

                    min_time = min(itt(i), itt(j))

                    cov[i, j] = min_time**(2*self.H)*G(max_time/min_time)

            elif i < self.N and j >= self.N:

                cov[i, j] = self.rho*D_H*(itt(j)**(self.H+1/2)-(itt(j) -
                                          min(itt(i), itt(j)))**(self.H+1/2))

            elif i >= self.N and j < self.N:

                cov[i, j] = self.rho*D_H*(itt(i)**(self.H+1/2)-(itt(i) -
                                          min(itt(i), itt(j)))**(self.H+1/2))

        # Store cholesky decomposition as attribute of the class.
        self.chol_matrix = np.linalg.cholesky(cov)

    def get_joint_paths(self, number_of_paths):
        """
        Returns required number of paths of joint process (Z, \hat{B})
        starting from spot.
        :param number_of_paths: number of required paths of joint process
                                (Z, \hat{B})
        :type number_of_paths: int
        :return: bm_paths, fbm_paths
        :rtype: numpy array
        """

        M = number_of_paths

        data = np.dot(self.chol_matrix,
                      np.random.randn(self.chol_matrix.shape[1], M))

        bm_paths = np.zeros((self.N + 1, M))

        bm_paths[1:, :] = data[:self.N, :]

        fbm_paths = np.zeros((self.N + 1, M))

        fbm_paths[1:, :] = data[self.N:, :]

        bm_paths = bm_paths.transpose()

        fbm_paths = fbm_paths.transpose()

        return bm_paths, fbm_paths

    def get_price_paths(self, number_of_paths):
        """
        Returns required number of stock price paths.
        :param number_of_paths: number of required samples from price paths
        :type number_of_paths: int
        """

        M = number_of_paths

        price_paths = np.zeros((M, self.N + 1))

        bm_paths, fbm_paths = self.get_joint_paths(M)

        vol_paths = np.sqrt(self.v0 * np.exp(self.eta * fbm_paths))

        price_paths[:, [0]] = self.S0 * np.ones((M, 1))

        S = price_paths[:, [0]]

        for j in range(1, self.N + 1):

            S += S * vol_paths[:, [j-1]] * (bm_paths[:, [j]] -
                                            bm_paths[:, [j-1]])

            price_paths[:, [j]] = S

        return price_paths

    def get_log_price_paths(self, number_of_paths):
        """
        Returns required number of log stock price paths.
        :param number_of_paths: number of required samples from log price paths
        :type number_of_paths: int
        """

        M = number_of_paths

        dt = self.T/self.N

        log_price_paths = np.zeros((M, self.N + 1))

        bm_paths, fbm_paths = self.get_joint_paths(M)

        vol_paths = np.sqrt(self.v0 * np.exp(self.eta * fbm_paths))

        log_price_paths[:, [0]] = log(self.S0) * np.ones((M, 1))

        S = log_price_paths[:, [0]]

        for j in range(1, self.N + 1):

            S += -0.5 * vol_paths[:, [j-1]]**2 * dt +  \
                 vol_paths[:, [j-1]] * (bm_paths[:, [j]] - bm_paths[:, [j-1]])

            log_price_paths[:, [j]] = S

        return log_price_paths

    def get_log_prices(self, number_of_prices):
        """
        Returns required number of log stock prices.
        :param number_of_prices: number of required samples from log price
        :type number_of_prices: int
        """

        M = number_of_prices

        dt = self.T/self.N

        bm_paths, fbm_paths = self.get_joint_paths(M)

        vol_paths = np.sqrt(self.v0 * np.exp(self.eta * fbm_paths))

        log_prices = log(self.S0) * np.ones((M, 1))

        for j in range(1, self.N + 1):

            log_prices += -0.5 * vol_paths[:, [j-1]]**2 * dt +  \
                 vol_paths[:, [j-1]] * (bm_paths[:, [j]] - bm_paths[:, [j-1]])

        return log_prices

    def get_option_price(self, flag, strike, mc_runs):
        """
        Returns European call option price based on Monte Carlo simulations.
        :param flag: either "c" for call or "p" for put
        :param strike: strike price of the European option
        :param mc_runs: number of Monte Carlo runs
        :return: option price
        :return: standard deviation of option price
        :type flag: str
        :type strike: float
        :type mc_runs: int
        :rtype: float
        :rtype: float
        """

        K = strike

        M = mc_runs

        stock_prices = np.exp(self.get_log_prices(M))

        if flag == "c":

            payoffs = np.maximum(stock_prices - K, 0)

        elif flag == 'p':

            payoffs = np.maximum(K - stock_prices, 0)

        price = np.average(payoffs)

        std_price = np.std(payoffs)/sqrt(M)

        return price, std_price


class Asymptotics:
    """
    Implementation of the various asymptotic formulae as given in the paper for
    the modified rough Bergomi model given by:
    $$  dS_t/S_t        = \sigma(B_t^H) d(Z_t) \\
        Z 				= \rho B + \sqrt{1-\rho^2} W \\
        \sigma(B_t^H)   = \sigma_0 * \exp(1/2 * \eta * B_t^H) \\
        B_t^H       	= \sqrt(2H) * \int_0^t (t-s)^{H-1/2} dB_t $$
    where $(W,B)$ a $2$-dim Brownian motion and $\rho \in (-1,1)$.
    :param hurst_index: hurst index of the fractional brownian motion
    :param log_strike_const: constant $k$ in $k_t = k * t^{1/2-H+\beta}$
    :param time_to_maturity: time to expiry of the option
    :param beta: scaling speed in $k_t = k * t^{1/2-H+\beta}$
    :param rho: correlation parameter between driving noises
    :param spot_vol: spot volatility
    :param spot_vol_prime: derivative of vol function at 0
    :type hurst_index: float
    :type log_strike_const: float
    :type time_to_maturity: float
    :type beta: float
    :type rho: float
    :type spot_vol: float
    :type spot_vol_prime: float
    """
    def __init__(self, hurst_index, log_strike_const, time_to_maturity,
                 beta, rho, spot_vol, spot_vol_prime):

        self.H = hurst_index
        self.k = log_strike_const
        self.t = time_to_maturity
        self.beta = beta
        self.rho = rho
        self.spot_vol = spot_vol
        self.spot_vol_prime = spot_vol_prime

        self.K11 = sqrt(2*self.H)/((self.H+1/2) * (self.H+3/2))

        # Computing derivatives of rate function at 0

        self.I_2d = 1/self.spot_vol**2

        self.I_3d = (-6) * self.rho * self.K11 * (self.spot_vol_prime /
                                                  (self.spot_vol**4))

    def get_log_strike(self):
        """
        Computes the formula $k_t = kt**(1/2-H+beta)$.
        """

        log_strike = self.k * self.t**(1/2 - self.H + self.beta)

        return log_strike

    def get_abs_log_call_price(self):

        if self.beta >= 2/3 * self.H and self.beta < self.H:

            print('Attention: 1st order approximation used.')

            result = 1/2 * (self.k**2) * self.I_2d * (self.t**(2*self.beta -
                                                               2*self.H))

        elif self.beta >= 0 and self.beta < 2/3 * self.H:

            secondorder = 1/2 * self.k**2 * self.I_2d * self.t**(2*self.beta -
                                                                 2*self.H)

            thirdorder = 1/6 * self.k**3 * self.I_3d * self.t**(3*self.beta -
                                                                2*self.H)

            result = secondorder + thirdorder

        else:

            print("Beta has to be smaller than H.")

        return result

    def get_implied_vol(self):

        imp_vol = self.spot_vol + self.rho * self.K11 * \
                  self.spot_vol_prime/self.spot_vol * self.get_log_strike() * \
                  self.t**(self.H-0.5)

        return imp_vol