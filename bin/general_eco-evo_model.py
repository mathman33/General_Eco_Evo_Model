from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, exp
from scipy.integrate import odeint
from timeit import default_timer
from mpl_toolkits.mplot3d import Axes3D
import gc as garbage
import os
import json
import argparse
import shlex
from time import sleep
from subprocess import Popen, PIPE
from datetime import datetime

LINESTYLES = ['--', '-', ':', '_']

AVG_TIME_PER_GRAPH = 0.583
SUPPORTED_DIMENSIONS = ["1x1", "1x2", "2x1"]
DIRECTORY = os.path.dirname(os.path.abspath(__file__))
if DIRECTORY.endswith("/bin"):
    DIRECTORY = DIRECTORY[:-4]
IMAGEMAGICK_SIDE_BY_SIDE_COMMAND = "convert +append %s %s %s"
IMAGEMAGICK_ON_TOP_COMMAND = "convert -append %s %s %s"
NO_IMAGEMAGICK_ERROR = "Unable to create dual graphs.  Please install 'ImageMagick'\n"
DATE_TIME_DIRECTORY_FORMAT = '%y%m%d_%H%M%S'
TIME_NEEDED_MESSAGE = "Approximate Time Needed: %.03f minutes\n\n"
NUMBER_OF_GRAPHS_MESSAGE = "\n\n%d graphs will be generated."
AVG_TIME_MESSAGE = "average time per graph: %.03f seconds"
TOT_TIME_MESSAGE = "total time taken: %.03f seconds"
NO_STABILITY_CHECK = """
NO STABILITY CHECK FILE EXISTS --- NO ANALYSIS WILL BE DONE
"""
F_OPTION_HELP = """\
Final values are displayed on the time graphs."""
K_OPTION_HELP = """\
Individual Density and Trait graphs are not deleted after being combined into a single .png
through ImageMagick.  This option is automatically enabled if -c/--no-combine is specified."""
N_OPTION_HELP = """\
ImageMagick is not called, and no combination .png's are created.  This option automatically
enables the -k/--keep-orignial-images option."""
P_OPTION_HELP = """\
No parameter values are printed on the graphs."""
L_OPTION_HELP = """\
No legends are included on the graphs."""
LOWER_LIMIT_HELP = """
Define lowest y-value on the trait graph.  Default is -5.
"""
UPPER_LIMIT_HELP = """
Define highest y-value on the trait graph.  Default is 5.
"""
O_OPTION_HELP = """
Include black horizontal lines indicating optimal prey traits in the traits vs time graphs.
"""
GRAPH_SAVED = """
GRAPH SAVED
-----------
%s"""
LaTeX_VARIABLE_FORMAT = {
    "M0":     "M_{0",
    "m0":     "m_{0",
    "sigma":  "\\sigma_{",
    "sigmaG": "\\sigma_{G",
    "d":      "d_{",

    "N0":     "N_{0",
    "n0":     "n_{0",
    "beta":   "\\beta_{",
    "betaG":  "\\beta_{G",
    "r":      "r_{",
    "K":      "K_{",

    "eff":    "e_{",
    "tau":    "\\tau_{",
    "alpha":  "\\alpha_{",
    "theta":  "\\theta_{",

    "rho":    "\\rho_{",
    "phi":    "\\phi_{",
    "gamma":  "\\gamma_{",

    "kappa":  "\\kappa_{",
    "xi":     "\\xi_{",
    "delta":  "\\delta_{"
}


def header(string):
    l = len(string)
    hashtags = "\n" + "#"*(3*l) + "\n"
    total = hashtags + " "*l + string + hashtags
    return total


def subheader(string):
    l = len(string)
    hashtags = "\n" + " "*4 + "#"*l + "\n"
    total = hashtags + " "*4 + string + hashtags
    return total


def remove_command(*items):
    command = "rm"
    for item in items:
        command += " %s" % item
    return command


def plot_densities(system, densities_file, text, args):
    plt.figure()

    if args.display_parameters:
        plt.axes([0.20, 0.1, 0.75, 0.8], axisbg="white", frameon=True)

    limit = 1.1*max([max(system.M[value]) for value in system.M] + [max(system.N[value]) for value in system.N])

    plt.ylim(-1., limit)
    plt.xlabel('Time')
    plt.ylabel('Population Density')

    if args.display_parameters:
        for index, text_line in enumerate(text):
            plt.text(-.25*system.tf, limit*(1-(.05*index)), text_line)

    LINESTYLE_NO = 0

    for i, value in enumerate(system.M):
        plt.plot(system.t, system.M[value], LINESTYLES[LINESTYLE_NO], label="Predator %d Density" % (i+1), lw=2)
        if args.final_values:
            plt.text(system.t[-1]*1.01, system.M[value][-1], "%.3f" % system.M[value][-1])
        LINESTYLE_NO += 1
    for i, value in enumerate(system.N):
        plt.plot(system.t, system.N[value], LINESTYLES[LINESTYLE_NO], label="Prey %d Density" % (i+1), lw=2)
        if args.final_values:
            plt.text(system.t[-1]*1.01, system.N[value][-1], "%.3f" % system.N[value][-1])
        LINESTYLE_NO += 1

    if args.legend:
        plt.legend(loc=0)

    plt.savefig(densities_file, format='png', dpi=200)
    plt.close()
    garbage.collect()
    print GRAPH_SAVED % densities_file


def plot_traits(system, traits_file, text, args):
    plt.figure()

    if args.display_parameters and not args.combine:
        plt.axes([0.20, 0.1, 0.75, 0.8], axisbg="white", frameon=True)

    plt.ylim(LOWER_LIMIT, UPPER_LIMIT)
    plt.xlabel('Time')
    plt.ylabel('Mean Trait Value')

    if args.display_parameters and not args.combine:
        for index, text_line in enumerate(text):
            plt.text(-.25*system.tf, UPPER_LIMIT*(1-(.05*index)), text_line)

    LINESTYLE_NO = 0

    for i, value in enumerate(system.m):
        plt.plot(system.t, system.m[value], LINESTYLES[LINESTYLE_NO], label="Predator %d Mean Trait Value" % (i+1), lw=2)
        if args.final_values:
            plt.text(system.t[-1]*1.01, system.m[value][-1], "%.3f" % system.m[value][-1])
        LINESTYLE_NO += 1
    for i, value in enumerate(system.n):
        plt.plot(system.t, system.n[value], LINESTYLES[LINESTYLE_NO], label="Prey %d Mean Trait Value" % (i+1), lw=2)
        if args.final_values:
            plt.text(system.t[-1]*1.01, system.n[value][-1], "%.3f" % system.n[value][-1])
        LINESTYLE_NO += 1

    if args.opt_traits_horiz_line:
        if system.__class__.__name__ == "VariableGrowthSystem":
            for prey, phi in system.phi.iteritems():
                plt.axhline(phi, color="k")
        elif system.__class__.__name__ == "VariableCapacitySystem":
            for prey, xi in system.xi.iteritems():
                plt.axhline(xi, color="k")
        elif system.__class__.__name__ == "FullSystem":
            for prey, phi in system.phi.iteritems():
                plt.axhline(phi, color="k")
            for prey, xi in system.xi.iteritems():
                plt.axhline(xi, color="k", ls="dashed")

    if args.legend:
        plt.legend(loc=0)

    plt.savefig(traits_file, format='png', dpi=200)
    plt.close()
    garbage.collect()
    print GRAPH_SAVED % traits_file


### 1x1 HARD-CODED PHASE PLANE FUNCTION
def plot_1x1_densities_phase_plane(system, phase_plane_file, text, args):
    plt.figure()

    if args.display_parameters and not args.combine:
        plt.axes([0.20, 0.1, 0.75, 0.8], axisbg="white", frameon=True)

    xlimit = 1.1*max(system.M["1"])
    ylimit = 1.1*max(system.N["1"])

    plt.xlim(-1., xlimit)
    plt.ylim(-1., ylimit)
    plt.xlabel('Predator 1 Density')
    plt.ylabel('Prey 1 Density')

    if args.display_parameters and not args.combine:
        for index, text_line in enumerate(text):
            plt.text(-.25*xlimit, ylimit*(1-(.05*index)), text_line)

    plt.plot(system.M["1"], system.N["1"], lw=1)
    plt.plot(system.M["1"][0], system.N["1"][0], 'g^', label="TIME=0.0")
    plt.plot(system.M["1"][-1], system.N["1"][-1], 'ro', label="TIME=%.1f" % system.tf)

    if args.legend:
        plt.legend(loc=0)

    plt.savefig(phase_plane_file, format='png', dpi=200)
    plt.close()
    garbage.collect()
    print GRAPH_SAVED % phase_plane_file


### 1x2 HARD-CODED PHASE PLANE FUNCTION
def plot_1x2_densities_phase_plane(system, phase_plane_file, text, args):
    fig = plt.figure()

    if args.display_parameters and not args.combine:
        plt.axes([0.20, 0.1, 0.75, 0.8], axisbg="white", frameon=True)

    ax = fig.gca(projection='3d')

    xlimit = 1.1*max(system.M["1"])
    ylimit = 1.1*max(system.N["1"])
    zlimit = 1.1*max(system.N["2"])

    if args.display_parameters and not args.combine:
        for index, text_line in enumerate(text):
            plt.text(-.25*xlimit, ylimit*(1-(.05*index)), text_line)

    ax.set_xlim3d(0, xlimit)
    ax.set_ylim3d(0, ylimit)
    ax.set_zlim3d(0, zlimit)

    ax.set_xlabel('Predator Density')
    ax.set_ylabel('Prey 1 Density')
    ax.set_zlabel('Prey 2 Density')

    ax.plot(system.M["1"], system.N["1"], system.N["2"])

    plt.savefig(phase_plane_file, format='png', dpi=200)
    plt.close()
    garbage.collect()
    print GRAPH_SAVED % phase_plane_file


### 2x1 HARD-CODED PHASE PLANE FUNCTION
def plot_2x1_densities_phase_plane(system, phase_plane_file, text, args):
    fig = plt.figure()

    if args.display_parameters and not args.combine:
        plt.axes([0.20, 0.1, 0.75, 0.8], axisbg="white", frameon=True)

    ax = fig.gca(projection='3d')

    xlimit = 1.1*max(system.M["1"])
    ylimit = 1.1*max(system.M["2"])
    zlimit = 1.1*max(system.N["1"])

    if args.display_parameters and not args.combine:
        for index, text_line in enumerate(text):
            plt.text(-.25*xlimit, ylimit*(1-(.05*index)), text_line)

    ax.set_xlim3d(0, xlimit)
    ax.set_ylim3d(0, ylimit)
    ax.set_zlim3d(0, zlimit)

    ax.set_xlabel('Predator 1 Density')
    ax.set_ylabel('Predator 2 Density')
    ax.set_zlabel('Prey 1 Density')

    ax.plot(system.M["1"], system.M["2"], system.N["1"])

    plt.savefig(phase_plane_file, format='png', dpi=200)
    plt.close()
    garbage.collect()
    print GRAPH_SAVED % phase_plane_file


### 1x1 HARD-CODED PHASE PLANE FUNCTION
def plot_1x1_traits_phase_plane(system, phase_plane_file, text, args):
    plt.figure()

    if args.display_parameters and not args.combine:
        plt.axes([0.20, 0.1, 0.75, 0.8], axisbg="white", frameon=True)

    m_min = min(system.m["1"])
    m_max = max(system.m["1"])

    m_range = m_max - m_min

    m_min -= (0.05*m_range)
    m_max += (0.05*m_range)

    n_min = min(system.n["1"])
    n_max = max(system.n["1"])

    n_range = n_max - n_min

    n_min -= (0.05*n_range)
    n_max += (0.05*n_range)

    plt.xlim(m_min, m_max)
    plt.ylim(n_min, n_max)

    plt.xlabel('Predator 1 Mean Trait Value')
    plt.ylabel('Prey 1 Mean Trait Value')

    if args.display_parameters and not args.combine:
        for index, text_line in enumerate(text):
            x_diff = m_max - m_min
            plt.text(m_min-.25*x_diff, n_max*(1-(.05*index)), text_line)

    plt.plot(system.m["1"], system.n["1"], lw=1)
    plt.plot(system.m["1"][0], system.n["1"][0], 'g^', label="TIME=0.0")
    plt.plot(system.m["1"][-1], system.n["1"][-1], 'ro', label="TIME=%.1f" % system.tf)

    if args.legend:
        plt.legend(loc=0)

    plt.savefig(phase_plane_file, format='png', dpi=200)
    plt.close()
    garbage.collect()
    print GRAPH_SAVED % phase_plane_file


### 1x2 HARD-CODED PHASE PLANE FUNCTION
def plot_1x2_traits_phase_plane(system, phase_plane_file, text, args):
    fig = plt.figure()

    if args.display_parameters and not args.combine:
        plt.axes([0.20, 0.1, 0.75, 0.8], axisbg="white", frameon=True)

    ax = fig.gca(projection='3d')

    m_min = min(system.m["1"])
    m_max = max(system.m["1"])

    m_range = m_max - m_min

    m_min -= (0.05*m_range)
    m_max += (0.05*m_range)

    n1_min = min(system.n["1"])
    n1_max = max(system.n["1"])

    n1_range = n1_max - n1_min

    n1_min -= (0.05*n1_range)
    n1_max += (0.05*n1_range)

    n2_min = min(system.n["2"])
    n2_max = max(system.n["2"])

    n2_range = n2_max - n2_min

    n2_min -= (0.05*n2_range)
    n2_max += (0.05*n2_range)

    if args.display_parameters and not args.combine:
        for index, text_line in enumerate(text):
            x_diff = m_max - m_min
            plt.text(m_min-.25*x_diff, n1_max*(1-(.05*index)), text_line)

    ax.set_xlim3d(m_min, m_max)
    ax.set_ylim3d(n1_min, n1_max)
    ax.set_zlim3d(n2_min, n2_max)

    ax.set_xlabel('Predator Mean Trait Value')
    ax.set_ylabel('Prey 1 Mean Trait Value')
    ax.set_zlabel('Prey 2 Mean Trait Value')

    ax.plot(system.m["1"], system.n["1"], system.n["2"])

    plt.savefig(phase_plane_file, format='png', dpi=200)
    plt.close()
    garbage.collect()
    print GRAPH_SAVED % phase_plane_file


### 2x1 HARD-CODED PHASE PLANE FUNCTION
def plot_2x1_traits_phase_plane(system, phase_plane_file, text, args):
    fig = plt.figure()

    if args.display_parameters and not args.combine:
        plt.axes([0.20, 0.1, 0.75, 0.8], axisbg="white", frameon=True)

    ax = fig.gca(projection='3d')

    m1_min = min(system.m["1"])
    m1_max = max(system.m["1"])

    m1_range = m1_max - m1_min

    m1_min -= (0.05*m1_range)
    m1_max += (0.05*m1_range)

    m2_min = min(system.m["2"])
    m2_max = max(system.m["2"])

    m2_range = m2_max - m2_min

    m2_min -= (0.05*m2_range)
    m2_max += (0.05*m2_range)

    n_min = min(system.n["1"])
    n_max = max(system.n["1"])

    n_range = n_max - n_min

    n_min -= (0.05*n_range)
    n_max += (0.05*n_range)

    if args.display_parameters and not args.combine:
        for index, text_line in enumerate(text):
            x_diff = m1_max - m1_min
            plt.text(m1_min-.25*x_diff, n_max*(1-(.05*index)), text_line)

    ax.set_xlim3d(m1_min, m1_max)
    ax.set_ylim3d(m2_min, m2_max)
    ax.set_zlim3d(n_min, n_max)

    ax.set_xlabel('Predator 1 Mean Trait Value')
    ax.set_ylabel('Predator 2 Mean Trait Value')
    ax.set_zlabel('Prey 1 Mean Trait Value')

    ax.plot(system.m["1"], system.m["2"], system.n["1"])

    plt.savefig(phase_plane_file, format='png', dpi=200)
    plt.close()
    garbage.collect()
    print GRAPH_SAVED % phase_plane_file


def combine_images(input1, input2, output, keep_original_images, COMMAND):
    command = COMMAND % (input1, input2, output)
    try:
        Popen(shlex.split(command))
        if not keep_original_images:
            while True:
                if os.path.isfile(output):
                    ts = default_timer()
                    while True:
                        output_size = float(os.path.getsize(output))
                        if default_timer() - ts > 10:
                            break
                        sleep(.005)
                        if output_size != 0:
                            Popen(shlex.split(remove_command(input1, input2)))
                            break
                    break
                else:
                    pass
    except:
        print(NO_IMAGEMAGICK_ERROR)


def irange(start, stop, step):
    r = start
    while r <= stop:
        yield r
        r += step


def give_params_avgattack(system, interaction_subscript):
    A = system.A[interaction_subscript]
    alpha = system.alpha[interaction_subscript]
    tau = system.tau[interaction_subscript]
    theta = system.theta[interaction_subscript]

    def avgattack(m, n):
        return ((alpha*tau)/sqrt(A))*exp(-(m - n - theta)**2/(2*A))

    return avgattack


def give_params_avg_pred_fitness(system, pred_subscript):
    d = system.d[pred_subscript]

    def avg_pred_fitness(m, N, n):
        fitness_source = 0
        for prey_subscript in system.N0:
            interaction_subscript = pred_subscript + prey_subscript
            eff = system.eff[interaction_subscript]
            avgattack = system.avgattack[interaction_subscript]
            fitness_source += eff*avgattack(m, n[int(prey_subscript)-1])*N[int(prey_subscript)-1]
        return fitness_source - d

    return avg_pred_fitness


def give_params_pred_trait_response(system, pred_subscript):

    def pred_trait_response(N, m, n):
        response = 0
        for prey_subscript in system.N0:
            interaction_subscript = pred_subscript + prey_subscript
            avgattack = system.avgattack[interaction_subscript]
            eff = system.eff[interaction_subscript]
            A = system.A[interaction_subscript]
            theta = system.theta[interaction_subscript]
            response += avgattack(m, n[int(prey_subscript)-1])*(eff*N[int(prey_subscript)-1]*(theta + n[int(prey_subscript)-1] - m))/(A)
        return response

    return pred_trait_response


def give_params_avg_prey_growth_rate(system, prey_subscript):
    rho = system.rho[prey_subscript]
    phi = system.phi[prey_subscript]
    gamma = system.gamma[prey_subscript]
    B = system.B[prey_subscript]

    def avg_prey_growth_rate(n):
        numerator = rho*gamma
        denominator = sqrt(B)
        exponent_num = -(n - phi)**2
        exponent_denom = 2*B
        return (numerator/denominator)*exp(exponent_num/exponent_denom)

    return avg_prey_growth_rate


class FullSystem:
    def __init__(self, parameters, tf):
        self.tf = tf
        self.t = np.linspace(0, self.tf, 100000)

        self.M0 = parameters["M0"]
        self.m0 = parameters["m0"]
        self.sigma = parameters["sigma"]
        self.sigmaG = parameters["sigmaG"]
        self.d = parameters["d"]
        self.num_preds = len(self.M0)

        self.N0 = parameters["N0"]
        self.n0 = parameters["n0"]
        self.beta = parameters["beta"]
        self.betaG = parameters["betaG"]
        self.num_preys = len(self.N0)

        self.eff = parameters["eff"]
        self.tau = parameters["tau"]
        self.alpha = parameters["alpha"]
        self.theta = parameters["theta"]

        self.rho = parameters["rho"]
        self.phi = parameters["phi"]
        self.gamma = parameters["gamma"]

        self.kappa = parameters["kappa"]
        self.xi = parameters["xi"]
        self.delta = parameters["delta"]

        self.y0 = []
        for i in xrange(0, self.num_preds):
            self.y0.append(self.M0[str(i+1)])
        for i in xrange(0, self.num_preys):
            self.y0.append(self.N0[str(i+1)])
        for i in xrange(0, self.num_preds):
            self.y0.append(self.m0[str(i+1)])
        for i in xrange(0, self.num_preys):
            self.y0.append(self.n0[str(i+1)])

        self.A = {}
        for pred_subscript in self.M0:
            for prey_subscript in self.N0:
                interaction_subscript = pred_subscript + prey_subscript
                self.A[interaction_subscript] = self.sigma[pred_subscript]**2 + self.beta[prey_subscript]**2 + self.tau[interaction_subscript]**2

        self.B = {}
        for prey_subscript in self.N0:
            self.B[prey_subscript] = self.beta[prey_subscript]**2 + self.gamma[prey_subscript]**2

        self.D = {}
        for prey_subscript in self.N0:
            self.D[prey_subscript] = (self.delta[prey_subscript]**2)*(self.beta[prey_subscript]**2)
            self.D[prey_subscript] += (self.gamma[prey_subscript]**2)*(self.delta[prey_subscript]**2)
            self.D[prey_subscript] -= (self.beta[prey_subscript]**2)*(self.gamma[prey_subscript]**2)

        self.E = {}
        for prey_subscript in self.N0:
            self.E[prey_subscript] = (self.delta[prey_subscript]**2) - (((self.beta[prey_subscript]*self.gamma[prey_subscript])**2)/(self.B[prey_subscript]))

        self.avgattack = {}
        for pred_subscript in self.M0:
            for prey_subscript in self.N0:
                interaction_subscript = pred_subscript + prey_subscript
                self.avgattack[interaction_subscript] = give_params_avgattack(self, interaction_subscript)

        self.avg_pred_fitness = {}
        self.pred_trait_response = {}
        for pred_subscript in self.M0:
            self.avg_pred_fitness[pred_subscript] = give_params_avg_pred_fitness(self, pred_subscript)
            self.pred_trait_response[pred_subscript] = give_params_pred_trait_response(self, pred_subscript)

        self.avg_prey_fitness = {}
        self.prey_trait_response = {}
        self.avg_prey_growth_rate = {}
        self.capacity_hat = {}
        for prey_subscript in self.N0:
            self.avg_prey_growth_rate[prey_subscript] = give_params_avg_prey_growth_rate(self, prey_subscript)
            self.capacity_hat[prey_subscript] = self.give_params_capacity_hat(prey_subscript)
            self.avg_prey_fitness[prey_subscript] = self.give_params_avg_prey_fitness(prey_subscript)
            self.prey_trait_response[prey_subscript] = self.give_params_prey_trait_response(prey_subscript)

        self.soln = odeint(self.f, self.y0, self.t)

        self.M = {}
        self.N = {}
        self.m = {}
        self.n = {}

        for i in xrange(0, self.num_preds):
            index = i
            self.M[str(i+1)] = self.soln[:, index]
        for i in xrange(0, self.num_preys):
            index = i + self.num_preds
            self.N[str(i+1)] = self.soln[:, index]
        for i in xrange(0, self.num_preds):
            index = i + self.num_preds + self.num_preys
            self.m[str(i+1)] = self.soln[:, index]
        for i in xrange(0, self.num_preys):
            index = i + self.num_preds + self.num_preys + self.num_preds
            self.n[str(i+1)] = self.soln[:, index]

    def f(self, y, t):
        M = [0] * self.num_preds
        for i in xrange(0, self.num_preds):
            index = i
            M[i] = y[index]

        N = [0] * self.num_preys
        for i in xrange(0, self.num_preys):
            index = i + self.num_preds
            N[i] = y[index]

        m = [0] * self.num_preds
        for i in xrange(0, self.num_preds):
            index = i + self.num_preds + self.num_preys
            m[i] = y[index]

        n = [0] * self.num_preys
        for i in xrange(0, self.num_preys):
            index = i + self.num_preds + self.num_preys + self.num_preds
            n[i] = y[index]

        f = [0] * 2*(self.num_preds + self.num_preys)
        for i in xrange(0, self.num_preds):
            index = i
            f[index] = M[i]*self.avg_pred_fitness[str(i+1)](m[i], N, n)
        for i in xrange(0, self.num_preys):
            index = i + self.num_preds
            f[index] = N[i]*self.avg_prey_fitness[str(i+1)](M, m, N[i], n[i])
        for i in xrange(0, self.num_preds):
            index = i + self.num_preds + self.num_preys
            f[index] = (self.sigmaG[str(i+1)]**2)*self.pred_trait_response[str(i+1)](N, m[i], n)
        for i in xrange(0, self.num_preys):
            index = i + self.num_preds + self.num_preys + self.num_preds
            f[index] = (self.betaG[str(i+1)]**2)*self.prey_trait_response[str(i+1)](M, m, N[i], n[i])

        return f

    def give_params_capacity_hat(self, prey_subscript):
        kappa = self.kappa[prey_subscript]
        xi = self.xi[prey_subscript]
        delta = self.delta[prey_subscript]
        beta = self.beta[prey_subscript]
        gamma = self.gamma[prey_subscript]
        phi = self.phi[prey_subscript]

        D = self.D[prey_subscript]
        E = self.E[prey_subscript]

        def capacity_hat(n):
            numerator = kappa*sqrt(E)
            denominator = delta
            exponent_num = -((beta*(phi - xi))**2 + (gamma*(n - xi))**2 + (E - (delta**2))*((n - phi)**2))
            exponent_denom = 2*D
            return (numerator/denominator)*exp(exponent_num/exponent_denom)

        return capacity_hat

    def give_params_avg_prey_fitness(self, prey_subscript):
        r = self.avg_prey_growth_rate[prey_subscript]
        K_hat = self.capacity_hat[prey_subscript]

        def avg_prey_fitness(M, m, N, n):
            fitness_sink = 0
            for pred_subscript in self.M0:
                interaction_subscript = pred_subscript + prey_subscript
                avgattack = self.avgattack[interaction_subscript]
                fitness_sink += avgattack(m[int(pred_subscript)-1], n)*M[int(pred_subscript)-1]
            return r(n)*(1 - (N/K_hat(n))) - fitness_sink

        return avg_prey_fitness

    def give_params_prey_trait_response(self, prey_subscript):
        r = self.avg_prey_growth_rate[prey_subscript]
        K_hat = self.capacity_hat[prey_subscript]

        xi = self.xi[prey_subscript]
        delta = self.delta[prey_subscript]
        phi = self.phi[prey_subscript]
        gamma = self.gamma[prey_subscript]

        B = self.B[prey_subscript]
        D = self.D[prey_subscript]

        def prey_trait_response(M, m, N, n):
            response = r(n)*(((phi - n)/(B)) - ((N*(n*(gamma**2 - delta**2) + ((delta**2)*phi - (gamma**2)*xi)))/(K_hat(n)*D)))
            for pred_subscript in self.M0:
                interaction_subscript = pred_subscript + prey_subscript
                avgattack = self.avgattack[interaction_subscript]
                A = self.A[interaction_subscript]
                theta = self.theta[interaction_subscript]
                response += avgattack(m[int(pred_subscript)-1], n)*(M[int(pred_subscript)-1]*(theta + n - m[int(pred_subscript)-1]))/(A)
            return response

        return prey_trait_response


class VariableGrowthSystem:
    def __init__(self, parameters, tf):
        self.tf = tf
        self.t = np.linspace(0, self.tf, 100000)

        self.M0 = parameters["M0"]
        self.m0 = parameters["m0"]
        self.sigma = parameters["sigma"]
        self.sigmaG = parameters["sigmaG"]
        self.d = parameters["d"]
        self.num_preds = len(self.M0)

        self.N0 = parameters["N0"]
        self.n0 = parameters["n0"]
        self.beta = parameters["beta"]
        self.betaG = parameters["betaG"]
        self.num_preys = len(self.N0)

        self.eff = parameters["eff"]
        self.tau = parameters["tau"]
        self.alpha = parameters["alpha"]
        self.theta = parameters["theta"]

        self.rho = parameters["rho"]
        self.phi = parameters["phi"]
        self.gamma = parameters["gamma"]
        self.K = parameters["K"]

        self.y0 = []
        for i in xrange(0, self.num_preds):
            self.y0.append(self.M0[str(i+1)])
        for i in xrange(0, self.num_preys):
            self.y0.append(self.N0[str(i+1)])
        for i in xrange(0, self.num_preds):
            self.y0.append(self.m0[str(i+1)])
        for i in xrange(0, self.num_preys):
            self.y0.append(self.n0[str(i+1)])

        self.A = {}
        for pred_subscript in self.M0:
            for prey_subscript in self.N0:
                interaction_subscript = pred_subscript + prey_subscript
                self.A[interaction_subscript] = self.sigma[pred_subscript]**2 + self.beta[prey_subscript]**2 + self.tau[interaction_subscript]**2

        self.B = {}
        for prey_subscript in self.N0:
            self.B[prey_subscript] = self.beta[prey_subscript]**2 + self.gamma[prey_subscript]**2

        self.avgattack = {}
        for pred_subscript in self.M0:
            for prey_subscript in self.N0:
                interaction_subscript = pred_subscript + prey_subscript
                self.avgattack[interaction_subscript] = give_params_avgattack(self, interaction_subscript)

        self.avg_pred_fitness = {}
        self.pred_trait_response = {}
        for pred_subscript in self.M0:
            self.avg_pred_fitness[pred_subscript] = give_params_avg_pred_fitness(self, pred_subscript)
            self.pred_trait_response[pred_subscript] = give_params_pred_trait_response(self, pred_subscript)

        self.avg_prey_fitness = {}
        self.prey_trait_response = {}
        self.avg_prey_growth_rate = {}
        for prey_subscript in self.N0:
            self.avg_prey_growth_rate[prey_subscript] = give_params_avg_prey_growth_rate(self, prey_subscript)
            self.avg_prey_fitness[prey_subscript] = self.give_params_avg_prey_fitness(prey_subscript)
            self.prey_trait_response[prey_subscript] = self.give_params_prey_trait_response(prey_subscript)

        self.soln = odeint(self.f, self.y0, self.t)

        self.M = {}
        self.N = {}
        self.m = {}
        self.n = {}

        for i in xrange(0, self.num_preds):
            index = i
            self.M[str(i+1)] = self.soln[:, index]
        for i in xrange(0, self.num_preys):
            index = i + self.num_preds
            self.N[str(i+1)] = self.soln[:, index]
        for i in xrange(0, self.num_preds):
            index = i + self.num_preds + self.num_preys
            self.m[str(i+1)] = self.soln[:, index]
        for i in xrange(0, self.num_preys):
            index = i + self.num_preds + self.num_preys + self.num_preds
            self.n[str(i+1)] = self.soln[:, index]

    def f(self, y, t):
        M = [0] * self.num_preds
        for i in xrange(0, self.num_preds):
            index = i
            M[i] = y[index]

        N = [0] * self.num_preys
        for i in xrange(0, self.num_preys):
            index = i + self.num_preds
            N[i] = y[index]

        m = [0] * self.num_preds
        for i in xrange(0, self.num_preds):
            index = i + self.num_preds + self.num_preys
            m[i] = y[index]

        n = [0] * self.num_preys
        for i in xrange(0, self.num_preys):
            index = i + self.num_preds + self.num_preys + self.num_preds
            n[i] = y[index]

        f = [0] * 2*(self.num_preds + self.num_preys)
        for i in xrange(0, self.num_preds):
            index = i
            f[index] = M[i]*self.avg_pred_fitness[str(i+1)](m[i], N, n)
        for i in xrange(0, self.num_preys):
            index = i + self.num_preds
            f[index] = N[i]*self.avg_prey_fitness[str(i+1)](M, m, N[i], n[i])
        for i in xrange(0, self.num_preds):
            index = i + self.num_preds + self.num_preys
            f[index] = (self.sigmaG[str(i+1)]**2)*self.pred_trait_response[str(i+1)](N, m[i], n)
        for i in xrange(0, self.num_preys):
            index = i + self.num_preds + self.num_preys + self.num_preds
            f[index] = (self.betaG[str(i+1)]**2)*self.prey_trait_response[str(i+1)](M, m, N[i], n[i])

        return f

    def give_params_avg_prey_fitness(self, prey_subscript):
        r = self.avg_prey_growth_rate[prey_subscript]
        K = self.K[prey_subscript]

        def avg_prey_fitness(M, m, N, n):
            fitness_sink = 0
            for pred_subscript in self.M0:
                interaction_subscript = pred_subscript + prey_subscript
                avgattack = self.avgattack[interaction_subscript]
                fitness_sink += avgattack(m[int(pred_subscript)-1], n)*M[int(pred_subscript)-1]
            return r(n)*(1 - (N/K)) - fitness_sink

        return avg_prey_fitness

    def give_params_prey_trait_response(self, prey_subscript):
        r = self.avg_prey_growth_rate[prey_subscript]
        phi = self.phi[prey_subscript]
        B = self.B[prey_subscript]
        K = self.K[prey_subscript]

        def prey_trait_response(M, m, N, n):
            response = ((phi - n)/B)*(1 - N/K)*r(n)
            for pred_subscript in self.M0:
                interaction_subscript = pred_subscript + prey_subscript
                avgattack = self.avgattack[interaction_subscript]
                A = self.A[interaction_subscript]
                theta = self.theta[interaction_subscript]
                response += avgattack(m[int(pred_subscript)-1], n)*(M[int(pred_subscript)-1]*(theta + n - m[int(pred_subscript)-1]))/(A)
            return response

        return prey_trait_response


class VariableCapacitySystem:
    def __init__(self, parameters, tf):
        self.tf = tf
        self.t = np.linspace(0, self.tf, 100000)

        self.M0 = parameters["M0"]
        self.m0 = parameters["m0"]
        self.sigma = parameters["sigma"]
        self.sigmaG = parameters["sigmaG"]
        self.d = parameters["d"]
        self.num_preds = len(self.M0)

        self.N0 = parameters["N0"]
        self.n0 = parameters["n0"]
        self.beta = parameters["beta"]
        self.betaG = parameters["betaG"]
        self.num_preys = len(self.N0)

        self.eff = parameters["eff"]
        self.tau = parameters["tau"]
        self.alpha = parameters["alpha"]
        self.theta = parameters["theta"]

        self.r = parameters["r"]
        self.kappa = parameters["kappa"]
        self.xi = parameters["xi"]
        self.delta = parameters["delta"]

        self.y0 = []
        for i in xrange(0, self.num_preds):
            self.y0.append(self.M0[str(i+1)])
        for i in xrange(0, self.num_preys):
            self.y0.append(self.N0[str(i+1)])
        for i in xrange(0, self.num_preds):
            self.y0.append(self.m0[str(i+1)])
        for i in xrange(0, self.num_preys):
            self.y0.append(self.n0[str(i+1)])

        self.A = {}
        for pred_subscript in self.M0:
            for prey_subscript in self.N0:
                interaction_subscript = pred_subscript + prey_subscript
                self.A[interaction_subscript] = self.sigma[pred_subscript]**2 + self.beta[prey_subscript]**2 + self.tau[interaction_subscript]**2

        self.C = {}
        for prey_subscript in self.N0:
            self.C[prey_subscript] = self.delta[prey_subscript]**2 - self.beta[prey_subscript]**2

        self.avgattack = {}
        for pred_subscript in self.M0:
            for prey_subscript in self.N0:
                interaction_subscript = pred_subscript + prey_subscript
                self.avgattack[interaction_subscript] = give_params_avgattack(self, interaction_subscript)

        self.avg_pred_fitness = {}
        self.pred_trait_response = {}
        for pred_subscript in self.M0:
            self.avg_pred_fitness[pred_subscript] = give_params_avg_pred_fitness(self, pred_subscript)
            self.pred_trait_response[pred_subscript] = give_params_pred_trait_response(self, pred_subscript)

        self.avg_prey_fitness = {}
        self.prey_trait_response = {}
        self.capacity_tilde = {}
        for prey_subscript in self.N0:
            self.capacity_tilde[prey_subscript] = self.give_params_capacity_tilde(prey_subscript)
            self.avg_prey_fitness[prey_subscript] = self.give_params_avg_prey_fitness(prey_subscript)
            self.prey_trait_response[prey_subscript] = self.give_params_prey_trait_response(prey_subscript)

        self.soln = odeint(self.f, self.y0, self.t)

        self.M = {}
        self.N = {}
        self.m = {}
        self.n = {}

        for i in xrange(0, self.num_preds):
            index = i
            self.M[str(i+1)] = self.soln[:, index]
        for i in xrange(0, self.num_preys):
            index = i + self.num_preds
            self.N[str(i+1)] = self.soln[:, index]
        for i in xrange(0, self.num_preds):
            index = i + self.num_preds + self.num_preys
            self.m[str(i+1)] = self.soln[:, index]
        for i in xrange(0, self.num_preys):
            index = i + self.num_preds + self.num_preys + self.num_preds
            self.n[str(i+1)] = self.soln[:, index]

    def f(self, y, t):
        M = [0] * self.num_preds
        for i in xrange(0, self.num_preds):
            index = i
            M[i] = y[index]

        N = [0] * self.num_preys
        for i in xrange(0, self.num_preys):
            index = i + self.num_preds
            N[i] = y[index]

        m = [0] * self.num_preds
        for i in xrange(0, self.num_preds):
            index = i + self.num_preds + self.num_preys
            m[i] = y[index]

        n = [0] * self.num_preys
        for i in xrange(0, self.num_preys):
            index = i + self.num_preds + self.num_preys + self.num_preds
            n[i] = y[index]

        f = [0] * 2*(self.num_preds + self.num_preys)
        for i in xrange(0, self.num_preds):
            index = i
            f[index] = M[i]*self.avg_pred_fitness[str(i+1)](m[i], N, n)
        for i in xrange(0, self.num_preys):
            index = i + self.num_preds
            f[index] = N[i]*self.avg_prey_fitness[str(i+1)](M, m, N[i], n[i])
        for i in xrange(0, self.num_preds):
            index = i + self.num_preds + self.num_preys
            f[index] = (self.sigmaG[str(i+1)]**2)*self.pred_trait_response[str(i+1)](N, m[i], n)
        for i in xrange(0, self.num_preys):
            index = i + self.num_preds + self.num_preys + self.num_preds
            f[index] = (self.betaG[str(i+1)]**2)*self.prey_trait_response[str(i+1)](M, m, N[i], n[i])

        return f

    def give_params_capacity_tilde(self, prey_subscript):
        kappa = self.kappa[prey_subscript]
        xi = self.xi[prey_subscript]
        delta = self.delta[prey_subscript]
        C = self.C[prey_subscript]

        def capacity_tilde(n):
            numerator = kappa*sqrt(C)
            denominator = delta
            exponent_num = -(n - xi)**2
            exponent_denom = 2*C
            return (numerator/denominator)*exp(exponent_num/exponent_denom)

        return capacity_tilde

    def give_params_avg_prey_fitness(self, prey_subscript):
        r = self.r[prey_subscript]
        K_tilde = self.capacity_tilde[prey_subscript]

        def avg_prey_fitness(M, m, N, n):
            fitness_sink = 0
            for pred_subscript in self.M0:
                interaction_subscript = pred_subscript + prey_subscript
                avgattack = self.avgattack[interaction_subscript]
                fitness_sink += avgattack(m[int(pred_subscript)-1], n)*M[int(pred_subscript)-1]
            return r*(1 - (N/K_tilde(n))) - fitness_sink

        return avg_prey_fitness

    def give_params_prey_trait_response(self, prey_subscript):
        r = self.r[prey_subscript]
        xi = self.xi[prey_subscript]
        C = self.C[prey_subscript]
        K_tilde = self.capacity_tilde[prey_subscript]

        def prey_trait_response(M, m, N, n):
            response = (r*N*(xi - n))/(K_tilde(n)*C)
            for pred_subscript in self.M0:
                interaction_subscript = pred_subscript + prey_subscript
                avgattack = self.avgattack[interaction_subscript]
                A = self.A[interaction_subscript]
                theta = self.theta[interaction_subscript]
                response += avgattack(m[int(pred_subscript)-1], n)*(M[int(pred_subscript)-1]*(theta + n - m[int(pred_subscript)-1]))/(A)
            return response

        return prey_trait_response


class NullSystem:
    def __init__(self, parameters, tf):
        self.tf = tf
        self.t = np.linspace(0, self.tf, 100000)

        self.M0 = parameters["M0"]
        self.m0 = parameters["m0"]
        self.sigma = parameters["sigma"]
        self.sigmaG = parameters["sigmaG"]
        self.d = parameters["d"]
        self.num_preds = len(self.M0)

        self.N0 = parameters["N0"]
        self.n0 = parameters["n0"]
        self.beta = parameters["beta"]
        self.betaG = parameters["betaG"]
        self.num_preys = len(self.N0)

        self.eff = parameters["eff"]
        self.tau = parameters["tau"]
        self.alpha = parameters["alpha"]
        self.theta = parameters["theta"]

        self.r = parameters["r"]
        self.K = parameters["K"]

        self.y0 = []
        for i in xrange(0, self.num_preds):
            self.y0.append(self.M0[str(i+1)])
        for i in xrange(0, self.num_preys):
            self.y0.append(self.N0[str(i+1)])
        for i in xrange(0, self.num_preds):
            self.y0.append(self.m0[str(i+1)])
        for i in xrange(0, self.num_preys):
            self.y0.append(self.n0[str(i+1)])

        self.A = {}
        for pred_subscript in self.M0:
            for prey_subscript in self.N0:
                interaction_subscript = pred_subscript + prey_subscript
                self.A[interaction_subscript] = self.sigma[pred_subscript]**2 + self.beta[prey_subscript]**2 + self.tau[interaction_subscript]**2

        self.avgattack = {}
        for pred_subscript in self.M0:
            for prey_subscript in self.N0:
                interaction_subscript = pred_subscript + prey_subscript
                self.avgattack[interaction_subscript] = give_params_avgattack(self, interaction_subscript)

        self.avg_pred_fitness = {}
        self.pred_trait_response = {}
        for pred_subscript in self.M0:
            self.avg_pred_fitness[pred_subscript] = give_params_avg_pred_fitness(self, pred_subscript)
            self.pred_trait_response[pred_subscript] = give_params_pred_trait_response(self, pred_subscript)

        self.avg_prey_fitness = {}
        self.prey_trait_response = {}
        for prey_subscript in self.N0:
            self.avg_prey_fitness[prey_subscript] = self.give_params_avg_prey_fitness(prey_subscript)
            self.prey_trait_response[prey_subscript] = self.give_params_prey_trait_response(prey_subscript)

        self.soln = odeint(self.f, self.y0, self.t)

        self.M = {}
        self.N = {}
        self.m = {}
        self.n = {}

        for i in xrange(0, self.num_preds):
            index = i
            self.M[str(i+1)] = self.soln[:, index]
        for i in xrange(0, self.num_preys):
            index = i + self.num_preds
            self.N[str(i+1)] = self.soln[:, index]
        for i in xrange(0, self.num_preds):
            index = i + self.num_preds + self.num_preys
            self.m[str(i+1)] = self.soln[:, index]
        for i in xrange(0, self.num_preys):
            index = i + self.num_preds + self.num_preys + self.num_preds
            self.n[str(i+1)] = self.soln[:, index]

    def f(self, y, t):
        M = [0] * self.num_preds
        for i in xrange(0, self.num_preds):
            index = i
            M[i] = y[index]

        N = [0] * self.num_preys
        for i in xrange(0, self.num_preys):
            index = i + self.num_preds
            N[i] = y[index]

        m = [0] * self.num_preds
        for i in xrange(0, self.num_preds):
            index = i + self.num_preds + self.num_preys
            m[i] = y[index]

        n = [0] * self.num_preys
        for i in xrange(0, self.num_preys):
            index = i + self.num_preds + self.num_preys + self.num_preds
            n[i] = y[index]

        f = [0] * 2*(self.num_preds + self.num_preys)
        for i in xrange(0, self.num_preds):
            index = i
            f[index] = M[i]*self.avg_pred_fitness[str(i+1)](m[i], N, n)
        for i in xrange(0, self.num_preys):
            index = i + self.num_preds
            f[index] = N[i]*self.avg_prey_fitness[str(i+1)](M, m, N[i], n[i])
        for i in xrange(0, self.num_preds):
            index = i + self.num_preds + self.num_preys
            f[index] = (self.sigmaG[str(i+1)]**2)*self.pred_trait_response[str(i+1)](N, m[i], n)
        for i in xrange(0, self.num_preys):
            index = i + self.num_preds + self.num_preys + self.num_preds
            f[index] = (self.betaG[str(i+1)]**2)*self.prey_trait_response[str(i+1)](M, m, N[i], n[i])

        return f

    def give_params_avg_prey_fitness(self, prey_subscript):
        r = self.r[prey_subscript]
        K = self.K[prey_subscript]

        def avg_prey_fitness(M, m, N, n):
            fitness_sink = 0
            for pred_subscript in self.M0:
                interaction_subscript = pred_subscript + prey_subscript
                avgattack = self.avgattack[interaction_subscript]
                fitness_sink += avgattack(m[int(pred_subscript)-1], n)*M[int(pred_subscript)-1]
            return r*(1 - (N/K)) - fitness_sink

        return avg_prey_fitness

    def give_params_prey_trait_response(self, prey_subscript):

        def prey_trait_response(M, m, N, n):
            response = 0
            for pred_subscript in self.M0:
                interaction_subscript = pred_subscript + prey_subscript
                avgattack = self.avgattack[interaction_subscript]
                A = self.A[interaction_subscript]
                theta = self.theta[interaction_subscript]
                response += avgattack(m[int(pred_subscript)-1], n)*(M[int(pred_subscript)-1]*(theta + n - m[int(pred_subscript)-1]))/(A)
            return response

        return prey_trait_response


def get_system_dimension(set_):
    number_of_predators = str(len(set_["predator"]["initial_values"]["densities"]))
    number_of_preys = str(len(set_["prey"]["initial_values"]["densities"]))
    dim = "%sx%s" % (number_of_predators, number_of_preys)
    return dim


def get_system_model(set_):
    if "max_growth_rates" in set_["prey"]:
        growth = "variable"
    elif "intrinsic_growth_rates" in set_["prey"]:
        growth = "constant"
    else:
        raise ValueError
    if "max_carrying_capacities" in set_["prey"]:
        capacity = "variable"
    elif "carrying_capacities" in set_["prey"]:
        capacity = "constant"
    else:
        raise ValueError

    return (growth, capacity)


def SET_TRAIT_GRAPH_LIMITS(args):
    global LOWER_LIMIT
    global UPPER_LIMIT

    default_limit = 5

    if args.trait_graph_lower_limit == 0.0:
        args.trait_graph_lower_limit = -0.000001
    if args.trait_graph_upper_limit == 0.0:
        args.trait_graph_upper_limit = 0.000001

    if not args.trait_graph_lower_limit and not args.trait_graph_upper_limit:
        LOWER_LIMIT = -default_limit
        UPPER_LIMIT = default_limit
    elif args.trait_graph_lower_limit and not args.trait_graph_upper_limit:
        LOWER_LIMIT = args.trait_graph_lower_limit
        if LOWER_LIMIT < default_limit:
            UPPER_LIMIT = default_limit
        else:
            UPPER_LIMIT = 2*default_limit + LOWER_LIMIT
    elif not args.trait_graph_lower_limit and args.trait_graph_upper_limit:
        UPPER_LIMIT = args.trait_graph_upper_limit
        if UPPER_LIMIT > -default_limit:
            LOWER_LIMIT = -default_limit
        else:
            LOWER_LIMIT = UPPER_LIMIT - 2*default_limit
    else:
        if args.trait_graph_upper_limit > args.trait_graph_lower_limit:
            LOWER_LIMIT = args.trait_graph_lower_limit
            UPPER_LIMIT = args.trait_graph_upper_limit
        else:
            print "Lower limit must be less than Upper limit ;)"
            raise ValueError


def save_OUT_and_ERR(OUT, ERR, current_directory):
    total = header("STD_OUT") + OUT + "\n"*10 + header("STD_ERR")
    if len(ERR) == 0:
        total += "\n<<< NONE >>>\n\n"
    else:
        total += ERR
    filename = "%s/stability_results.txt" % current_directory
    with open(filename, "w") as stab_res_file_object:
        stab_res_file_object.write(total)


def PARSE_ARGS():
    parser = argparse.ArgumentParser()
    parser.add_argument("config")
    parser.add_argument("-f", "--final-values", action="store_true", dest="final_values", default=False, help=F_OPTION_HELP)
    parser.add_argument("-k", "--keep-orignial-images", action="store_true", dest="keep_original_images", default=False, help=K_OPTION_HELP)
    parser.add_argument("-n", "--no-combine", action="store_false", dest="combine", default=True, help=N_OPTION_HELP)
    parser.add_argument("-p", "--no-parameters", action="store_false", dest="display_parameters", default=True, help=P_OPTION_HELP)
    parser.add_argument("-L", "--no-legend", action="store_false", dest="legend", default=True, help=L_OPTION_HELP)
    parser.add_argument("-l", "--lower-limit", dest="trait_graph_lower_limit", type=float, help=LOWER_LIMIT_HELP)
    parser.add_argument("-u", "--upper-limit", dest="trait_graph_upper_limit", type=float, help=UPPER_LIMIT_HELP)
    parser.add_argument("-o", "--opt_traits_horiz_line", action="store_true", dest="opt_traits_horiz_line", default=False, help=O_OPTION_HELP)
    return parser.parse_args()


def plot_time_graphs(current_directory, step, system, text, args, time, dimension):
    ts = default_timer()
    densities_file = "%s/densities_%03d.png" % (current_directory, step)
    plot_densities(system, densities_file, text, args)
    data_time = default_timer() - ts
    print "Time taken for this data: %.03f\n" % (data_time)
    time += data_time

    ts = default_timer()
    traits_file = "%s/traits_%03d.png" % (current_directory, step)
    plot_traits(system, traits_file, text, args)
    data_time = default_timer() - ts
    print "Time taken for this data: %.03f\n" % (data_time)
    time += data_time

    combined_time_graphs = "%s/combined_time_graphs_%.3d.png" % (current_directory, step)
    if args.combine:
        combine_images(densities_file, traits_file, combined_time_graphs, args.keep_original_images, IMAGEMAGICK_SIDE_BY_SIDE_COMMAND)


def plot_phase_planes(current_directory, step, system, text, args, time, dimension):
    if dimension == "1x1":
        (density_phase_plane_plot, trait_phase_plane_plot) = (plot_1x1_densities_phase_plane, plot_1x1_traits_phase_plane)
    elif dimension == "1x2":
        (density_phase_plane_plot, trait_phase_plane_plot) = (plot_1x2_densities_phase_plane, plot_1x2_traits_phase_plane)
    elif dimension == "2x1":
        (density_phase_plane_plot, trait_phase_plane_plot) = (plot_2x1_densities_phase_plane, plot_2x1_traits_phase_plane)

    ts = default_timer()
    density_phase_plane_file = "%s/density_phase_plane_%03d.png" % (current_directory, step)
    density_phase_plane_plot(system, density_phase_plane_file, text, args)
    data_time = default_timer() - ts
    print "Time taken for this data: %.03f\n" % (data_time)

    ts = default_timer()
    trait_phase_plane_file = "%s/trait_phase_plane_%03d.png" % (current_directory, step)
    trait_phase_plane_plot(system, trait_phase_plane_file, text, args)
    data_time = default_timer() - ts
    print "Time taken for this data: %.03f\n" % (data_time)

    combined_phase_planes = "%s/combined_phase_planes_%.3d.png" % (current_directory, step)
    if args.combine:
        combine_images(density_phase_plane_file, trait_phase_plane_file, combined_phase_planes, args.keep_original_images, IMAGEMAGICK_SIDE_BY_SIDE_COMMAND)


def plot_graphs(current_directory, step, system, text, args, time, dimension):
    plot_time_graphs(current_directory, step, system, text, args, time, dimension)
    if dimension in SUPPORTED_DIMENSIONS:
        plot_phase_planes(current_directory, step, system, text, args, time, dimension)


def main():
    args = PARSE_ARGS()

    SET_TRAIT_GRAPH_LIMITS(args)

    data = json.loads(open(args.config).read())

    for set_ in data["system_parameters"]:
        now = datetime.now()
        date_time_stamp = now.strftime(DATE_TIME_DIRECTORY_FORMAT)
        dimension = get_system_dimension(set_)
        (growth, capacity) = get_system_model(set_)

        if growth == "variable" and capacity == "variable":
            system_type = "full_system"
        elif growth == "variable" and capacity == "constant":
            system_type = "variable_growth"
        elif growth == "constant" and capacity == "variable":
            system_type = "variable_capacity"
        elif growth == "constant" and capacity == "constant":
            system_type = "null_system"

        current_directory = "%s/graphs/%s/%s/%s" % (DIRECTORY, system_type, dimension, date_time_stamp)
        os.system("mkdir -p %s" % current_directory)

        relevant_data_file = "%s/relevant_data.json" % current_directory
        pretty_data = json.dumps(set_, indent=4, sort_keys=True)
        with open(relevant_data_file, "w") as rel_data_file_object:
            rel_data_file_object.write(pretty_data)

        stability_check = "%s/bin/stability_checks/%s_check.py" % (DIRECTORY, dimension)
        if os.path.isfile(stability_check):
            command_line = [stability_check, relevant_data_file]
            p = Popen(command_line, stdout=PIPE, stderr=PIPE)
            (OUT, ERR) = p.communicate()
            save_OUT_and_ERR(OUT, ERR, current_directory)
        else:
            print(NO_STABILITY_CHECK)

        config_descriptions_to_variables = {
            "M0":     set_["predator"]["initial_values"]["densities"],
            "m0":     set_["predator"]["initial_values"]["traits"],
            "sigma":  set_["predator"]["trait_variances"]["total"],
            "sigmaG": set_["predator"]["trait_variances"]["genetic"],
            "d":      set_["predator"]["death_rates"],

            "N0":     set_["prey"]["initial_values"]["densities"],
            "n0":     set_["prey"]["initial_values"]["traits"],
            "beta":   set_["prey"]["trait_variances"]["total"],
            "betaG":  set_["prey"]["trait_variances"]["genetic"],

            "eff":    set_["interaction_parameters"]["efficiencies"],
            "tau":    set_["interaction_parameters"]["specialization"],
            "alpha":  set_["interaction_parameters"]["max_attack_rates"],
            "theta":  set_["interaction_parameters"]["optimal_trait_differences"],
        }
        if growth == "variable":
            config_descriptions_to_variables["rho"] = set_["prey"]["max_growth_rates"]
            config_descriptions_to_variables["phi"] = set_["prey"]["optimum_growth_rate_trait_values"]
            config_descriptions_to_variables["gamma"] = set_["prey"]["growth_rate_stabilizing_selection_variances"]
        elif growth == "constant":
            config_descriptions_to_variables["r"] = set_["prey"]["intrinsic_growth_rates"]
        if capacity == "variable":
            config_descriptions_to_variables["kappa"] = set_["prey"]["max_carrying_capacities"]
            config_descriptions_to_variables["xi"] = set_["prey"]["optimum_carrying_capacity_trait_values"]
            config_descriptions_to_variables["delta"] = set_["prey"]["carrying_capacity_stabilizing_selection_variances"]
        elif capacity == "constant":
            config_descriptions_to_variables["K"] = set_["prey"]["carrying_capacities"]

        # Parameter Step
        steps = int(set_["steps"])

        def get_step(dictionary):
            return (dictionary["stop"] - dictionary["start"])/steps

        number_of_graphs = (steps+1)*4
        time_needed = AVG_TIME_PER_GRAPH*number_of_graphs/60.
        print NUMBER_OF_GRAPHS_MESSAGE % number_of_graphs
        print TIME_NEEDED_MESSAGE % time_needed

        starts_and_steps = {}

        for variable, dict_ in config_descriptions_to_variables.iteritems():
            starts_and_steps[variable] = {}
            for subscript, parameter_range in dict_.iteritems():
                starts_and_steps[variable][subscript] = {}
                step = get_step(parameter_range)
                starts_and_steps[variable][subscript]["start"] = parameter_range["start"]
                starts_and_steps[variable][subscript]["step"] = step

        time = 0
        # Final time
        tf = set_["final_time"]

        for step in irange(0, steps, 1):

            parameters = {}
            for variable, dict_ in starts_and_steps.iteritems():
                parameters[variable] = {}
                for subscript, parameter_range in dict_.iteritems():
                    parameters[variable][subscript] = parameter_range["start"] + (step*parameter_range["step"])

            if system_type == "full_system":
                system = FullSystem(parameters, tf)
            elif system_type == "variable_growth":
                system = VariableGrowthSystem(parameters, tf)
            elif system_type == "variable_capacity":
                system = VariableCapacitySystem(parameters, tf)
            elif system_type == "null_system":
                system = NullSystem(parameters, tf)

            ### Get parameters in text format for the graphs
            text = []
            for variable in parameters:
                for subscript, value in parameters[variable].iteritems():
                    if starts_and_steps[variable][subscript]["step"] != 0:
                        text.append(r"$%s%s}= %.03f$" % (LaTeX_VARIABLE_FORMAT[variable], subscript, value))

            plot_graphs(current_directory, step, system, text, args, time, dimension)

        ### This separate section is so that the images are guaranteed to exist before attempting to combine them
        print "\n\nUsing ImageMagick to combine graphs.....\n\n"
        for step in irange(0, steps, 1):
            combined_time_graphs = "%s/combined_time_graphs_%.3d.png" % (current_directory, step)
            combined_phase_planes = "%s/combined_phase_planes_%.3d.png" % (current_directory, step)
            output = "%s/output_%.3d.png" % (current_directory, step)
            if args.combine and dimension in SUPPORTED_DIMENSIONS:
                combine_images(combined_time_graphs, combined_phase_planes, output, args.keep_original_images, IMAGEMAGICK_ON_TOP_COMMAND)

        print date_time_stamp
        print "\a\a\a"

if __name__ == "__main__":
    main()
