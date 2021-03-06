#!/usr/bin/python

from __future__ import division
import numpy as np
import math
import matplotlib.pyplot as plt
import argparse
import os
# import random
import scipy.stats
import pylab
from time import sleep

RATIO_OPTION_HELP = """
Use division and a comparison against 1, rather than subtraction and a comparison against 0.
"""
KEEP_OPTION_HELP = """
Keep all individual images.  Do not remove them after combining.
"""

VERTICAL_RATIO = {
    "ratio": r"$\frac{\sigma_G}{\beta_G}$",
    "d": r"$d$",
    "kappa": r"$\kappa$",
    "delta": r"$\delta$",
    "sigma": r"$\sigma$",
    "beta": r"$\beta$",
    "alpha": r"$\alpha$",
    "tau": r"$\tau$",
    "e": r"$e$",
    "r": r"$r$"
}

PLAIN_TO_LATEX = {
    "ratio": r"$\sigma_G/\beta_G$",
    "d": r"$d$",
    "kappa": r"$\kappa$",
    "delta": r"$\delta$",
    "sigma": r"$\sigma$",
    "beta": r"$\beta$",
    "alpha": r"$\alpha$",
    "tau": r"$\tau$",
    "e": r"$e$",
    "r": r"$r$"
}


def S_3(ratio, d, kappa, delta, sigma, beta, alpha, tau, e, r):
    A = sigma**2 + beta**2 + tau**2
    C = delta**2 - beta**2
    max_K = (kappa*math.sqrt(C))/(delta)
    max_a = (alpha*tau)/(math.sqrt(A))
    N = (d)/(e*max_a)

    return ratio**2 - (r/d)*(1 - (N/max_K)*(1 + A/C))


def S_3_RATIO(ratio, d, kappa, delta, sigma, beta, alpha, tau, e, r):
    A = sigma**2 + beta**2 + tau**2
    C = delta**2 - beta**2
    max_K = (kappa*math.sqrt(C))/(delta)
    max_a = (alpha*tau)/(math.sqrt(A))
    N = (d)/(e*max_a)

    return ratio**2/((r/d)*(1 - (N/max_K)*(1 + A/C)))


def PARSE_ARGS():
    parser = argparse.ArgumentParser()
    parser.add_argument("--RATIO", action="store_true", dest="RATIO", default=False, help=RATIO_OPTION_HELP)
    parser.add_argument("--KEEP", action="store_true", dest="KEEP", default=False, help=KEEP_OPTION_HELP)
    return parser.parse_args()


def main():
    args = PARSE_ARGS()

    # random_grid = np.random.randint(resolution, size=(samples, variables))
    parameter_list = ["ratio", "d", "kappa", "delta", "sigma", "beta", "alpha", "tau", "e", "r"]

    parameter_ranges = {
        "ratio": (0.25, 4),
        "d": (0.01, 1),
        "r": (0.01, 1),
        "delta": (0.1, 3),
        "sigma": (0.01, 0.1),
        "beta": (0.01, 0.1),
        "alpha": (0.01, 1),
        "tau": (0.01, 3),
        "e": (0.01, 0.5),
        "kappa": (50, 500)
    }

    parameters = len(parameter_list)
    resolution = 10000

    parameter_grid = {}
    for parameter in parameter_list:
        min_ = parameter_ranges[parameter][0]
        max_ = parameter_ranges[parameter][1]
        parameter_grid[parameter] = list(np.linspace(min_, max_, resolution, endpoint=True))

    print "Parameter Grid Formed"

    hyper_cube_sample = {}
    for parameter in parameter_list:
        sample = [0]*resolution
        for i in xrange(0, resolution):
            rand_num = np.random.randint(resolution-i)
            sample[i] = parameter_grid[parameter][rand_num]
            parameter_grid[parameter].pop(rand_num)
        hyper_cube_sample[parameter] = sample

    print "Hypercube Sample Formed"

    function_values = [0]*resolution
    for i in xrange(0, resolution):
        kwargs = {}
        for key, values in hyper_cube_sample.iteritems():
            kwargs[key] = values[i]
        if args.RATIO:
            function_values[i] = S_3_RATIO(**kwargs)
        else:
            function_values[i] = S_3(**kwargs)

    print "Function Values Obtained"

    statistical_values = {}
    for key, values in hyper_cube_sample.iteritems():
        spearman, p_value = scipy.stats.spearmanr(values, function_values)

        statistical_values[key] = {
            "SPEARMAN": spearman,
            "P_VALUE": p_value
        }

    print "Statistical Values Computed"

    directory = "Figures_Model_3"
    if os.path.isdir(directory):
        overwrite = str(raw_input("overwrite? [y/n] --> "))
        if overwrite.lower() == "y" or overwrite.lower() == "yes":
            os.system("rm -r %s" % directory)
            os.system("mkdir %s" % directory)
        else:
            return(0)
    else:
        os.system("mkdir %s" % directory)
    for key, values in hyper_cube_sample.iteritems():
        plt.figure()
        plt.scatter(values, function_values, color="brown")
        plt.ylabel(r"$S_3$", fontsize=20, rotation=0)
        plt.xlabel(PLAIN_TO_LATEX[key], fontsize=20)
        plt.ylim(-20, 50)
        if args.RATIO:
            plt.axhline(1)
        else:
            plt.axhline(0)
        plt.title(PLAIN_TO_LATEX[key] + " --- SPEARMAN: %.6f" % statistical_values[key]["SPEARMAN"] + " --- P-VALUE: %.6f" % statistical_values[key]["P_VALUE"])
        plt.savefig(directory + "/" + key+".png", format="png", dpi=200)
        plt.close()

    print "Parameter Graphs Built"

    SP_VALS = []
    P_VALS = []
    for parameter in parameter_list:
        SP_VALS.append(statistical_values[parameter]["SPEARMAN"])
        P_VALS.append(statistical_values[parameter]["P_VALUE"])

    pos = np.arange(parameters)+.5    # the bar centers on the y axis
    pylab.figure()
    pylab.barh(pos, SP_VALS, align='center', color="brown")
    LATEX_PARAMETER_LIST = [VERTICAL_RATIO[p] for p in parameter_list]
    pylab.yticks(pos, LATEX_PARAMETER_LIST, fontsize=20)
    pylab.xlabel('Spearman Constant', fontsize=20)
    pylab.ylabel("Parameter", fontsize=20)
    pylab.savefig(directory + "/bar_plot.png", format="png", dpi=200)
    pylab.close()

    print "Bar Graph Built"

    plt.figure()
    plt.axes(frameon=False)
    frame = pylab.gca()
    frame.axes.get_xaxis().set_ticks([])
    frame.axes.get_yaxis().set_ticks([])
    plt.text(0.25, 0.5, "MODEL 3", fontsize=50)
    plt.savefig(directory + "/TITLE_3.png", format="png", dpi=200)
    plt.close()

    print "Title Page Built"

    os.system("convert %s/ratio.png %s/d.png %s/kappa.png +append %s/output1.png" % tuple([directory]*4))
    os.system("convert %s/delta.png %s/sigma.png %s/beta.png +append %s/output2.png" % tuple([directory]*4))
    os.system("convert %s/alpha.png %s/tau.png %s/e.png +append %s/output3.png" % tuple([directory]*4))
    os.system("convert %s/bar_plot.png %s/r.png %s/TITLE_3.png -append %s/output4.png" % tuple([directory]*4))
    os.system("convert %s/output1.png %s/output2.png %s/output3.png -append %s/output5.png" % tuple([directory]*4))
    FILENAME = "%s/Sensitivity_Analysis_Model_3" % directory
    if args.RATIO:
        FILENAME += "_RATIO"
    FILENAME += ".png"
    os.system("convert %s/output5.png %s/output4.png +append %s" % tuple([directory]*2 + [FILENAME]))
    sleep(1)
    if not args.KEEP:
        for param in parameter_list + ["output1", "output2", "output3", "output4", "output5", "bar_plot", "TITLE_3"]:
            os.system("rm %s/%s.png" % param)

    print "Graphs Merged"
    print "\n\nCOMPLETE\n\n"

if __name__ == "__main__":
    main()
