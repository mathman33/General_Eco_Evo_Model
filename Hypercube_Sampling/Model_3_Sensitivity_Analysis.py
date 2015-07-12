from __future__ import division
import numpy as np
import math
import matplotlib.pyplot as plt
import os
# import random
import scipy.stats
import pylab
from time import sleep

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


def f_3(ratio, d, kappa, delta, sigma, beta, alpha, tau, e, r):
    A = sigma**2 + beta**2 + tau**2
    C = delta**2 - beta**2
    max_K = (kappa*math.sqrt(C))/(delta)
    max_a = (alpha*tau)/(math.sqrt(A))
    N = (d)/(e*max_a)

    return ratio**2 - (r/d)*(1 - (N/max_K)*(1 + A/C))


def main():
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
        "kappa": (1, 500)
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
        function_values[i] = f_3(**kwargs)

    print "Function Values Obtained"

    statistical_values = {}
    for key, values in hyper_cube_sample.iteritems():
        spearman, p_value = scipy.stats.spearmanr(values, function_values)

        statistical_values[key] = {
            "SPEARMAN": spearman,
            "P_VALUE": p_value
        }

    print "Statistical Values Computed"

    for key, values in hyper_cube_sample.iteritems():
        plt.figure()
        plt.scatter(values, function_values, color="brown")
        plt.ylabel(r"$f_3$", fontsize=20, rotation=0)
        plt.xlabel(PLAIN_TO_LATEX[key], fontsize=20)
        plt.ylim(-20, 50)
        plt.axhline(0)
        plt.title(PLAIN_TO_LATEX[key] + " --- SPEARMAN: %.6f" % statistical_values[key]["SPEARMAN"] + " --- P-VALUE: %.6f" % statistical_values[key]["P_VALUE"])
        plt.savefig(key+".png", format="png", dpi=200)
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
    pylab.savefig("bar_plot.png", format="png", dpi=200)
    pylab.close()

    print "Bar Graph Built"

    plt.figure()
    plt.axes(frameon=False)
    frame = pylab.gca()
    frame.axes.get_xaxis().set_ticks([])
    frame.axes.get_yaxis().set_ticks([])
    plt.text(0.25, 0.5, "MODEL 3", fontsize=50)
    plt.savefig("TITLE_3.png", format="png", dpi=200)
    plt.close()

    print "Title Page Built"

    os.system("convert ratio.png d.png kappa.png +append output1.png")
    os.system("convert delta.png sigma.png beta.png +append output2.png")
    os.system("convert alpha.png tau.png e.png +append output3.png")
    os.system("convert bar_plot.png r.png TITLE_3.png -append output4.png")
    os.system("convert output1.png output2.png output3.png -append output5.png")
    os.system("convert output5.png output4.png +append Sensitivity_Analysis_Model_3.png")
    sleep(1)
    for param in parameter_list + ["output1", "output2", "output3", "output4", "output5", "bar_plot", "TITLE_3"]:
        os.system("rm %s.png" % param)

    print "Graphs Merged"
    print "\n\nCOMPLETE\n\n"

if __name__ == "__main__":
    main()
