#!/usr/bin/python

from __future__ import division
import gc as garbage
import matplotlib.pyplot as plt
import numpy as np
from pylab import meshgrid, cm
import json
from datetime import datetime
import os


DATE_TIME_DIRECTORY_FORMAT = '%y%m%d_%H%M%S'
LaTeX_VARIABLE_FORMAT = {
    "sigma" : "\\sigma",
    "d"     : "d",
    "beta"  : "\\beta",
    "K"     : "K",
    "e"     : "e",
    "tau"   : "\\tau",
    "alpha" : "\\alpha",
    "rho"   : "\\rho",
    "gamma" : "\\gamma",
    "ratio" : "\\frac{\\sigma_G}{\\beta_G}"
}


def make_stability_checks(data):
    dictionary = data

    def coex_stability_check(x, y):
        for key, value in dictionary.iteritems():
            if type(value) == list:
                if value[0] == "x":
                    dictionary[key] = x
                elif value[0] == "y":
                    dictionary[key] = y

        A = dictionary["sigma"]**2 + dictionary["beta"]**2 + dictionary["tau"]**2
        B = dictionary["beta"]**2 + dictionary["gamma"]**2
        r = (dictionary["rho"]*dictionary["gamma"])/np.sqrt(B)
        N_star = (dictionary["d"]*np.sqrt(A))/(dictionary["e"]*dictionary["alpha"]*dictionary["tau"])

        LHS = dictionary["ratio"]**2
        RHS = (r/dictionary["d"])*(1 - N_star/dictionary["K"])*(1 - (A/B))

        return (LHS - RHS)

    def excl_stability_check(x, y):
        for key, value in dictionary.iteritems():
            if type(value) == list:
                if value[0] == "x":
                    dictionary[key] = x
                elif value[0] == "y":
                    dictionary[key] = y

        A = dictionary["sigma"]**2 + dictionary["beta"]**2 + dictionary["tau"]**2
        LHS = dictionary["d"]**2
        RHS = (dictionary["e"]*dictionary["alpha"]*dictionary["tau"]*dictionary["K"])/np.sqrt(A)

        return (LHS - RHS)

    return (coex_stability_check, excl_stability_check)


def make_title(data):
    Title = ""
    for key, value in data.iteritems():
        if type(value) == float:
            Title += r"$%s = %.1f$, " % (LaTeX_VARIABLE_FORMAT[str(key)], value)
    return Title


def main():
    NOW = datetime.now()
    date_time_stamp = NOW.strftime(DATE_TIME_DIRECTORY_FORMAT)

    data = json.loads(open("data.json").read())

    (coex_stability_check, excl_stability_check) = make_stability_checks(data)

    mesh_refinement = 100
    for key, value in data.iteritems():
        if type(value) == list:
            if value[0] == "x":
                x_var = key
                x_start = value[1][0]
                x_stop = value[1][1]
                x_delta = (x_stop - x_start)/mesh_refinement
            elif value[0] == "y":
                y_var = key
                y_start = value[1][0]
                y_stop = value[1][1]
                y_delta = (y_stop - y_start)/mesh_refinement

    x = np.arange(x_start, x_stop, x_delta)
    y = np.arange(y_start, y_stop, y_delta)
    X, Y = meshgrid(x, y)
    Z = coex_stability_check(X, Y)

    plt.figure()
    plt.xlabel(r"$%s$" % str(LaTeX_VARIABLE_FORMAT[x_var]), fontsize=15, rotation=0)
    plt.ylabel(r"$%s$" % str(LaTeX_VARIABLE_FORMAT[y_var]), fontsize=15, rotation=0)
    CS = plt.contour(X, Y, Z, np.arange(-2, 2, 0.5), cmap=cm.RdBu)
    plt.clabel(CS, inline=1, fontsize=10)

    Title = make_title(data)
    plt.title(Title, fontsize=20)

    plots = "plots"
    cwd = os.path.dirname(os.path.realpath(__file__))
    if not os.path.isdir(os.path.join(cwd, plots)):
        os.system("mkdir %s" % plots)

    file_name = os.path.join(plots, "%s_%s_%s" % (x_var, y_var, date_time_stamp))
    plt.savefig(file_name, format="png")

    print "\nFile saved: %s\n" % file_name

    plt.close()
    garbage.collect()

if __name__ == "__main__":
    main()
