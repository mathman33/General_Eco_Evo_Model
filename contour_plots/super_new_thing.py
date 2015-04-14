#!/usr/bin/python

from __future__ import division
import gc as garbage
import matplotlib.pyplot as plt
import numpy as np
from pylab import meshgrid, cm
import json
from datetime import datetime


DATE_TIME_DIRECTORY_FORMAT = '%y%m%d_%H%M%S'
LaTeX_VARIABLE_FORMAT = {
    "M0"    : "M_0",
    "m0"    : "m_0",
    "sigma" : "\\sigma",
    "sigmaG": "\\sigma_G",
    "d"     : "d",

    "N0"    : "N_0",
    "n0"    : "n_0",
    "beta"  : "\\beta",
    "betaG" : "\\beta_G",
    "r"     : "r",
    "K"     : "K",

    "e"     : "e",
    "tau"   : "\\tau",
    "alpha" : "\\alpha",
    "theta" : "\\theta",

    "rho"   : "\\rho",
    "phi"   : "\\phi",
    "gamma" : "\\gamma",

    "ratio" : "\\frac{\\sigma_G}{\\beta_G}"
}


def make_stability_checks(p):

    def coex_stability_check(x, y):
        for k, v in p.iteritems():
            if v == "var_1":
                p[k] = x
            elif v == "var_2":
                p[k] = y

        A = p["sigma"]**2 + p["beta"]**2 + p["tau"]**2
        B = p["beta"]**2 + p["gamma"]**2
        r = (p["rho"]*p["gamma"])/np.sqrt(B)
        N_star = (p["d"]*np.sqrt(A))/(p["e"]*p["alpha"]*p["tau"])

        LHS = p["ratio"]**2
        RHS = (r/p["d"])*(1 - N_star/p["K"])*(1 - (A/B))

        return (LHS - RHS)

    def excl_stability_check(x, y):
        for k, v in p.iteritems():
            if v == "var_1":
                p[k] = x
            elif v == "var_2":
                p[k] = y

        A = p["sigma"]**2 + p["beta"]**2 + p["tau"]**2
        LHS = p["d"]**2
        RHS = (p["e"]*p["alpha"]*p["tau"]*p["K"])/np.sqrt(A)

        return (LHS - RHS)

    return (coex_stability_check, excl_stability_check)


def make_title(data, changing_vars):
    Title = ""
    count = 0
    for k, v in data.iteritems():
        if k not in changing_vars:
            Title += r"$%s = %.1f$, " % (LaTeX_VARIABLE_FORMAT[str(k)], v[0])
            if count == 4:
                Title += "\n"
            count += 1
    return Title


def main():
    NOW = datetime.now()
    date_time_stamp = NOW.strftime(DATE_TIME_DIRECTORY_FORMAT)

    data = json.loads(open("data.json").read())

    p = {}
    changing_vars = []
    for k, v in data.iteritems():
        if len(v) == 1:
            p[k] = v[0]
        else:
            changing_vars.append(k)
            p[k] = "var_%d" % len(changing_vars)
    if len(changing_vars) != 2:
        raise ValueError

    (coex_stability_check, excl_stability_check) = make_stability_checks(p)

    x = np.arange(data[changing_vars[0]][0], data[changing_vars[0]][1], ((data[changing_vars[0]][1] - data[changing_vars[0]][0])/100))
    y = np.arange(data[changing_vars[1]][0], data[changing_vars[1]][1], ((data[changing_vars[1]][1] - data[changing_vars[1]][0])/100))
    X, Y = meshgrid(x, y)
    Z = coex_stability_check(X, Y)

    plt.figure()
    plt.xlabel(r"$%s$" % str(LaTeX_VARIABLE_FORMAT[changing_vars[0]]), fontsize=20, rotation=0)
    plt.ylabel(r"$%s$" % str(LaTeX_VARIABLE_FORMAT[changing_vars[1]]), fontsize=20, rotation=0)
    CS = plt.contour(X, Y, Z, np.arange(-2, 2, 0.5), cmap=cm.RdBu)
    plt.clabel(CS, inline=1, fontsize=10)

    Title = make_title(data, changing_vars)
    plt.title(Title)

    file_name = "%s_%s_%s" % (changing_vars[0], changing_vars[1], date_time_stamp)
    plt.savefig(file_name, format="png")
    print file_name
    plt.close()
    garbage.collect()

if __name__ == "__main__":
    main()
