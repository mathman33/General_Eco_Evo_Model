#!/usr/bin/python

from __future__ import division
import gc as garbage
import matplotlib.pyplot as plt
import numpy as np
from pylab import meshgrid, cm
import json
from datetime import datetime
import os
import argparse


DATE_TIME_DIRECTORY_FORMAT = '%y%m%d_%H%M%S'
LaTeX_VARIABLE_FORMAT = {
    "sigma": "\\sigma",
    "d": "d",
    "beta": "\\beta",
    "kappa": "\\kappa",
    "e": "e",
    "tau": "\\tau",
    "alpha": "\\alpha",
    "r": "r",
    "delta": "\\delta",
    "ratio": "\\sigma_G/\\beta_G"
}


def make_stability_checks(data):

    def coexistence_stability_check(x, y):
        for key, value in data.iteritems():
            if type(value) == list:
                if value[0] == "x":
                    data[key] = x
                elif value[0] == "y":
                    data[key] = y

        A = data["sigma"]**2 + data["beta"]**2 + data["tau"]**2
        C = data["delta"]**2 - data["beta"]**2
        K_tilde = (data["kappa"]*np.sqrt(C))/(data["delta"])
        N_star = (data["d"]*np.sqrt(A))/(data["e"]*data["alpha"]*data["tau"])

        LHS = data["ratio"]**2
        RHS = (data["r"]/data["d"])*(1 - ((N_star/K_tilde)*(1 + (A/C))))

        return (LHS - RHS)

    def exclusion_stability_check(x, y):
        for key, value in data.iteritems():
            if type(value) == list:
                if value[0] == "x":
                    data[key] = x
                elif value[0] == "y":
                    data[key] = y

        A = data["sigma"]**2 + data["beta"]**2 + data["tau"]**2
        C = data["delta"]**2 - data["beta"]**2
        K_tilde = (data["kappa"]*np.sqrt(C))/(data["delta"])

        LHS = data["d"]
        RHS = (data["e"]*data["alpha"]*data["tau"]*K_tilde)/np.sqrt(A)

        return (LHS - RHS)

    return (coexistence_stability_check, exclusion_stability_check)


def make_title(data):
    Title = ""
    for key, value in data.iteritems():
        if type(value) == float:
            Title += r"$%s = %.3f$, " % (LaTeX_VARIABLE_FORMAT[str(key)], value)
    return Title


def PARSE_ARGS():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", dest="contour_line_delta", type=float, default=0.5)
    parser.add_argument("-r", dest="contour_range", type=float, nargs=2, default=[-1., 1.])
    parser.add_argument("-p", dest="print_parameters", action="store_true", default=False)
    return parser.parse_args()


def main():
    args = PARSE_ARGS()

    NOW = datetime.now()
    date_time_stamp = NOW.strftime(DATE_TIME_DIRECTORY_FORMAT)

    cwd = os.path.dirname(os.path.realpath(__file__))
    data_file = "variable_capacity_1x1_data.json"
    data = json.loads(open(os.path.join(cwd, data_file)).read())

    (coexistence_stability_check, exclusion_stability_check) = make_stability_checks(data)

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
    Z = {
        "coexistence": coexistence_stability_check(X, Y),
        "exclusion": exclusion_stability_check(X, Y)
    }

    plots = "plots"
    plots_dir = os.path.join(cwd, plots)
    if not os.path.isdir(plots_dir):
        os.system("mkdir %s" % plots)

    for type_ in ["coexistence", "exclusion"]:
        direc = os.path.join(plots_dir, type_, "variable_capacity")
        if not os.path.isdir(direc):
            os.system("mkdir %s" % direc)

        plt.figure()
        plt.xlabel(r"$%s$" % str(LaTeX_VARIABLE_FORMAT[x_var]), fontsize=15, rotation=0)
        plt.ylabel(r"$%s$" % str(LaTeX_VARIABLE_FORMAT[y_var]), fontsize=15, rotation=0)
        try:
            l = len(Z[type_])
            warning = False
        except:
            warning = True
            constant_value = float(Z[type_])
            Z[type_] = np.asarray([[float(Z[type_])]*mesh_refinement]*mesh_refinement)
        CS = plt.contour(X, Y, Z[type_], np.arange(min(args.contour_range), max(args.contour_range), args.contour_line_delta), cmap=cm.RdBu)
        plt.clabel(CS, inline=1, fontsize=10)

        if args.print_parameters:
            Title = make_title(data)
            plt.title(Title, fontsize=15)

        file_name = "%s_%s_%s.png" % (x_var, y_var, date_time_stamp)
        file_path = os.path.join(plots, type_, "variable_capacity", file_name)
        save_location = os.path.join(direc, file_name)
        try:
            plt.savefig(save_location, format="png")
            print "\nFile saved: %s" % file_path
            if warning:
                print "Warning: %s function is constant with respect to %s and %s" % (type_, x_var, y_var)
                print "  Value: %.3f" % constant_value
            else:
                pass
        except:
            print "Unable to save figure"
            plt.show()

        plt.close()
        garbage.collect()

    print

if __name__ == "__main__":
    main()
