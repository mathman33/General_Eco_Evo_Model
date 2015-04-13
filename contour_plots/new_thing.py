#!/usr/bin/python

from __future__ import division
from math import sqrt
import gc as garbage
import matplotlib.pyplot as plt
import json
from datetime import datetime


DATE_TIME_DIRECTORY_FORMAT = '%y%m%d_%H%M%S'


def coexistence_stability(p):
    A = p["sigma"]**2 + p["beta"]**2 + p["tau"]**2
    B = p["beta"]**2 + p["gamma"]**2
    r = (p["rho"]*p["gamma"])/sqrt(B)
    N_star = (p["d"]*sqrt(A))/(p["e"]*p["alpha"]*p["tau"])

    LHS = ((p["ratio"])**2)
    RHS = (r/p["d"])*(1 - N_star/p["K"])*(1 - (A/B))

    stable = LHS > RHS
    return stable


def exclusion_stability(p):
    A = p["sigma"]**2 + p["beta"]**2 + p["tau"]**2
    LHS = p["d"]
    RHS = (p["e"]*p["alpha"]*p["tau"]*p["K"])/(sqrt(A))

    stable = LHS > RHS
    return stable


def irange(start, stop, step):
    r = start
    while r <= stop:
        yield r
        r += step


def main():
    NOW = datetime.now()
    date_time_stamp = NOW.strftime(DATE_TIME_DIRECTORY_FORMAT)
    data = json.loads(open("data.json").read())

    changing_vars = []
    constant_vars = []
    for key, value in data.iteritems():
        if len(value) == 2:
            changing_vars.append(key)
        elif len(value) == 1:
            constant_vars.append(key)

    if len(changing_vars) != 2:
        raise ValueError
    else:
        var_1 = changing_vars[0]
        var_2 = changing_vars[1]

    p = {}
    for var in constant_vars:
        p[var] = data[var][0]

    plt.figure()
    for i in irange(data[var_1][0], data[var_1][1], ((data[var_1][1] - data[var_1][0])/100)):
        p[var_1] = i
        for j in irange(data[var_2][0], data[var_2][1], ((data[var_2][1] - data[var_2][0])/100)):
            p[var_2] = j

            plt.xlabel(var_1)
            plt.ylabel(var_2)

            if exclusion_stability(p):
                plt.plot(p[var_1], p[var_2], 'rs')
            elif coexistence_stability(p):
                plt.plot(p[var_1], p[var_2], 'gs')
            else:
                plt.plot(p[var_1], p[var_2], 'ys')

    file_name = "%s_%s_%s" % (var_1, var_2, date_time_stamp)
    plt.savefig(file_name, format="png")
    plt.close()
    garbage.collect()


if __name__ == "__main__":
    main()
