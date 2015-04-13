#!/usr/bin/python

from __future__ import division
from math import sqrt
import gc as garbage
import matplotlib.pyplot as plt
import pylab
import json
from datetime import datetime


def make_stability_check(data):
    data["d"]
    data["sigma"]
    data["beta"]
    data["ratio"]
    data["tau"]
    data["gamma"]
    data["rho"]
    data["e"]
    data["alpha"]
    data["K"]


    # A = sigma**2 + beta**2 + tau**2
    # B = beta**2 + gamma**2
    # r = (rho*gamma)/sqrt(B)
    # N_star = (d*sqrt(A))/(e*alpha*tau)

    if len(data["ratio"]) == 1:
        LHS = lambda x: data["ratio"][0]**2
    else:
        LHS = lambda x: x**2

        if len(data["beta"]) == 1:
            if len(data):
                pass

    ## NNOOOOOOOOO





    def stability_check(var_1, var_2):
        pass



    return stability_check



def main():
    NOW = datetime.now()
    date_time_stamp = NOW.strftime(DATE_TIME_DIRECTORY_FORMAT)

    data = json.loads(open("data.json").read())

    stability_check = make_stability_check(data)



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
