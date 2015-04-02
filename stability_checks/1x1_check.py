from __future__ import division
from math import sqrt
import json
import argparse
import os

DIRECTORY = os.path.dirname(os.path.abspath(__file__))

def irange(start, stop, step):
    r = start
    while r <= stop:
        yield r
        r += step

def get_N_star(d, sigma, beta, tau, e, alpha):
    A = sigma**2 + beta**2 + tau**2
    num = d*sqrt(A)
    denom = e*alpha*tau
    return num/denom

def check_condition(p):
    eff    = p["eff"]
    alpha  = p["alpha"]
    theta  = p["theta"]
    tau    = p["tau"]
    d      = p["d"]
    sigma  = p["sigma"]
    sigmaG = p["sigmaG"]
    K      = p["K"]
    gamma  = p["gamma"]
    rho    = p["rho"]
    phi    = p["phi"]
    beta   = p["beta"]
    betaG  = p["betaG"]

    A = sigma**2 + beta**2 + tau**2
    B = beta**2 + gamma**2

    r_hat = (rho*gamma)/(sqrt(B))
    N_star = get_N_star(d, sigma, beta, tau, eff, alpha)

    LHS = d*(sigmaG**2)
    RHS = r_hat*(betaG**2)*(1 - (N_star/K))*(1 - A/B)

    analysis = "unstable"
    if LHS > RHS:
        sign = ">"
        analysis = "stable"
    elif LHS < RHS:
        sign = "<"
    else:
        sign = "="

    
    print "%.07f %s %.07f\n%s\n" % (LHS, sign, RHS, analysis)

def PARSE_ARGS():
    parser = argparse.ArgumentParser()
    parser.add_argument("relavent_data")
    return parser.parse_args()

def main():
    args = PARSE_ARGS()

    data = json.loads(open("%s/%s" % (DIRECTORY, args.relavent_data)).read())

    config_descriptions_to_variables = {
        "M0"    : data["predator"]["initial_values"]["densities"],
        "m0"    : data["predator"]["initial_values"]["traits"],
        "sigma" : data["predator"]["trait_variances"]["total"],
        "sigmaG": data["predator"]["trait_variances"]["genetic"],
        "d"     : data["predator"]["death_rates"],

        "N0"    : data["prey"]["initial_values"]["densities"],
        "n0"    : data["prey"]["initial_values"]["traits"],
        "beta"  : data["prey"]["trait_variances"]["total"],
        "betaG" : data["prey"]["trait_variances"]["genetic"],
        "K"     : data["prey"]["carrying_capacities"],
        "rho"   : data["prey"]["max_growth_rates"],
        "phi"   : data["prey"]["optimum_trait_values"],
        "gamma" : data["prey"]["cost_variances"],

        "eff"   : data["interaction_parameters"]["efficiencies"],
        "tau"   : data["interaction_parameters"]["specialization"],
        "alpha" : data["interaction_parameters"]["max_attack_rates"],
        "theta" : data["interaction_parameters"]["optimal_trait_differences"],
    }
    steps = int(data["steps"])
    def get_step(dictionary):
        return (dictionary["stop"] - dictionary["start"])/steps

    starts_and_steps = {}
    for variable, dict_ in config_descriptions_to_variables.iteritems():    
        starts_and_steps[variable] = {}
        for subscript, parameter_range in dict_.iteritems():
            starts_and_steps[variable][subscript] = {}
            step = get_step(parameter_range)
            starts_and_steps[variable][subscript]["start"] = parameter_range["start"]
            starts_and_steps[variable][subscript]["step"] = step

    tf = data["final_time"]
    print steps
    for step in irange(0, steps, 1):
        parameters = {}
        for variable, dict_ in starts_and_steps.iteritems():
            for subscript, parameter_range in dict_.iteritems():
                parameters[variable] = parameter_range["start"] + (step*parameter_range["step"])
        print "graph %.03d" % step
        check_condition(parameters)


if __name__ == "__main__":
    main()