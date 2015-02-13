from math import sqrt

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
	d = p["d"]
	sigma = p["sigma"]
	sigmaG = p["sigmaG"]
	beta = p["beta"]
	betaG = p["betaG"]
	K = p["K"]
	tau = p["tau"]
	e = p["e"]
	alpha = p["alpha"]
	r = p["r"]

	N_star = get_N_star(d, sigma, beta, tau, e, alpha)

	LHS = d*(sigmaG**2)
	RHS = r*(betaG**2)*(1 - (N_star/K))

	analysis = "unstable"
	if LHS > RHS:
		sign = ">"
		analysis = "stable"
	elif LHS < RHS:
		sign = "<"
	else:
		sign = "="

	
	print "%.04f %s %.04f\n%s\n" % (LHS, sign, RHS, analysis)

count = 1
for i in irange(0.10, 0.25, 0.01):
	p = {
		"d" : 0.05,
		"sigma" : 0.25,
		"sigmaG" : i,
		"beta" : 0.25,
		"betaG" : 0.1,
		"K" : 225,
		"tau" : 0.1,
		"e" : 0.5,
		"alpha" : 0.05,
		"r" : 0.2,
	}
	print "graph %.03d" % count
	check_condition(p)
	count += 1