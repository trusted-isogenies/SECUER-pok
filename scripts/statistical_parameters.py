from math import *

sec_levels = [128, 128, 192, 256]
e = [(216, 137), (250, 159), (305, 192), (372, 239)]


def stat_distance_two(p, k):
	return -2 + log(p-1, 2)/2 + log(k + 1/3, 2) - k/2 

def stat_distance_three(p, two_len, k):
	x = -2 + log(p-1, 2)/2 + log(1 + 2**(two_len/2)*(1+1/3)**(1/2), 2)
	sec = x + log(k + 2/4, 2) - k/2*log(3, 2) 
	return sec

def find_k_two(p, security=128):
	for k in range(100, 2000):
		if stat_distance_two(p, k) < -security:
			return k

def find_k_three(p, two_len, security=128):
	for k in range(100, 2000):
		if stat_distance_three(p, two_len, k) < -security:
			return k


scaling_constant = 1/log(3/2, 2)

print("log(p)\treps\ttwo-length\tthree-length\tw\th\n")



for i, (m, n) in enumerate(e):
	p = 2**m * 3**n - 1
	logp = ceil(log(p, 2))


	two_len = logp
	three_len = find_k_three(p, two_len, sec_levels[i])
	reps = ceil(sec_levels[i]*scaling_constant)
	w = ceil(two_len / m)
	h = ceil(three_len / n)

	print(" {:} \t {:} \t  {:4} \t\t   {:4} \t{:} \t{:}".format(logp, reps, two_len, three_len, w, h))


print()

for i, (m, n) in enumerate(e):
	p = 2**m * 3**n - 1
	logp = ceil(log(p, 2))


	two_len = find_k_two(p, sec_levels[i])
	three_len = find_k_three(p, two_len, sec_levels[i])
	reps = ceil(sec_levels[i]*scaling_constant)
	w = ceil(two_len / m)
	h = ceil(three_len / n)

	print(" {:} \t {:} \t  {:4} \t\t   {:4} \t{:} \t{:}".format(logp, reps, two_len, three_len, w, h))