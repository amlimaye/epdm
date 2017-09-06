#!/usr/local/bin/python

from __future__ import print_function
import sys
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def main(sargs):
	outfile = sargs[0]
	with open(outfile,'r') as f:
		data = json.load(f)

	data_list = data['epdm']
	times = np.array([elem['current_time'] for elem in data_list])

	species = [[pop['species_name'] for pop in elem['populations']] for elem in data_list]
	distinct_species = set(sum(species,[])) - set(['void'])
	print(distinct_species)

	species_map = {s:[] for s in distinct_species}

	for elem in data_list:
		pops = {pop['species_name']:pop['num_molecules'] for pop in elem['populations']}
		for species in species_map.keys():
			species_map[species].append(pops.get(species,0))

	#the average number of species must be computed with a _time_ average!
	for k,v in species_map.items():
		mean_value = (1.0/times[-1])*np.trapz(v,x=times)
		print('<n_{%s}(t)> = %0.4f' % (k,mean_value))
		ax = sns.tsplot(data=v,time=times)
		ax.set_xlabel(r'$t$')
		ax.set_ylabel(r'$n(t)$')

		plt.show()

if __name__ == '__main__':
	main(sys.argv[1:])