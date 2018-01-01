import sys
import json
import numpy as np
import multiprocessing
import pickle
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import seaborn

def main(sargs):
    recompute = True

    if recompute:
        outfile = sargs[0]
        with open(outfile,'r') as f:
            data = json.load(f)
        data_list = data['epdm']

        ncores = multiprocessing.cpu_count()
        print('pooling with %d cores' % ncores)
        pool = multiprocessing.Pool(ncores)
        parsed = pool.map(parse_one_elem,data_list)
        pool.close()

        times,full_species_map = collate_all_elements(parsed)
        species_statistics = get_species_statistics(times,full_species_map)
        cache_these = dict(times=times,full_species_map=full_species_map,species_statistics=species_statistics)
        make_cache(cache_these)

    else:
        cache = load_cache()
        times,full_species_map,species_statisics = cache['times'],cache['full_species_map'],cache['species_statistics']

    make_plots(species_statistics,full_species_map)

def make_cache(cache_this,fname='cache.pkl'):
    with open(fname,'wb') as f:
        pickle.dump(cache_this,f)

def load_cache(fname='cache.pkl'):
    with open(fname,'rb') as f:
        cache = pickle.load(f)
    return cache

def collate_all_elements(parsed):
    nsteps = len(parsed)
    times = np.zeros(nsteps)

    full_species_map = {}
    for idx,p in enumerate(parsed):
        times[idx] = p[0]
        for species in p[1].keys():
            if species not in full_species_map.keys():
                full_species_map[species] = {}
                full_species_map[species]['num'] = np.zeros(nsteps)
                full_species_map[species]['frac'] = np.zeros(nsteps)

            full_species_map[species]['num'][idx] = p[1][species]['num']
            full_species_map[species]['frac'][idx] = p[1][species]['frac']

        if idx % 1000 == 0:
            print('packed record %d/%d' % (idx,nsteps))

    return times,full_species_map

def get_species_statistics(times,full_species_map):
    species_statistics = {}

    for k,v in full_species_map.items():
        inv_time = (1.0/times[-1])
        mean_value = inv_time*np.trapz(v['num'],x=times)
        mean_fraction = inv_time*np.trapz(v['frac'],x=times)
        species_statistics[k] = (mean_value,mean_fraction)

    return species_statistics

def make_plots(species_statistics,full_species_map):
    population_vs_nresidues = np.zeros(max([len(species) for species in full_species_map.keys()]))
    for species,(pop,frac) in species_statistics.items():
        population_vs_nresidues[len(species)-1] += frac

    max_length = population_vs_nresidues.size
    plt.figure()
    plt.plot(np.arange(1,max_length+1),np.log(population_vs_nresidues),'o-')
    plt.xlabel(r'$n$')
    plt.ylabel(r'$\log\,\,\mathrm{Pr}(N_{\mathrm{residues}} = n)$')
    plt.savefig('flory_distribution.png')

def parse_one_elem(elem):
    species_map = {}
    time = elem['current_time']

    for species,pop_num in elem['populations'].items():
        species_map[species] = {}
        species_map[species]['num'] = pop_num

    tot_num_species = np.sum([species_map[species]['num'] for species in species_map.keys()])

    for species in elem['populations'].keys():
        species_map[species]['frac'] = float(species_map[species]['num'])/float(tot_num_species)

    return time,species_map

if __name__ == '__main__':
    main(sys.argv[1:])
