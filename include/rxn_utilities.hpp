#include "types.hpp"

namespace rxn_utilities {
	reaction_t spawn_rxn(const species_t spawned_species);
	reaction_t decay_rxn(const species_t decay_species);
	reaction_t elongation_rxn(const species_t s1, const species_t s2);
	reaction_t splitting_rxn(const species_t s1, size_t position);
    reaction_t folding_rxn(species_t in, int num_contacts);
	reaction_t unfolding_rxn(species_t in);
	double compute_partial_propensity(const reaction_t& reaction, long int num_molecules);
	std::list<reaction_t> get_reactions_for_species_pair(const population_t& p1, const population_t& p2);
	void print_reaction(const reaction_t& in);
};