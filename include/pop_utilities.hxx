#include "types.hxx"

namespace pop_utilities {
	population_t make_void_population();
	void update_molecule_count(population_t* pop, long int delta);
	void print_population(const population_t& in);
} 
