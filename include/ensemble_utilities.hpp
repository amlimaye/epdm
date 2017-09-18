#include "types.hpp"
#include <json/json.h>

namespace ensemble_utilities {
	//initialize an empty ensemble
	ensemble_t initialize_ensemble(const uint32_t seed);

	//add a new population to ensemble
	void add_population_to_ensemble(ensemble_t* ensemble, population_t population);

	//sample a reaction according to a weighted distribution over reaction propensities
	const reaction_t sample_reaction(const ensemble_t& ensemble, double rand);

	//check if population is already tabulated
	std::tuple<bool,size_t> 
		is_in_tabulated_populations(const ensemble_t& ensemble, const species_t& species);

	//take a timestep
	void take_timestep(ensemble_t& ensemble);

	//function to print out ensemble to std::cout
	void print_ensemble(const ensemble_t& in);

	//function to serialize ensemble to a Json object
	Json::Value serialize_to_json(const ensemble_t& in);
}