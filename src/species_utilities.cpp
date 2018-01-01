#include "species_utilities.hpp"

species_t species_utilities::make_void_species() {
    species_t void_species;
    void_species.name = "void";
    void_species.folded = false;
    void_species.foldable = false;
    return void_species;
}

species_t species_utilities::make_h_species() {
    species_t h_species;
    h_species.name = "H";
    h_species.folded = false;
    h_species.foldable = false;
    return h_species;
}

species_t species_utilities::make_p_species() {
    species_t p_species;
    p_species.name = "P";
    p_species.folded = false;
    p_species.foldable = false;
    return p_species;
}	

species_t species_utilities::make_arbitrary_species(std::string name) {
	species_t new_species;
	new_species.name = name;
    new_species.folded = false;
    new_species.foldable = false;
	return new_species;
}
