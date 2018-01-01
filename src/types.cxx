#include "types.hxx"

bool operator==(const relation_address_t& x, const relation_address_t& y) {
    return ((x.relation_ptr->owner_population_ptr->species.name == y.relation_ptr->owner_population_ptr->species.name) &&
            (x.relation_ptr->owner_population_ptr->species.folded == y.relation_ptr->owner_population_ptr->species.folded));
}

bool operator==(const relation_t& x, const relation_t& y) {
    return ((x.relation_address_ptr->owner_population_ptr->species.name == y.relation_address_ptr->owner_population_ptr->species.name) &&
            (x.relation_address_ptr->owner_population_ptr->species.folded == y.relation_address_ptr->owner_population_ptr->species.folded));
}
