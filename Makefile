SRC_DIR=src
INCL_DIR=include
BUILD_DIR=.
JSONLIBDIR=/usr/local/Cellar/jsoncpp/1.8.3/lib/

CXXFLAGS=-O2 -g -std=c++11
FOLDING_FLAG=-D__INCLUDE_FOLDING

types:
	g++ -c -o $(BUILD_DIR)/$@.o $(CXXFLAGS) $(FOLDING_FLAG) -I$(INCL_DIR)\
		$(SRC_DIR)/$@.cxx

folding:
	g++ -c -o $(BUILD_DIR)/$@.o $(CXXFLAGS) $(FOLDING_FLAG) -I$(INCL_DIR) \
		$(SRC_DIR)/$@.cxx

ensemble_utilities:
	g++ -c -o $(BUILD_DIR)/$@.o $(CXXFLAGS) $(FOLDING_FLAG) -I$(INCL_DIR) \
		$(SRC_DIR)/$@.cxx

rxn_utilities:
	g++ -c -o $(BUILD_DIR)/$@.o $(CXXFLAGS) $(FOLDING_FLAG) -I$(INCL_DIR) \
		$(SRC_DIR)/$@.cxx

pop_utilities:
	g++ -c -o $(BUILD_DIR)/$@.o $(CXXFLAGS) $(FOLDING_FLAG) -I$(INCL_DIR) \
		$(SRC_DIR)/$@.cxx

species_utilities:
	g++ -c -o $(BUILD_DIR)/$@.o $(CXXFLAGS) $(FOLDING_FLAG) -I$(INCL_DIR) \
		$(SRC_DIR)/$@.cxx

binary: types folding ensemble_utilities rxn_utilities pop_utilities species_utilities
	g++ -o main $(CXXFLAGS) $(SRC_DIR)/toplevel.cxx $(addsuffix .o, $^) \
		-I$(INCL_DIR) -L$(JSONLIBDIR) -ljsoncpp

binary-debug: types folding ensemble_utilities rxn_utilities pop_utilities species_utilities
	g++ -o main $(CXXFLAGS) $(SRC_DIR)/toplevel.cxx $(addsuffix .o, $^) \
		-I$(INCL_DIR) -L$(JSONLIBDIR) -ljsoncpp

clean:
	rm -rf $(BUILD_DIR)/*.o $(BUILD_DIR)/main $(BUILD_DIR)/main.dSYM
	rm *.pkl *.json
