CXX = g++
CXXFLAGS = -std=c++20 -O2
INCLUDES = -I../universal/include
LDFLAGS = -lm

all: ladder cg_compare static_conditioning cg_compare_seeded

static_conditioning: src/static_conditioning.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LDFLAGS)

ladder: src/ladder.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LDFLAGS)

cg_compare: src/cg_compare.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LDFLAGS)

clean:
	rm -f ladder cg_compare cg_compare_seeded static_conditioning

cg_compare_cond_probe: src/cg_compare_cond_probe.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LDFLAGS)

cg_compare_seeded: src/cg_compare_seeded.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LDFLAGS)
