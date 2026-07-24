# biq reference implementation.
#
#   make          build biqconv
#   make test     build and run the self-test (includes the Rust golden vectors)
#   make golden   regenerate golden/ from the Rust encoder (needs the SpectraWeb tree)

CXX      ?= c++
CXXFLAGS ?= -std=c++11 -O2 -Wall -Wextra

all: biqconv

biqconv: biqconv.cpp biq.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<

biq_test: biq_test.cpp biq.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<

test: biq_test
	./biq_test golden
	@command -v python3 >/dev/null && python3 biq.py --test || \
	  echo "(python3 not found — skipped the Python reference check)"

# The golden vectors come from the independent Rust implementation in
# crates/spectra-dsp. Regenerating them is only needed if the codec changes.
golden:
	cd ../.. && cargo run -q -p spectra-dsp --example biq_golden -- biq-format/reference/golden

clean:
	rm -f biqconv biq_test

.PHONY: all test golden clean
