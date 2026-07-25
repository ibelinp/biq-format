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

# -B keeps Python from writing __pycache__ next to the sources. Cosmetic for a
# build, but this tree gets published by copying the folder, and a stray .pyc
# would ride along.
test: biq_test
	./biq_test golden
	@command -v python3 >/dev/null && python3 -B biq.py --test || \
	  echo "(python3 not found — skipped the Python reference check)"

# The golden vectors come from an independent Rust implementation that lives in
# the author's tree, not here. They are checked into golden/ precisely so you
# don't need it — this target exists only for regenerating them if the codec
# ever changes.
golden:
	@test -f ../../crates/spectra-dsp/Cargo.toml || { \
	  echo "make golden needs the author's Rust tree, which isn't part of this repo."; \
	  echo "The vectors in golden/ are checked in — that's what 'make test' uses."; \
	  exit 1; }
	cd ../.. && cargo run -q -p spectra-dsp --example biq_golden -- biq-format/reference/golden

clean:
	rm -f biqconv biq_test
	rm -rf __pycache__

.PHONY: all test golden clean
