
SIZES = 434 503 610 751

################################################################

SIDH = PQCrypto-SIDH

################################################################

TARGETS_PROVE = $(foreach t, $(SIZES), prove_$(t))
TARGETS_VERIFY = $(foreach t, $(SIZES), verify_$(t))
TARGETS_TEST = $(foreach t, $(SIZES), test_$(t))
TARGETS = $(TARGETS_PROVE) $(TARGETS_VERIFY) $(TARGETS_TEST)

.DELETE_ON_ERROR:

.PHONY: all test clean

all: $(TARGETS)

test:
	$(foreach t, $(SIZES), $(MAKE) runtest$(t); )

runtest%: test_% prove_% verify_%
	time ./test_$*
	time ./prove_$* --initial > /tmp/proof0.txt
	time ./verify_$* < /tmp/proof0.txt > /tmp/curve0.txt
	cat /tmp/curve0.txt
	time tail -n1 /tmp/curve0.txt | ./prove_$* > /tmp/proof1.txt
	time ./verify_$* < /tmp/proof1.txt > /tmp/curve1.txt
	cat /tmp/curve1.txt

clean:
	rm -f $(TARGETS)
	$(MAKE) -C $(SIDH) clean

################################################################

THREAD_CFLAGS = $(shell echo 'int main() {}' | $(CC) -lpthread -x c - -o /dev/null \
                        && printf "%s %s" -DTHREADING -pthread)

EXTRA_CFLAGS = -O2 -g -fsanitize=address -fsanitize=undefined       # sanitizer
EXTRA_CFLAGS = -O3 -g -pg -DNDEBUG                                  # profiling
EXTRA_CFLAGS = -O2 -g                                               # debug
EXTRA_CFLAGS = -O3 -DNDEBUG                                         # optimized
EXTRA_CFLAGS = -O3 -DNDEBUG $(THREAD_CFLAGS)                        # w/ threads

.PRECIOUS: $(SIDH)/%/libsidh.a
$(SIDH)/%/libsidh.a: $(shell find $(SIDH)/src/)
	$(MAKE) -C $(SIDH) EXTRA_CFLAGS="-fPIE $(EXTRA_CFLAGS)" $*

ADDITIONAL_SETTINGS=-march=native
ifeq "$(ARCH)" "M1"
	ADDITIONAL_SETTINGS=-mcpu=apple-a14
endif

RECIPE = $(CC) \
		-std=c11 -pedantic \
		-Wall -Wextra -Werror \
		$(ADDITIONAL_SETTINGS) \
		$(EXTRA_CFLAGS) \
		-I $(SIDH)/src/ \
		-D POIK_BITS=$* \
		<NAME>.c \
		$(SIDH)/lib$(*)poik/libsidh.a \
		-o <NAME>_$*

prove_%: prove.c $(SIDH)/lib%poik/libsidh.a
	$(subst <NAME>, prove, $(RECIPE))

verify_%: verify.c $(SIDH)/lib%poik/libsidh.a
	$(subst <NAME>, verify, $(RECIPE))

test_%: test.c $(SIDH)/lib%poik/libsidh.a
	$(subst <NAME>, test, $(RECIPE))
