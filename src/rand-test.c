#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>
#include "rand.h"
#include "nifty.h"

int
main(int argc, char *argv[])
{
	static uint_fast32_t buckets[256U];
	const float lambda = (float)(argc > 1 ? strtod(argv[1], NULL) : 1.f);

	for (size_t i = 0; i < 1000000U; i++) {
		float x = dr_rand_poiss(lambda);

		if (x < (float)countof(buckets)) {
			uint_fast8_t xu = (uint_fast8_t)x;
			buckets[xu]++;
		}
	}
	for (size_t i = 0; i < countof(buckets); i++) {
		printf("%zu\t%" PRIuFAST32 "\n", i, buckets[i]);
	}
	return 0;
}

/* rand-test.c ends here */
