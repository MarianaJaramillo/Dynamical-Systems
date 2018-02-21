//
// Thanks to Patrick Roberts
// https://stackoverflow.com/users/1541563/patrick-roberts
// https://stackoverflow.com/questions/48653228/in-c-how-to-write-a-n-dimensional-nested-for-loop-where-n-is-variable
//

#include <stdio.h>
#include <assert.h>

void loopn_recurse (int *min, int *max, int *counters, size_t N, size_t n, void (*func)(int*, size_t)) {
  for (int i = min[n]; i < max[n]; ++i) {
    counters[n] = i;

    if (N - n > 1) {
      loopn_recurse(min, max, counters, N, n + 1, func);
    } else {
      // innermost algorithm
      func(counters, N);
    }
  }
}

void loopn (int *min, int *max, int *counters, size_t N, void (*func)(int*, size_t)) {
  assert(N > 0);
  loopn_recurse(min, max, counters, N, 0, func);
}

// example usage

void test (int *counters, size_t N) {
  putchar('a');

  for (size_t j = 0; j < N; ++j) {
    printf("[%d] ", counters[j]);
  }

  putchar('\n');
}

int main () {
  int min[3] = {1, 2, 3};
  int max[3] = {2, 4, 6};
  int counters[3];
  size_t N=3;

  loopn(min, max, counters, N, test);
}
