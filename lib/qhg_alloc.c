#include <stdlib.h>
#include <stdio.h>
#include <qhg_defs.h>

/*
 * Malloc, with minimal error handling
 */
void *
qhg_alloc(size_t size)
{
  void *ptr;
  posix_memalign(&ptr, QHG_MEMALIGN, size);
  if(ptr == NULL) {
    fprintf(stderr, " posix_memalign(&ptr, %d, %lu) returned NULL. Out of memory?\n", QHG_MEMALIGN, size);
    exit(-1);
  }
  return ptr;
}
