#ifndef __mem_alloc_h
#define __mem_alloc_h
#include <stddef.h>
#include <stdlib.h>
#ifdef __linux__
#include <malloc.h>
#endif
#include "pll.h"

#define rax_memalign memalign
#define rax_malloc malloc
#define rax_realloc realloc
#define rax_free free
#define rax_posix_memalign posix_memalign
#define rax_calloc calloc

#endif
