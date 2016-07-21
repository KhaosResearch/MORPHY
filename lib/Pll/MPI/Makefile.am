CC = @MPICC@

lib_LTLIBRARIES =
lib_LIBRARIES =

if BUILD_MPI 

include $(top_srcdir)/sources.am

if BUILD_SSE3
lib_LTLIBRARIES += ../src/libpll-sse3-mpi.la
lib_LIBRARIES += ../src/libpll-sse3-mpi.a
___src_libpll_sse3_mpi_la_SOURCES = $(ALL_SOURCES)
___src_libpll_sse3_mpi_la_CFLAGS = -D_FINE_GRAIN_MPI -D__SSE3 -msse3 -c -O2 -fomit-frame-pointer -funroll-loops -Wall -Wredundant-decls  -Wreturn-type  -Wswitch-default -Wimplicit  -Wimplicit-function-declaration  -Wimplicit-int -Wimport -Wunused-label -Wno-int-to-pointer-cast -Wbad-function-cast  -Wmissing-declarations -Wmissing-prototypes  -Wnested-externs  -Wold-style-definition -Wstrict-prototypes -Wpointer-sign -Wextra -Wredundant-decls -Wunused -Wunused-function -Wunused-parameter -Wunused-value  -Wunused-variable -Wformat  -Wformat-nonliteral -Wparentheses -Wsequence-point -Wuninitialized -Wundef 
___src_libpll_sse3_mpi_la_LDFLAGS = -version-info @LIBPLL_MAJOR@:@LIBPLL_MINOR@:@LIBPLL_REV@
___src_libpll_sse3_mpi_a_SOURCES = $(___src_libpll_sse3_mpi_la_SOURCES)
___src_libpll_sse3_mpi_a_CFLAGS = $(___src_libpll_sse3_mpi_la_CFLAGS)
endif

if BUILD_AVX
lib_LTLIBRARIES += ../src/libpll-avx-mpi.la
lib_LIBRARIES += ../src/libpll-avx-mpi.a
___src_libpll_avx_mpi_la_SOURCES = $(AVX_SOURCES)
___src_libpll_avx_mpi_la_CFLAGS = -D_FINE_GRAIN_MPI -D__SSE3 -D__AVX -mavx -c -O2 -fomit-frame-pointer -funroll-loops -Wall -Wredundant-decls  -Wreturn-type  -Wswitch-default -Wimplicit  -Wimplicit-function-declaration  -Wimplicit-int -Wimport -Wunused-label -Wno-int-to-pointer-cast -Wbad-function-cast  -Wmissing-declarations -Wmissing-prototypes  -Wnested-externs  -Wold-style-definition -Wstrict-prototypes -Wpointer-sign -Wextra -Wredundant-decls -Wunused -Wunused-function -Wunused-parameter -Wunused-value  -Wunused-variable -Wformat  -Wformat-nonliteral -Wparentheses -Wsequence-point -Wuninitialized -Wundef 
___src_libpll_avx_mpi_la_LDFLAGS = -version-info @LIBPLL_MAJOR@:@LIBPLL_MINOR@:@LIBPLL_REV@
___src_libpll_avx_mpi_a_SOURCES = $(___src_libpll_avx_mpi_la_SOURCES)
___src_libpll_avx_mpi_a_CFLAGS = $(___src_libpll_avx_mpi_la_CFLAGS)
endif

#pkgincludedir=$(includedir)/pll
#pkginclude_HEADERS = $(top_buildir)/src/pll.h $(top_buildir)/src/newick.h $(top_buildir)/src/stack.h $(top_buildir)/src/hash.h $(top_buildir)/src/errcodes.h $(top_buildir)/src/globalVariables.h $(top_buildir)/src/lexer.h $(top_buildir)/src/genericParallelization.h $(top_buildir)/src/treeIO.h $(top_buildir)/src/queue.h $(top_buildir)/src/parsePartition.h $(top_buildir)/src/mem_alloc.h $(top_buildir)/src/cycle.h
#
endif
