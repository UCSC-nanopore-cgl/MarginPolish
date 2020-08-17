cmake_minimum_required(VERSION 3.5)
project(marginPhase VERSION 1.1.0)

#set(CMAKE_CXX_STANDARD 14)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -DNDEBUG -std=gnu99 -D_GNU_SOURCE -D_XOPEN_SOURCE=500 -D_POSIX_C_SOURCE=200809L -mpopcnt ")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -DNDEBUG -std=gnu99 -D_GNU_SOURCE -D_XOPEN_SOURCE=500 -D_POSIX_C_SOURCE=200809L -mpopcnt ")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O0 -fno-inline -UNDEBUG -fsanitize=address -static-libasan")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O0 -fno-inline -UNDEBUG")

project(marginPhase LANGUAGES C)

# configure a header file to pass some of the CMake settings
set(POLISH_VERSION_MAJOR 1)
set(POLISH_VERSION_MINOR 3)
execute_process(
        COMMAND git log -1 --format=%h
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
)
configure_file(
        "${PROJECT_SOURCE_DIR}/inc/marginVersion.h.in"
        "${PROJECT_SOURCE_DIR}/inc/marginVersion.h"
)

include(GNUInstallDirs)
set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY true)

###############################
##### HEADERS AND SOURCES #####
###############################

set(HTSLIB_HEADERS
        externalTools/htslib
        externalTools/htslib/htslib
        externalTools/htslib/cram
        externalTools/htslib/os
        )
set(SONLIB_HEADERS
        externalTools/sonLib/C/inc
        externalTools/sonLib/C/impl
        externalTools/sonLib/C/tests
        externalTools/sonLib/externalTools/cutest
        externalTools/sonLib/externalTools/quicktree_1.1/include
        )
set(SSW_HEADERS
        externalTools/ssw/
        )

set(LP_HEADERS
        externalTools/lp_solve/
        externalTools/lp_solve/bfp/
        externalTools/lp_solve/bfp/bfp_LUSOL/
        externalTools/lp_solve/bfp/bfp_LUSOL/LUSOL/
        externalTools/lp_solve/colamd/
        externalTools/lp_solve/shared/)


include_directories(inc)
include_directories(externalTools/samtools)
include_directories(${HTSLIB_HEADERS})
include_directories(${SONLIB_HEADERS})
include_directories(${SSW_HEADERS})
include_directories(${LP_HEADERS})


set(HTSLIB_SOURCES
        externalTools/htslib/cram/cram_codecs.c
        externalTools/htslib/cram/cram_decode.c
        externalTools/htslib/cram/cram_encode.c
        externalTools/htslib/cram/cram_external.c
        externalTools/htslib/cram/cram_index.c
        externalTools/htslib/cram/cram_io.c
        externalTools/htslib/cram/cram_samtools.c
        externalTools/htslib/cram/cram_stats.c
        externalTools/htslib/cram/files.c
        externalTools/htslib/cram/mFILE.c
        externalTools/htslib/cram/open_trace_file.c
        externalTools/htslib/cram/pooled_alloc.c
        externalTools/htslib/cram/rANS_static.c
        externalTools/htslib/cram/sam_header.c
        externalTools/htslib/cram/string_alloc.c
        externalTools/htslib/os/rand.c
        externalTools/htslib/test/fieldarith.c
        externalTools/htslib/test/hfile.c
        externalTools/htslib/test/hts_endian.c
        externalTools/htslib/test/sam.c
        externalTools/htslib/test/test-bcf-sr.c
        externalTools/htslib/test/test-regidx.c
        externalTools/htslib/test/test-vcf-api.c
        externalTools/htslib/test/test-vcf-sweep.c
        externalTools/htslib/test/test_bgzf.c
        externalTools/htslib/test/test_view.c
        externalTools/htslib/test/thrash_threads1.c
        externalTools/htslib/test/thrash_threads2.c
        externalTools/htslib/test/thrash_threads3.c
        externalTools/htslib/test/thrash_threads4.c
        externalTools/htslib/test/thrash_threads5.c
        externalTools/htslib/test/thrash_threads6.c
        externalTools/htslib/bcf_sr_sort.c
        externalTools/htslib/bgzf.c
        externalTools/htslib/bgzip.c
        externalTools/htslib/errmod.c
        externalTools/htslib/faidx.c
        externalTools/htslib/hfile.c
        externalTools/htslib/hfile_gcs.c
        externalTools/htslib/hfile_libcurl.c
        externalTools/htslib/hfile_net.c
        #externalTools/htslib/hfile_s3.c
        externalTools/htslib/hts.c
        externalTools/htslib/htsfile.c
        externalTools/htslib/kfunc.c
        externalTools/htslib/knetfile.c
        externalTools/htslib/kstring.c
        externalTools/htslib/md5.c
        externalTools/htslib/multipart.c
        #externalTools/htslib/plugin.c
        externalTools/htslib/probaln.c
        externalTools/htslib/realn.c
        externalTools/htslib/regidx.c
        externalTools/htslib/sam.c
        externalTools/htslib/synced_bcf_reader.c
        externalTools/htslib/tabix.c
        externalTools/htslib/tbx.c
        externalTools/htslib/textutils.c
        externalTools/htslib/thread_pool.c
        externalTools/htslib/vcf.c
        externalTools/htslib/vcf_sweep.c
        externalTools/htslib/vcfutils.c
        )

set(SONLIB_SOURCES
        externalTools/sonLib/C/impl/avl.c
        externalTools/sonLib/C/impl/bioioC.c
        externalTools/sonLib/C/impl/commonC.c
        externalTools/sonLib/C/impl/fastCMaths.c
        externalTools/sonLib/C/impl/hashTableC.c
        externalTools/sonLib/C/impl/hashTableC_itr.c
        externalTools/sonLib/C/impl/lz4.c
        externalTools/sonLib/C/impl/lz4hc.c
        externalTools/sonLib/C/impl/pairwiseAlignment.c
        externalTools/sonLib/C/impl/sonLibAlign.c
        externalTools/sonLib/C/impl/sonLibCache.c
        externalTools/sonLib/C/impl/sonLibCommon.c
        externalTools/sonLib/C/impl/sonLibCompression.c
        externalTools/sonLib/C/impl/sonLibConnectivity.c
        externalTools/sonLib/C/impl/sonLibEulerTour.c
        externalTools/sonLib/C/impl/sonLibExcept.c
        externalTools/sonLib/C/impl/sonLibFile.c
        externalTools/sonLib/C/impl/sonLibHash.c
        externalTools/sonLib/C/impl/sonLibKVDatabase.c
        externalTools/sonLib/C/impl/sonLibKVDatabase_BigRecordFile.c
        #externalTools/sonLib/C/impl/sonLibKVDatabase_KyotoTycoon.cpp
        externalTools/sonLib/C/impl/sonLibKVDatabase_MySql.c
        externalTools/sonLib/C/impl/sonLibKVDatabase_TokyoCabinet.c
        externalTools/sonLib/C/impl/sonLibKVDatabaseConf.c
        externalTools/sonLib/C/impl/sonLibList.c
        externalTools/sonLib/C/impl/sonLibMath.c
        externalTools/sonLib/C/impl/sonLibNaiveConnectivity.c
        externalTools/sonLib/C/impl/sonLibRandom.c
        externalTools/sonLib/C/impl/sonLibSet.c
        externalTools/sonLib/C/impl/sonLibSortedSet.c
        externalTools/sonLib/C/impl/sonLibString.c
        externalTools/sonLib/C/impl/sonLibTreap.c
        externalTools/sonLib/C/impl/sonLibTree.c
        externalTools/sonLib/C/impl/sonLibTuples.c
        externalTools/sonLib/C/impl/stGraph.c
        externalTools/sonLib/C/impl/stMatrix.c
        externalTools/sonLib/C/impl/stPhylogeny.c
        externalTools/sonLib/C/impl/stPosetAlignment.c
        externalTools/sonLib/C/impl/stSafeC.c
        externalTools/sonLib/C/impl/stThreadPool.c
        externalTools/sonLib/C/impl/stUnionFind.c
        externalTools/sonLib/C/impl/stJson.c
        externalTools/sonLib/C/impl/jsmn.c
        externalTools/sonLib/C/tests/allTests.c
        externalTools/sonLib/C/tests/cigarsTest.c
        externalTools/sonLib/C/tests/fastaCTest.c
        #externalTools/sonLib/C/tests/kt_connect_test.cpp
        externalTools/sonLib/C/tests/kvDatabaseTest.c
        externalTools/sonLib/C/tests/kvDatabaseTestCommon.c
        externalTools/sonLib/C/tests/sonLibAlignTest.c
        externalTools/sonLib/C/tests/sonLibCacheTest.c
        externalTools/sonLib/C/tests/sonLibCommonTest.c
        externalTools/sonLib/C/tests/sonLibCompressionTest.c
        externalTools/sonLib/C/tests/sonLibConnectivityTests.c
        externalTools/sonLib/C/tests/sonLibDoubleTuplesTest.c
        externalTools/sonLib/C/tests/sonLibEdgeContainerTests.c
        externalTools/sonLib/C/tests/sonLibEulerTest.c
        externalTools/sonLib/C/tests/sonLibExceptTest.c
        externalTools/sonLib/C/tests/sonLibFileTest.c
        externalTools/sonLib/C/tests/sonLibGraphTest.c
        externalTools/sonLib/C/tests/sonLibHashTest.c
        externalTools/sonLib/C/tests/sonLibIntTuplesTest.c
        externalTools/sonLib/C/tests/sonLibListTest.c
        externalTools/sonLib/C/tests/sonLibPosetAlignmentTest.c
        externalTools/sonLib/C/tests/sonLibRandomTest.c
        externalTools/sonLib/C/tests/sonLibSetTest.c
        externalTools/sonLib/C/tests/sonLibSortedSetTest.c
        externalTools/sonLib/C/tests/sonLibStringTest.c
        externalTools/sonLib/C/tests/sonLibTreapTest.c
        externalTools/sonLib/C/tests/sonLibTreeTest.c
        externalTools/sonLib/C/tests/sonLibTuplesTest.c
        externalTools/sonLib/C/tests/stMatrixTest.c
        externalTools/sonLib/C/tests/stPhylogenyTest.c
        externalTools/sonLib/C/tests/stThreadPoolTest.c
        externalTools/sonLib/C/tests/stUnionFindTest.c
        externalTools/sonLib/externalTools/cutest/AllTests.c
        externalTools/sonLib/externalTools/cutest/CuTest.c
        externalTools/sonLib/externalTools/cutest/CuTestTest.c
        externalTools/sonLib/externalTools/quicktree_1.1/src/align.c
        externalTools/sonLib/externalTools/quicktree_1.1/src/aln2dist.c
        externalTools/sonLib/externalTools/quicktree_1.1/src/aln2tree.c
        externalTools/sonLib/externalTools/quicktree_1.1/src/buildtree.c
        externalTools/sonLib/externalTools/quicktree_1.1/src/cluster.c
        externalTools/sonLib/externalTools/quicktree_1.1/src/dist2tree.c
        externalTools/sonLib/externalTools/quicktree_1.1/src/distancemat.c
        externalTools/sonLib/externalTools/quicktree_1.1/src/distancemat_merops.c
        externalTools/sonLib/externalTools/quicktree_1.1/src/options.c
        externalTools/sonLib/externalTools/quicktree_1.1/src/quicktree.c
        externalTools/sonLib/externalTools/quicktree_1.1/src/sequence.c
        externalTools/sonLib/externalTools/quicktree_1.1/src/tree.c
        externalTools/sonLib/externalTools/quicktree_1.1/src/util.c
        externalTools/sonLib/C/inc/sort_r.h)


set(SSW_SOURCES
        externalTools/ssw/ssw.c)

#file(GLOB LP_SOURCES
#        externalTools/lp_solve/lp_lib.c
#        )


file(GLOB LP_SOURCES
        externalTools/lp_solve/*.c
        externalTools/lp_solve/bfp/*.c
        externalTools/lp_solve/bfp_LUSOL/*.c
        externalTools/lp_solve/bfp_LUSOL/LUSOL/*.c
        externalTools/lp_solve/lp_solve/*.c
        externalTools/lp_solve/shared/*.c
        )

set(LP_SOURCES
        externalTools/lp_solve/lp_MDO.c
        externalTools/lp_solve/shared/commonlib.c
        externalTools/lp_solve/colamd/colamd.c
        externalTools/lp_solve/shared/mmio.c
        externalTools/lp_solve/shared/myblas.c
        externalTools/lp_solve/ini.c
        externalTools/lp_solve/lp_rlp.c
        externalTools/lp_solve/lp_crash.c
        externalTools/lp_solve/bfp/bfp_LUSOL/lp_LUSOL.c
        externalTools/lp_solve/bfp/bfp_LUSOL/LUSOL/lusol.c
        externalTools/lp_solve/lp_Hash.c
        externalTools/lp_solve/lp_lib.c
        externalTools/lp_solve/lp_wlp.c
        externalTools/lp_solve/lp_matrix.c
        externalTools/lp_solve/lp_mipbb.c
        externalTools/lp_solve/lp_MPS.c
        externalTools/lp_solve/lp_params.c
        externalTools/lp_solve/lp_presolve.c
        externalTools/lp_solve/lp_price.c
        externalTools/lp_solve/lp_pricePSE.c
        externalTools/lp_solve/lp_report.c
        externalTools/lp_solve/lp_scale.c
        externalTools/lp_solve/lp_simplex.c
        externalTools/lp_solve/lp_SOS.c
        externalTools/lp_solve/lp_utils.c
        externalTools/lp_solve/yacc_read.c
        )

set(SOURCE_FILES
        impl/chunker.c
        impl/column.c
        impl/coordination.c
        impl/emissions.c
        impl/genomeFragment.c
        impl/helenFeatures.c
        impl/hmm.c
        impl/mergeColumn.c
        impl/misc.c
        impl/pairwiseAligner.c
        impl/parser.c
        impl/partitions.c
        impl/poa.c
        impl/profileSeq.c
        impl/stateMachine.c
        impl/view.c
        impl/rle.c
        impl/repeatSubMatrix.c
        impl/stitching.c
        impl/htsIntegration.c
        impl/vcf.c
        impl/profileSeq.c
        impl/bubbleGraph.c
        impl/randomSequences.c
        impl/poa.c
        externalTools/samtools/bedidx.c
        )

###############################
########## LIBRARIES ##########
###############################

# so we can generate margin.so
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Dependencies

find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()

find_package(HDF5)
if (HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIR})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_HDF5 ")
else ()
    message(WARNING "Couldn't find HDF5 library")
endif ()

# MISC
add_library(son ${SONLIB_SOURCES})
add_library(lpsolve ${LP_SOURCES})
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DYY_NEVER_INTERACTIVE -DPARSER_LP -DINVERSE_ACTIVE=INVERSE_LUSOL -DRoleIsExternalInvEngine ")

# HTSLIB
add_custom_command(
        OUTPUT "${PROJECT_SOURCE_DIR}/externalTools/htslib/config.h"
        COMMAND autoconf
        COMMAND autoheader
        COMMAND ./configure --disable-lzma --disable-s3 --disable-plugins --disable-bz2 --disable-libcurl
        COMMAND make
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/externalTools/htslib/
)
add_custom_target(HTSLIB_CONFIGURED DEPENDS "${PROJECT_SOURCE_DIR}/externalTools/htslib/config.h")
add_library(hts ${HTSLIB_SOURCES})
add_dependencies(hts HTSLIB_CONFIGURED)
target_link_libraries(hts pthread z)

# Build margin library
add_library(marginLib ${SOURCE_FILES} ${SSW_SOURCES})
target_link_libraries(marginLib son hts lpsolve)

if (OPENMP_FOUND)
    target_link_libraries(marginLib OpenMP::OpenMP_CXX)
endif ()

if (HDF5_FOUND)
    target_link_libraries(marginLib ${HDF5_C_LIBRARIES})
endif ()

###############################
######### EXECUTABLES #########
###############################

add_executable(marginPolish marginPolish.c)
target_link_libraries(marginPolish marginLib)

add_executable(margin margin.c)
target_link_libraries(margin marginLib)

add_executable(marginVcfPhase marginVcfPhase.c)
target_link_libraries(marginVcfPhase marginLib)

add_executable(marginBamTag tools/marginBamTag.c)
target_link_libraries(marginBamTag marginLib)

enable_testing()

if (HDF5_FOUND)
    add_executable(allTests
            tests/allTests.c
            tests/pairwiseAlignerTest.c
            tests/stRPHmmTest.c
            tests/parserTest.c
            tests/chunkingTest.c
            tests/polisherTest.c
            tests/marginTest.c
            tests/viewTest.c
            tests/featureTest.c
            tests/stitchingTest.c
            tests/vcfTest.c
            )
else ()
    add_executable(allTests
            tests/allTests.c
            tests/chunkingTest.c
            tests/marginTest.c
            tests/pairwiseAlignerTest.c
            tests/stRPHmmTest.c
            tests/parserTest.c
            tests/polisherTest.c
            tests/stRPHmmTest.c
            tests/viewTest.c
            tests/stitchingTest.c
            tests/vcfTest.c)
endif ()
target_link_libraries(allTests marginLib)
