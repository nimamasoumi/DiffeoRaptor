ADD_EXECUTABLE(ImageMatchingTest ImageMatchingTest.cxx)
TARGET_LINK_LIBRARIES(ImageMatchingTest ${PyCA_LIB} GeodesicShooting ${FFTWF_LIB} ${FFTWF_THREADS_LIB} DiffeoRaptor CalcWarpingIndex)

ADD_EXECUTABLE(DiffRapTest2D diffeoRapTest2D.cxx)
TARGET_LINK_LIBRARIES(DiffRapTest2D ${PyCA_LIB} DiffeoRaptor ${FFTWF_LIB} ${FFTWF_THREADS_LIB} CalcWarpingIndex)

ADD_EXECUTABLE(DiffRapIXImulti diffeoRaptorIXImulti.cxx)
TARGET_LINK_LIBRARIES(DiffRapIXImulti ${PyCA_LIB} DiffeoRaptor ${FFTWF_LIB} ${FFTWF_THREADS_LIB} CalcWarpingIndex)

ADD_EXECUTABLE(DiffRapO3 diffeoraptorO3test.cxx)
TARGET_LINK_LIBRARIES(DiffRapO3 ${PyCA_LIB} DiffeoRaptor ${FFTWF_LIB} ${FFTWF_THREADS_LIB} CalcWarpingIndex)

ADD_EXECUTABLE(ProcessO3 processO3results.cxx)
TARGET_LINK_LIBRARIES(ProcessO3 ${PyCA_LIB} DiffeoRaptor ${FFTWF_LIB} ${FFTWF_THREADS_LIB} CalcWarpingIndex)

ADD_EXECUTABLE(DiffRapL2R diffeoraptorL2Rtest.cxx)
TARGET_LINK_LIBRARIES(DiffRapL2R ${PyCA_LIB} DiffeoRaptor ${FFTWF_LIB} ${FFTWF_THREADS_LIB} CalcWarpingIndex)

ADD_EXECUTABLE(ProcessL2R processL2Rresults.cxx)
TARGET_LINK_LIBRARIES(ProcessL2R ${PyCA_LIB} DiffeoRaptor ${FFTWF_LIB} ${FFTWF_THREADS_LIB} CalcWarpingIndex)
