SET(Base_SRC
    Base/FftOper.cxx
    Base/FieldComplex3D.cxx
    Base/Armalib.cxx
    Base/MPIlib.cxx
    )

ADD_LIBRARY(GeodesicShooting ${Base_SRC} Alg/GeodesicShooting.cxx)
TARGET_LINK_LIBRARIES(GeodesicShooting ${PyCA_LIB} ${FFTWF_LIB} ${FFTWF_THREADS_LIB})

ADD_LIBRARY(Tpca ${Base_SRC} Alg/Tpca.cxx)
TARGET_LINK_LIBRARIES(Tpca ${PyCA_LIB} ${ARMADILLO_LIBRARIES} ${FFTWF_LIB} ${FFTWF_THREADS_LIB})


ADD_LIBRARY(CalcWarpingIndex ${Base_SRC} CalcWarpingIndex.cxx)
TARGET_LINK_LIBRARIES(CalcWarpingIndex ${PyCA_LIB} ${FFTWF_LIB} ${FFTWF_THREADS_LIB})

ADD_LIBRARY(DiffeoRaptor ${Base_SRC} DiffeoRaptor.cxx)
TARGET_LINK_LIBRARIES(DiffeoRaptor ${PyCA_LIB} ${FFTWF_LIB} ${FFTWF_THREADS_LIB})
