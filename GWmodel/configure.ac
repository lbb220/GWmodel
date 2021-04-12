AC_INIT([GWmodel], 2.2-0)

AC_LANG(C++)
AC_REQUIRE_CPP

echo "Checking for C++ Compiler"
AC_PROG_CXX

ENABLE_CUDA=0


AC_ARG_ENABLE([cuda], [AS_HELP_STRING([--enable-cuda],[install cuda (default no)])],
[ case "${enableval}" in
	yes) ENABLE_CUDA=1 ;;
	no)  ENABLE_CUDA=0 ;;
	*)   AC_MSG_ERROR(bad value ${enableval} for --enable-cuda) ;;
esac])

if test ${ENABLE_CUDA} -eq 1; then
	AC_MSG_CHECKING("Checking environment variable CUDA_HOME")
	if test -z "${CUDA_HOME}"; then
		AC_MSG_RESULT("CUDA_HOME not set, stop finding nvcc")
	else
		AC_MSG_RESULT("using CUDA_HOME=${CUDA_HOME}")
	fi

	if test -z "${CUDA_HOME}"; then
		if test ${ENABLE_CUDA} -eq 1; then
			echo "
				It seams that you don't have CUDA environment and the --enable-cuda is yes

				Please make sure your cuda environment is under the /usr/local directory

				If you don't  want to install cuda function, please set `--enable-cuda=no`
			"
			exit 1
		fi
	fi

	if test ${ENABLE_CUDA} -eq 1; then
		AC_MSG_CHECKING("whether this is the 64 bit linux version of CUDA")
		AC_CHECK_FILE([${CUDA_HOME}/lib64/libcudart.so], [BIT64="yes"])
		if test -z "${BIT64}"; then
				AC_MSG_RESULT("no -- using ${CUDA_HOME}/lib for CUDA libs")
				cu_libdir="${CUDA_HOME}/lib"
				AC_SUBST(CUDA_LIBDIR,["$cu_libdir"])
		else
				AC_MSG_RESULT("yes -- using ${CUDA_HOME}/lib64 for CUDA libs")
				cu_libdir="${CUDA_HOME}/lib64"
				AC_SUBST(CUDA_LIBDIR,["$cu_libdir"])
		fi
	else
	AC_SUBST(CUDA_LIBDIR,[""])
	fi
fi

if test ${ENABLE_CUDA} -eq 1; then
    AC_MSG_RESULT(You have CUDA environment and we will compile the cuda function for you)
	cu_sources="GWmodelKernel.cu"
	cu_sharedlibs="GWmodelKernel.o"
	cpp_sources_cuda="CGWmodelCUDA.cpp IGWmodelCUDA.cpp"
	cpp_sharedlibs_cuda="CGWmodelCUDA.o IGWmodelCUDA.o"
	OBJECTS_CUDA="${cu_sharedlibs} ${cpp_sharedlibs_cuda}"

	CDEFS="-DCUDA_ACCE"
	TARGET_CUDA_LIB="libgwmodelcuda.so"

	CUDA_LIBS="-L${CUDA_LIBDIR} -Wl,-rpath=${CUDA_LIBDIR} -lcudart -lcublas "
	PKG_LIBS=${PKG_LIBS}"${CUDA_LIBS} -L${HOME}/.cache/GWmodel -Wl,-rpath=${HOME}/.cache/GWmodel  -lgwmodelcuda"
	AC_SUBST(cu_sources)
	AC_SUBST(cu_sharedlibs)
	AC_SUBST(cpp_sources_cuda)
	AC_SUBST(cpp_sharedlibs_cuda)
	AC_SUBST(OBJECTS_CUDA)
	AC_SUBST(CDEFS)
	AC_SUBST(TARGET_CUDA_LIB)
	AC_SUBST(CUDA_LIBS)
fi


AC_MSG_RESULT(building the nvcc command line)
NVCC="${CUDA_HOME}/bin/nvcc"
AC_MSG_RESULT(building the cuda include path)
CUDA_INCL="${CUDA_HOME}/include"

AC_SUBST(GWmodel_CUDA_LIB_PATH,["~/.cache/GWmodel"])

AC_MSG_NOTICE([Building Makevars])
AC_SUBST(CU_INCL)
AC_SUBST(NVCC)
AC_SUBST(PKG_LIBS)
AC_SUBST(CUDA_INCL)
AC_SUBST(ENABLE_CUDA)
AC_OUTPUT(src/Makevars)
