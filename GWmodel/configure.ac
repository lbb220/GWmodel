AC_INIT([GWmodel], 2.2-0)

AC_LANG(C++)
AC_REQUIRE_CPP

echo "Checking for C++ Compiler"
AC_PROG_CXX

ENABLE_CUDA=0

AC_MSG_CHECKING("Checking environment variable CUDA_HOME")
if test -z "${CUDA_HOME}"; then
    CUDA_HOME=`find /usr/local/ -maxdepth 1 -type d -name "cuda*" | sort -V | tail -1`
    AC_MSG_RESULT("CUDA_HOME not set; using highest version found ${CUDA_HOME}")
else
    AC_MSG_RESULT("using CUDA_HOME=${CUDA_HOME}")
fi

AC_CHECK_FILE([${CUDA_HOME}/bin/nvcc], [HAS_NVCC="yes"])
if test -z "${HAS_NVCC}"; then
    AC_MSG_RESULT("NVCC not found")
else
    AC_MSG_RESULT("NVCC found")
	ENABLE_CUDA=1
fi

AC_ARG_ENABLE([cuda], [AS_HELP_STRING([--enable-cuda],[install cuda (default yes)])],
[ case "${enableval}" in
  yes) ENABLE_CUDA=1 ;;
  no)  ENABLE_CUDA=0 ;;
  *)   AC_MSG_ERROR(bad value ${enableval} for --enable-cuda) ;;
esac]
)

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

# find R home and set correct compiler + flags
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  AC_MSG_ERROR([cannot determine R_HOME. Make sure you use R CMD INSTALL!])
fi
AC_MSG_NOTICE([R_HOME: ${R_HOME}])
RBIN="${R_HOME}/bin/R"

R_INCL=`"${R_HOME}/bin/R" CMD config --cppflags`
R_LIBS=`"${R_HOME}/bin/R" CMD config --ldflags`
R_CPIC=`"${R_HOME}/bin/R" CMD config CPICFLAGS`

## look for Rscript, but use the one found via R_HOME to allow for multiple installations
AC_DEFUN([AC_PROG_RSCRIPT], [AC_CHECK_PROG(RSCRIPT,Rscript,yes)])
AC_PROG_RSCRIPT

if test x"${RSCRIPT}" == x"yes" ; then
    rcpp_incl=`${R_HOME}/bin/Rscript -e 'Rcpp:::CxxFlags()'`
    rcpp_armadillo_incl=`${R_HOME}/bin/Rscript -e "RcppArmadillo:::CxxFlags()"`
		AC_SUBST(RCPP_INCL,["$rcpp_incl"])
		AC_SUBST(RCPP_ARMA_INCL,["$rcpp_armadillo_incl"])
else
    echo "
		Your installation does not appear to have Rscript installed.
		
		Please make sure that you have a working and complete R installation.
		"
    exit 1
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

if test -z `echo "${R_PACKAGE_DIR}" | egrep '00LOCK-.*?/00new'`; then
	echo "You not use staged install"
	AC_SUBST(R_PACKAGE_DIR_REAL,["$R_PACKAGE_DIR"])
	AC_SUBST(R_PACKAGE_DIR)
else
	echo "You use staged install"
	R_PACKAGE_DIR_REAL=`echo "${R_PACKAGE_DIR}" | sed -e 's|/00LOCK-.*/00new||gi'`
	echo "Your real libray is ${R_PACKAGE_DIR_REAL}"
	AC_SUBST(R_PACKAGE_DIR_REAL)
	AC_SUBST(R_PACKAGE_DIR)
fi


AC_SUBST(GWmodel_CUDA_LIB_PATH,["~/.cache/GWmodel"])

AC_MSG_NOTICE([Building Makevars])
AC_SUBST(CU_INCL)
AC_SUBST(NVCC)
AC_SUBST(R_INCL)
AC_SUBST(R_CPIC)
AC_SUBST(R_LIBS)
AC_SUBST(PKG_LIBS)
AC_SUBST(CUDA_INCL)
AC_SUBST(ENABLE_CUDA)
AC_OUTPUT(src/Makevars)
