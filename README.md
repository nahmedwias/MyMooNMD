ParMooN - Parallel Mathematics and object oriented Numerics
==========

Introduction
------------

This is the introduction.


Quick start with ParMooN
------------

You are eager to get started and have no no time to read any verbose documentation?

This is the minimal out-of-source build instruction.
Write the following into a shell 

    mkdir parmoon
    cd parmoon
    mkdir code
    hg clone ssh://username@repos.wias-berlin.de/ParMooN code
    mkdir build
    cd build/
    cmake ../code
    make check doc -j

The '-j' for using all available cores is optional.


Recommended start with ParMooN
----------

ParMooN depends on and uses several external libraries. These partly depend on
each other in somewhat difficult ways. In order to facilitate the process as
much as possible you can use PETSc to handle external libraries completely. This
is the recommended way even though it is also possible to use libraries
installed on your system. The following guidelines create a certain directory
structure which can of course be changed. First download and configure PETSc as
given below

    cd <path to software directory>
    mkdir petsc
    git clone -b maint https://bitbucket.org/petsc/petsc petsc
    cd petsc
    ./configure --with-cc=gcc --with-cxx=g++ --with-fc=gfortran --with-debugging=0 COPTFLAGS='-O3 -march=native -mtune=native' CXXOPTFLAGS='-O3 -march=native -mtune=native' FOPTFLAGS='-O3 -march=native -mtune=native' --download-fblaslapack --download-openmpi --download-scalapack --download-mumps --download-metis --download-suitesparse --download-parmetis

This will take time and requires internet access because many external libraries
are downloaded and compiled. Then you have to compile PETSc itself, the exact
command is written at the end of the configure output. Something like

    make PETSC_DIR=<path to software directory>/petsc PETSC_ARCH=arch-linux2-c-opt all

This again will take time. Finally you can (but need not to) check the
installation, again the exact command is given at the end of the previous output

    make PETSC_DIR=<path to software directory>/petsc PETSC_ARCH=arch-linux2-c-opt check

Next ParMooN can be set up using the PETSC_DIR and PETSC_ARCH arguments:

    cd <path to software directory>
    mkdir parmoon
    cd parmoon
    mkdir code
    hg clone ssh://username@repos.wias-berlin.de/ParMooN code
    mkdir documentation
    mkdir build
    cd build
    cmake  -DPETSC_DIR=<path to software directory>/petsc -DPETSC_ARCH=arch-linux2-c-opt -DPARMOON_DOCUMENTATION_DIRECTORY=<path to software directory>/parmoon/documentation ../code

This will configure ParMooN. Then in the build directory you can compile. As a
start it makes sense to compile and run all the unit and regression tests:

    make check

Use the option '-j4' to use 4 parallel builds to speed up this process. Next

    make doc

will build the Doxygen documentation, open the file 
<path to software directory>/parmoon/documentation/html/index.html in a web
browser to see it. There is more information there for you to read.

