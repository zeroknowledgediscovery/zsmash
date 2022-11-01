#!/bin/bash
skip_install_deps=0
testing=0
OPTIND=1
while getopts ft opt; do
    case $opt in
        f) skip_install_deps=1 ;;
        t) testing=1 ;;
    esac
done
shift "$((OPTIND-1))"

# Check dependencies, use "-f" to skip
if [[ ${skip_install_deps} == 0 ]]; then
    if [[ $(yum list installed | grep gsl-devel | wc -l) == "0" ]]; then
        yum -y install gsl-devel || (echo "Failed to install gsl"; exit)
    fi

    if [[
        $(ld -lboost_system -L./boost/lib --verbose | grep succeed | wc -l) == "0"
        || $(ld -lboost_timer -L./boost/lib --verbose | grep succeed | wc -l) == "0"
        || $(ld -lboost_program_options -L./boost/lib --verbose | grep succeed | wc -l) == "0"
        || $(ld -lboost_thread -L./boost/lib --verbose | grep succeed | wc -l) == "0"
        || $(ld -lboost_chrono -L./boost/lib --verbose | grep succeed | wc -l) == "0"
    ]]; then
        echo "Compiling Boost..."
        cd ~/
        cp -r /opt/src/boost_src/ . || (echo "Boost not found. Download from boost.org and place in directory as \"boost_src\""; exit)
        cd boost_src
        ./bootstrap.sh --with-libraries=program_options,system,chrono,timer,thread
        ./b2 -j8 cxxflags="-fPIC" runtime-link=static variant=release link=static
        cp stage/lib/libboost_* /usr/local/lib
        cd ~/
        cp -r boost_src/boost /usr/local/include
        rm -rf boost
    fi
fi

PKGNAME="zedsuite"
VERSION=$(cat version.py | awk '{print $3}' | sed -e 's/^.//' -e 's/.$//')

# Remove old
# rm -rf dist && echo "Removed dist/"

# Standard paths to python on manylinux docker image
PIP310=/opt/python/cp310-cp310/bin/pip
PIP39=/opt/python/cp39-cp39/bin/pip
PIP38=/opt/python/cp38-cp38/bin/pip
PIP37=/opt/python/cp37-cp37m/bin/pip
PIP36=/opt/python/cp36-cp36m/bin/pip
PYTHON310=/opt/python/cp310-cp310/bin/python
PYTHON39=/opt/python/cp39-cp39/bin/python
PYTHON38=/opt/python/cp38-cp38/bin/python
PYTHON37=/opt/python/cp37-cp37m/bin/python
PYTHON36=/opt/python/cp36-cp36m/bin/python

if [[ ${testing} == 1 ]]; then
    ${PIP310} install cython
    ${PYTHON310} setup.py sdist bdist_wheel
    auditwheel repair "dist/${PKGNAME}-${VERSION}-cp310-cp310-linux_x86_64.whl"
    rm "dist/${PKGNAME}-${VERSION}-cp310-cp310-linux_x86_64.whl"
else
    # Compile wheel for each version and install cython if not already installed
    declare -a PYTHON_BINS=(${PYTHON310} ${PYTHON39} ${PYTHON38} ${PYTHON37} ${PYTHON36})
    declare -a PIP_BINS=(${PIP310} ${PIP39} ${PIP38} ${PIP37} ${PIP36})
    declare -a PY_VERSIONS=("cp310-cp310" "cp39-cp39" "cp38-cp38" "cp37-cp37m" "cp36-cp36m")
    NUM_PY_VERSIONS=${#PY_VERSIONS[@]}
    for (( i=0; i<${NUM_PY_VERSIONS}; i++ ));
    do
        ${PYTHON_BINS[i]} setup.py sdist bdist_wheel || (${PIP_BINS[i]} install cython && ${PYTHON_BINS[i]} setup.py sdist bdist_wheel)
        auditwheel repair "dist/${PKGNAME}-${VERSION}-${PY_VERSIONS[i]}-linux_x86_64.whl"
        rm "dist/${PKGNAME}-${VERSION}-${PY_VERSIONS[i]}-linux_x86_64.whl"
    done
fi

# Replace wheels with repaired ones
cp -r wheelhouse/* dist && rm -rf wheelhouse build ${PKGNAME}/*.cpp ${PKGNAME}.egg-info a.out && echo "Copied wheels to dist/"
