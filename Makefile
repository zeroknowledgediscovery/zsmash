OBJ = cythonized
all: $(OBJ)

# compile binaries --------------------------------

cythonized: setup.py zedsuite/*.pyx
	python3 setup.py build_ext --inplace; rm -rf zedsuite/*.cpp build zedsuite.egg-info zedsuite/__pycache__

wheels: setup.py zedsuite/*.pyx
	python3 setup.py sdist bdist_wheel; rm -rf zedsuite/*.cpp build zedsuite.egg-info zedsuite/__pycache__

# utility commands --------------------------------

clear:
	rm -rf *.o *.so
clean:
	rm -rf *.o *~ $(OBJ) *.tgz *.so

 # EOF --------------------------------
