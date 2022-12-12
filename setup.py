from os import path
from setuptools import setup, find_packages
from distutils.extension import Extension
from Cython.Build import cythonize

here = path.abspath(path.dirname(__file__))
package_name = 'zsmash'

version = {}
with open('version.py') as fp:
    exec(fp.read(), version)

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

extensions = [Extension(
    name=f'{package_name}.bootstrap',
    sources=[
        f'{package_name}/bootstrap.pyx',
        f'{package_name}/smash.pyx',
    ],
    extra_compile_args=[
        '-std=c++11',
        '-O3',
        '-static',
        '-fopenmp',
        '-I./zbase',
        '-I./boost',
        '-Wl,--as-needed',
        '-DCYTHON_PEP489_MULTI_PHASE_INIT=0', # switch multi-phase initialization off which is default in Cython 0.29 for Python>=3.5
    ],
    libraries=[
        'semcrct', 'configfile', # recompiled from zbase/lib with -fPIC
        'gsl', 'gslcblas',
        'gomp',
        'm',
        'boost_system',
        'boost_thread',
        'boost_program_options',
        'boost_timer',
        'boost_chrono',
    ],
    extra_link_args=[
        '-std=c++11',
        '-fopenmp',
        '-static-libgcc',
        '-static-libstdc++',
        '-O3',
        '-Wl,--as-needed',
    ],
    library_dirs=[
        '.',
        f'./{package_name}',
        './boost/lib',
        '/usr/lib64', # use gsl libs here instead of in zbase
        './zbase/lib',
    ],
    language='c++'
)]

setup(
    name=package_name,
    author='zed.uchicago.edu',
    author_email='ishanu@uchicago.edu',
    version = str(version['__version__']),
    packages=find_packages(),
    scripts=[],
    url='https://github.com/zeroknowledgediscovery/python_implementations_',
    license='LICENSE',
    description='Cythonized Libraires of ZeD Code',
    keywords=[
        'cython',
        'machine learning',
    ],
    # download_url='https://github.com/zeroknowledgediscovery/python_implementations_/archive/'+str(version['__version__'])+'.tar.gz',
    long_description=long_description,
    long_description_content_type='text/x-rst',
    ext_modules=cythonize(
        extensions,
        language_level='3' # set python3
    ),
    install_requires=[
        'Cython>=0.29.23'
        'numpy',
        'pandas',
    ],
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Software Development :: Libraries',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6'],
    include_package_data=True,
)
