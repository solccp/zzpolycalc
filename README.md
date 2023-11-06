ZZPolyCalc: An open-source code with structure caching for determination of Zhang-Zhang polynomials of carbon nanostructures

ZZPolyCalc calculates Zhang-Zhang polynomials for carbon structures. The code is improved version
of the non-cached ZZCalculator <https://github.com/solccp/zzcalculator>. It preserves most
of the ZZCalculator options, significantly improves the speed due to caching and
adds support for five-membered rings.  

ZZPolyCalc is released under the GNU General Public License ver 3. Please consult
the included `LICENSE <LICENSE>`_ file for the detailed licensing conditions.


Installation requirements
=========================

* A Fortran 2003 compliant compiler (recommended Intel compiler)

* A C-compiler 

* GNU make

Additionally, optional but recommended:

* CMake (version 2.8 or newer)

* `xxHash <https://github.com/Cyan4973/xxHash> library`_ for XXH128 hash

or

* libSSL library for MD5 or SHA256 hash. 

The code will compile with neither of these libraries but with a less optimized code for hashing. 

Building
========

# Recommended to use CMake

# Compile the code:

    mkdir build
    cd build
    cmake ..
    make
    make test
    make install

The code will install to bin/ZZPolyCalc

# Changing the compiler
    CC=icc FC=ifortran cmake ..

# Specifying non-standard location of librares

    CC=icc FC=ifortran CMAKE_PREFIX_PATH=/path/to/xxHash cmake ..

# 

# Alternatively, the code can be compiled with issuing

    make

But manual modification of compiler and/or libraries may be needed.

Usage
======

For basic usage, only a single input with XYZ file format is required

```bash
$ bin/ZZPolyCalc test/benzene.xyz

 2 + 1 x
 total: 3
     ```

calculates the ZZ polynomial equal to 2+x for benzene. The `total' printed corresponds to the
Clar number equal to value of the polynomal for x=1.

XYZ can contain any type of atoms but all atoms apart from carbons will be ignored. 

Other more advanced options can be printed with 

```bash
$ bin/ZZPolyCalc -h

 Usage: ZZPolyCalc [options] input
 Options:
     -a                Input file specifies adjacency instead of Cartesian geometry
     -m number         Maximum {number} of structures in cache database
     -p                Print intermediate bondlevel structures
     -Q                Print the ZZ polynomial in XML format
     -r file           Read cached structures from {file}
     -s number         Use {number} of buckets in cache database
     -u                Use unsorted geometry sorted by default
     -v                Verbose printing
     -w file           Write cached structures to {file}
     -h                Show this message
     ```

-a option specifies that the input geometry file contains adjacency (connectivity) matrix. The first line
should contain the total number of atoms and the following lines contain

linenumber  atom1 atom2 atom3

where atom1, atom2, and atom3 respresent connected atom numbers. If only two atoms are connected, atom3 should be 0. 

This option is useful for systems, where Cartesian geometry is difficult to obtain, such as fullerenes.  

-m specifies the number of structures in the cache database. The structures beyond this number will replace previously computed.

This option determines the memory usage, although exact memory will depend on size of the fragments and the size of the
ZZ polynomial. Since newer structures may be larger than the older ones, after achiving the limit, the memory
usage may still grow and it is recommended that this limit should be specified lower than the total estimated usage. 

-p prints intermediate structures that are obtained by removing bonds of the subseuent structures. Apart from 
showing the progress, for some reqular systems, these intermediate fragments represent all the systems smaller than calculated.

-Q prints the polynomial in XML format, instead of ASCII. It may be usefull for automatic postprocessing. 
 
-r and -w options specify files for reading and writing the cache file. It is usefull for checkpointing
or for calculating families of structures that share similarities. It is usefull only for very large systems. 

-s specifies the number of buckets in the cache database. The larger number increases slightly the speed of the program at the
cost of larger memory. In case the database is smaller than the calculated number of structures, it is recommended that
the databse size is at least 50 times larger than the number of buckets. 

-u makes the ZZPolyCalc read the XYZ specification without any modifications. By default, the code sorts the structure by z, y, and
x coordinates, respectively that usually makes the most optimal selection of order in the algorithm. 

-v prints some additionally diagnostics and prints progress 
