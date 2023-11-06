ZZPolyCalc: An Open-Source Code with Structure Caching for the Determination of Zhang-Zhang Polynomials of Carbon Nanostructures

ZZPolyCalc computes Zhang-Zhang polynomials for carbon nanostructures. This code is an improved version of the non-cached ZZCalculator (<https://github.com/solccp/zzcalculator>). It retains most of the ZZCalculator's options, significantly enhances speed due to caching, and adds support for five-membered rings.

ZZPolyCalc is released under the GNU General Public License ver. 3. Please consult the included `LICENSE <LICENSE>`_ file for detailed licensing conditions.

Installation Requirements
=========================

* A Fortran 2003 compliant compiler (Intel compiler recommended)
* A C compiler
* GNU make

Additionally, though optional, the following are recommended:

* CMake (version 2.8 or newer)
* `xxHash <https://github.com/Cyan4973/xxHash> library`_ for the XXH128 hash

or

* libSSL library for the MD5 or SHA256 hash.

The code will compile without these hashing libraries but with less optimized hashing code.

Building
========

# It is recommended to use CMake.

# To compile the code:

    mkdir build
    cd build
    cmake ..
    make
    make test
    make install

The code will be installed in bin/ZZPolyCalc.

# To change the compiler:
    CC=icc FC=ifort cmake ..

# To specify a non-standard location of libraries:
    CC=icc FC=ifort CMAKE_PREFIX_PATH=/path/to/xxHash cmake ..

# Alternatively, the code can be compiled by simply issuing:
    make

However, manual modification of the compiler and/or libraries may be necessary.

Usage
=====

For basic usage, a single input with an XYZ file format is required:

```bash
$ bin/ZZPolyCalc test/benzene.xyz

 2 + 1 x
 total: 3
```

This command calculates the ZZ polynomial equal to 2+x for benzene. The total printed corresponds to the Clar number, which is the value of the polynomial when x=1.

The XYZ file can contain any type of atoms, but all atoms apart from carbon will be ignored.

More advanced options can be displayed with:

```bash
$ bin/ZZPolyCalc -h

Usage: ZZPolyCalc [options] input
Options:
    -a                Specifies that the input file contains an adjacency matrix instead of XYZ format
    -m number         Sets the maximum {number} of structures in the cache database
    -p                Prints intermediate bond-level structures
    -Q                Prints the ZZ polynomial in XML format
    -r file           Reads cached structures from a {file}
    -s number         Sets the {number} of buckets in the cache database
    -u                Uses unmodified input XYZ geometry (sorted by default)
    -v                Enables verbose printing
    -w file           Writes cached structures to a {file}
    -h                Displays this help message
```

The -a option indicates that the input geometry file contains an adjacency (connectivity) matrix. The first line should contain the total number of atoms, and the following lines contain:

linenumber atom1 atom2 atom3

where atom1, atom2, and atom3 represent connected atom numbers. If only two atoms are connected, atom3 should be 0.

This option is useful for systems where Cartesian geometry is difficult to obtain, such as fullerenes.

The -m option specifies the number of structures in the cache database. Structures beyond this number will replace previously computed ones.

This option determines memory usage, although the exact memory footprint will depend on the size of the fragments and the size of the ZZ polynomial. Since newer structures may be larger than older ones, after reaching the limit, memory usage may still increase. Therefore, it is recommended to set this limit lower than the total estimated usage.

The -p option prints intermediate structures that are obtained by removing bonds from subsequent structures. This not only shows progress but, for some regular systems, these intermediate fragments represent all systems smaller than the one being calculated.

The -Q option outputs the polynomial in XML format, which may be useful for automatic post-processing.

The -r and -w options are for specifying files to read from and write to the cache, respectively. This is useful for checkpointing or calculating families of structures with similarities. It is only useful for very large systems.

The -s option sets the number of buckets in the cache database. A larger number slightly increases the speed of the program at the expense of more memory. If the database is smaller than the number of calculated structures, it is advisable that the database size is at least 50 times larger than the number of buckets.

The -u option makes ZZPolyCalc read the XYZ file without any modifications. By default, the code sorts the structure by z, y, and x coordinates, respectively, which usually allows the most optimal selection of order for the algorithm.

The -v option enables additional diagnostics and progress printing.


