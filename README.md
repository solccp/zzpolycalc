ZZPolyCalc: An Open-Source Code with Fragment Caching for the Determination of Zhang-Zhang Polynomials of Carbon Nanostructures

ZZPolyCalc computes Zhang-Zhang polynomials for carbon nanostructures. This code is an overhauled and improved version of the 
[non-cached ZZCalculator](https://github.com/solccp/zzcalculator). While retaining most of the ZZCalculator's features and command-line options, it is significantly faster due to caching, and can also work for fullerenes, owing to the added support for five-membered rings. ZZPolyCalc achieves polynomial scaling for quasi-1D systems, such as elongated graphene flakes or carbon nanotubes, thereby justifying the term 'Poly' in its name.

ZZPolyCalc is released under the GNU General Public License ver. 3. Please consult the included [LICENSE](LICENSE) file for detailed licensing conditions.

Installation Requirements
=========================

* A Fortran 2008 compliant compiler (Intel compiler recommended)
* A C compiler
* GNU make

Additionally, though optional, the following are recommended:

* CMake (version 2.8 or newer)
* [xxHash library](https://github.com/Cyan4973/xxHash) for the XXH128 hash or [OpenSSL library](https://github.com/openssl/openssl) for the MD5 or SHA256 hash.

If these hashing libraries are not provided, the code will compile with a less optimized MD5 hashing algorithm.
The XXH128 hash is the fastest, followed by MD5. XXH128 is recommended for regular use. To employ SHA256, an appropriate selection must be manually activated in CMakeLists.txt. (See: Optional modifications at compile time)

Obtaining the Source
====================

To clone the repository:

```bash
git clone https://github.com/quantumint/zzpolycalc
cd zzpolycalc
```

Building
========

It is recommended to use CMake.

## To compile the code:

```bash
    mkdir build
    cd build
    cmake ..
    make
    make test
    make install
```

The code will be installed as bin/ZZPolyCalc.

## To change the compiler:

```bash
    CC=icc FC=ifort cmake ..
```

## To specify a non-standard location of libraries:

```bash
    CC=icc FC=ifort CMAKE_PREFIX_PATH=/path/to/xxHash cmake ..
```

## Alternatively, the code can be compiled by simply issuing:

```bash
    make
```

However, manual modification of the compiler and/or libraries may be necessary.

Usage
=====

For basic usage, a single input with an XYZ file format is required:

```bash
$ bin/ZZPolyCalc test/benzene.xyz

 2 + 1 x
 total: 3
```

This command calculates the ZZ polynomial of benzene equal to 2+x. The term 'total:' corresponds to the number of Clar covers, which is equal to the value of the polynomial at x=1.

The XYZ file can contain any type of atoms, but all atoms apart from carbon will be ignored.

More advanced options can be displayed with:

```
$ bin/ZZPolyCalc -h

Usage: ZZPolyCalc [options] input
Options:
    -a                Specifies that the input file contains an adjacency matrix instead of XYZ format
    -c number         Changes cache status printing at verbose mode to every {number} million steps
    -f number         Sets the frequency of cache writes at {number} of million of structures. Requires -w
    -m number         Sets the maximum {number} of structures in the cache database
    -p                Prints intermediate bond-level structures
    -Q                Prints the ZZ polynomial in XML format
    -r file           Reads cached structures from a {file}
    -s number         Sets the {number} of buckets in the cache database
    -u                Uses unmodified input XYZ geometry (sorted by default)
    -v                Enables verbose printing
    -w file           Writes cached structures to a {file}
    -X                Reads connection table from the bottom of the XYZ file
    -h                Displays this help message
```

The -a option indicates that the input geometry file contains an adjacency (connectivity) matrix. The first line should contain the total number of atoms $`n_{at}`$ present in the system, and each of the following lines contains:

`atom_label` `atom1` `atom2` `atom3`

where `atom1`, `atom2`, and `atom3` represent the label of atoms connected to the atom `atom_label`. The labels `atom1`, `atom2`, `atom3`, and `atom_label` assume numerical values $`\in \left\{1,2,...,n_{at}\right\}`$. The labels `atom_label` need to be ordered from 1 to $`n_{at}`$. If only two atoms are connected to the atom `atom_label`, `atom3` should be set to 0. If only one atom is connected to the atom `atom_label`, `atom2` and `atom3` should be set to 0.

This option is useful for systems where Cartesian geometry is not readily available (e.g., for fullerenes).

The -c option allows users to modify the frequency of cache status printouts. The frequency is measured in millions of processed fragments and can be a positive non-integer value.

The -f option sets the frequency of cache checkpointing to a file indicated by the -w flag. This is particularly useful for very long jobs and for jobs that may need restarting before completion.

The -m option specifies the maximal number of structures kept in the cache database. Structures beyond this number will replace previously computed ones.

This option determines memory usage, although the exact memory footprint will depend on the size of the fragments and the size of the ZZ polynomial. Since newer structures may be larger than older ones, after reaching the limit, memory usage may still increase. Therefore, it is recommended to set this limit lower than the total estimated usage.

The -p option prints intermediate structures that are obtained by removing bonds from subsequent structures. This not only shows progress but, for some regular systems, these intermediate fragments represent all systems smaller than the one being calculated.

The -Q option outputs the polynomial in XML format, which may be useful for automatic post-processing.

The -r and -w options are for specifying files which the cache is read from and/or written to, respectively. This is useful for checkpointing or calculating families of similar structures. It is only useful for very large systems.

The -s option sets the number of buckets in the cache database. A larger number slightly increases the speed of the program at the expense of memory. If the database is smaller than the number of calculated structures, it is advisable that the database size is at least 50 times larger than the number of buckets.

The -u option makes ZZPolyCalc read the XYZ file without any modifications. By default, the code sorts the structure by the z, y, and x coordinates, respectively, which usually allows the most optimal selection of order for the algorithm. Note that proper orientation of an elongated benzenoid is crucial for efficient execution of ZZPolyCalc; such a structure should be elongated in the z direction.

The -v option enables additional diagnostics and progress printing.

The -X option causes ZZPolyCalc to bypass the default generation of connectivity tables and instead read a connectivity table appended to the bottom of the XYZ file. This feature ensures compatibility with ZZCalculator. The format of the connectivity tables is as follows:

  number_of_bonds
  `atom1` `atom2`
  ...

where the atoms `atom1` and `atom2` are connected by a bond. Activating this option implicitly enables the `-u` option. The provided bond list influences the traversal path taken during the ZZ polynomial calculation and can significantly affect the computation time.

Optional modifications at compile time
======================================

Most options are selected at runtime. However, during compilation, there is an automatic selection process for a hashing library. If xxHash is available, it will be chosen first, followed by MD5 from the OpenSSL library. If neither is available, an MD5 implementation is compiled in. To select the SHA256 hash when OpenSSL is available, set the option in CMakeLists.txt as follows:

option(USE_SHA256 "Enable SHA256" ON)

The use of SHA256 hash is not recommended, as already the 128-bit xxHash or MD5 results in an astronomically small chance of collision.

In src/types_module.f90, the parameter `vlongmax` (with a default value of 51) specifies the maximum size of the big integer available during the calculations. If the verbose option (-v) is enabled, the size (less or equal to `vlongmax`) that was actually used in the calculation process will be printed at the end of the program execution. For very large systems, it may be necessary to increase this parameter and recompile the program.

