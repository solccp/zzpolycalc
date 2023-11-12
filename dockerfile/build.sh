#!/bin/bash
apt-get update
apt-get install -y curl zip unzip tar cmake
cd /tmp
git clone  https://github.com/Microsoft/vcpkg.git
cd vcpkg/
./bootstrap-vcpkg.sh
./vcpkg integrate install
./vcpkg install xxhash

cd /tmp
git clone https://github.com/solccp/zzpolycalc.git
cd zzpolycalc/

mkdir build
cd build/

CC=icc FC=ifort CMAKE_PREFIX_PATH=/tmp/vcpkg/installed/x64-linux/ cmake ../
link_file="src/CMakeFiles/ZZPolyCalc.dir/link.txt"

sed -i 's/-lirng //g' $link_file
sed -i 's/-ldecimal //g' $link_file
sed -i 's/-lcilkrts //g' $link_file
sed -i 's/-lstdc++ //g' $link_file

make
mkdir /app
cp src/ZZPolyCalc /app/
rm -rf /tmp


