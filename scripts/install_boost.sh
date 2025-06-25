#!/bin/bash


if [ "$#" -ne 2 ]; then
 echo -e "\033[1;31mWrong number of input argumens. \033[1;m"
 echo -e "\033[1;31muse as $0 boost_version dir\033[1;m"
 echo -e "\033[1;31mexample: $0 1_45_0 $HOME/local\033[1;m"
 exit
fi

#===============================
boost_version=$1

# In case the last character in $2 is a slash (/), remove it:
boost_install_dir=${2%\/}/boost
#===============================

boost_basename=boost_$boost_version
boost_archive=$boost_basename.tar.gz

echo "----------------------------------------"
echo "boost version: $boost_version"
echo "boost install dir: $boost_install_dir"
echo "boost basename: $boost_basename"
echo "boost archive name: $boost_archive"
echo "----------------------------------------"

echo -e "\033[1;34mUnpacking the archive ... \033[1;m"
tar xzf $boost_archive
echo -e "\033[1;34m ... unpacking finished \033[1;m"

cd $boost_basename
echo -e "\033[1;34mBootstrapping ... \033[1;m"
./bootstrap.sh --prefix=$boost_install_dir --with-toolset=gcc
echo -e "\033[1;34m... bootstrapping finished \033[1;m"

echo -e "\033[1;34mCompiling ... \033[1;m"
./bjam --build-dir=build toolset=gcc threading=multi --layout=system
echo -e "\033[1;34m... compilation finshed\033[1;m"
#./bjam --build-dir=build toolset=gcc threading=multi -j6
#   or ./bjam -j6 --build-dir=build toolset=gcc --build-type=complete --layout=versioned
#   ( another option is to use threading=single)

echo -e "\033[1;34mInstalling ... \033[1;m"
./bjam install
echo -e "\033[1;34m... installation complete.\033[1;m"


