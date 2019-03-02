#!/bin/bash

root_dir=../
build_dir=../build

mkdir $build_dir

cmake -H$root_dir -B$build_dir
cmake --build $build_dir
