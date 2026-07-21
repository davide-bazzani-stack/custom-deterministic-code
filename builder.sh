#!/bin/bash

cmake --preset debug
cmake --build --preset debug --parallel 2

cmake --preset release
cmake --build --preset release --parallel 2

