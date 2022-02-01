# WaveDeck

A command deck for ROOT-based digital pulse-processing

## Installation

DISCLAIMER: this code hasn't been tested at all yet. 
There are reasonable expecteation for it to work on Linux systems mathcing 
the following requirements:

**Prerequisite**: `cmake` (> 3.4), `ROOT` (> 6.4) compiled with `C++11` standards, 
typical C++ pysicist's tools

**Optional**: `doxygen` for building documentations 
```
sudo apt install doxygen doxygen-doc doxygen-gui graphviz
```

To install this package simply clone it and built it using `cmake`
```
$ cd build/
$ cmake ..
$ make 
$ make install
```

## Usage

See the **Examples** pages in the documentation buil by `doxygen`

## Documentation

If `cmake` has been able to find `doxygen` on your system, 
the documentation of the packaged should have been automatically generated
at build time. 
To access the documentation, simply open the file `docs/html/index.html`
with your favourite browser. 

