# WaveDeck

A command deck for ROOT-based digital pulse-processing

## Installation

DISCLAIMER: this code hasn't been tested at all yet. 
There are reasonable expecteation for it to work on Linux systems mathcing 
the following requirements:

**Prerequisite**: `cmake` (> 3.4), `ROOT` (> 6.4) compiled with `C++11` standards, 
typical C++ physicist's tools. It is important that **`ROOT` supports for the 
`fftw3` libraries**: this typically requires that `fftw3` was installed _before_ 
`ROOT` was compiled from source. You can check if your `ROOT` supports
the `fftw3` Fourier transform libriaries by running 
```
root-config --has-fftw3
```
and/or by checking if you can create an instance of `TVirtualFFT`
by running in the `ROOT` prompt
```
$ int n = 1000;
$ TVirtualFFT* fft = TVirtualFFT::FFT(1, &n, "R2C M");
```
If your `ROOT` installation lacks `fftw3` support, then install the 
fftw3 from your software manager, the go back in your ROOT build 
directory and re-run `cmake -DCMAKE_INSTALL_PREFIX=<installdir>`.
You should see in the `cmake` output a message concerning the `fftw3`
library. Alternatively, you can use a version of `fftw3` built-in 
in ROOT whitout installing it on your system. To do so, run 
```
cmake -Dbuiltin_fftw3=ON -DCMAKE_INSTALL_PREFIX=<installdir>
```
(a working internet connection is required).

**Optional**: `doxygen` for building documentations 
```
sudo apt install doxygen doxygen-doc doxygen-gui graphviz
```

To install this package simply clone it and build it using `cmake`
```
$ mkdir build install
$ cd build/
$ cmake -DCMAKE_INSTALL_PREFIX=/path/to/install/directory ..
$ make 
$ make install
```
The install step will copty the `WDeck` dictionaries and libraries in 
the `/lib` folder, while the examples' executables are installed in `/bin`.

## Usage

See the **Examples** pages in the documentation built by `doxygen`. 
The examples source code is available in the `/examples` folder.

## Documentation

If `cmake` has been able to find `doxygen` on your system, 
the documentation of the packaged should have been automatically generated
at build time. 
To access the documentation, simply open the file `docs/html/index.html`
with your favourite browser. 

## Integration

The WaveDeck package is hopefully painless to integrate with other 
ROOT-based code. Depending on your project structure, the user can
pursue different strategies (in growing order of complexity)

### Simple ROOT macro
In case one wants to use WaveDeck in a simple root macro, 
the simples approach to include WaveDeck libraries is the following:
1. Copy the file `src/include/TWDeckPATH.h` into the working directory
2. Add to the `rootlogon.C` script the following lines (create a new one in case no `rootlogon.C` is present)
  ```
  #include "TWDeckPATH.h"

  void rootlogon() {
    gInterpreter->AddIncludePath(WDECK_INCLUDE);
    gSystem->Load(Form("%s/libTWDeckWfm.so", WDECK_LIB_DIR));
    gSystem->Load(Form("%s/libTWDeck.so", WDECK_LIB_DIR));
  }
  ```
  Then in any macro it should be possible to `#include` the WaveDeck classes
  
### Makefile

I've never learnt how to write a proper `Makefile`. Please feel free to help. 

### CMake 

The WaveDeck package comes with a handy `WaveDeckConfig.cmake` file that is 
installed in `bin/WaveDeck/cmake`. If you want your `cmake` project to 
find and source the definition of `WaveDeckConfig.cmake`, then
append to the `CMAKE_PREFIX_PATH` the `path/to/WaveDeck/bin`.
In your project `CMakeLists.txt` you can then use the `find_package` function.


## Bugs and report

If you feel like contributing to this project, report a bug, give any suggestion
or just show some moral support, `daniele<dot>guffanti<at>mib.infn.it`

