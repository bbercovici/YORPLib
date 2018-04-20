## YORPLib

A library of functions enabling one to compute the coefficients of the Fourier decomposition of the force, torque induced by solar insulation over a polyhedral surface.

OpenMP-compliant compilers enable faster performance. If GCC has been downloaded from Homebrew, OpenMP can be enabled like so:

- `brew install yorplib --with-gcc` if using Homebrew
- `cmake .. -DUSE_GCC:BOOL=TRUE` otherwise

## Installation: 

### Mac users

YORPLib can be retrieved from Homebrew:

    brew tap bbercovici/self
    brew update
    brew install yorplib

### Linux & Mac users

    git clone https://github.com/bbercovici/YORPLib
    cd YORPLib/build
    cmake ..
    make
    make install

## Getting updates

### Mac users

Assuming that YORPLib was installed with Homebrew

<pre>
brew update
brew upgrade yorplib
</pre>

### Linux & Mac users

    git pull
    cd build
    cmake ..
    make
    make install


## Credits

* Originally created by Jay McMahon on 5/19/14 
* Adapted by Benjamin Bercovici on 03/10/2018
* Copyright (c) 2014 Jay McMahon. All rights reserved.

## License

[This software is distributed under the MIT License](https://choosealicense.com/licenses/mit/)
