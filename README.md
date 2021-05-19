# A. Gierzkiewicz, P. Zgliczyński, "Periodic orbits in Rössler system"

## An installation instruction of the programs realizing numerical part of the proofs of periodic orbits existence/nonexistence and covering relations between h-sets on a Poincaré section.

### Requirements:
The program is written in C++ and has been tested under Linux Mint 18.1 with gcc (Ubuntu 5.4.0-6ubuntu1~16.04.12) compiler. The program uses the CAPD library ver. 5.0.6 (see also [sourceforge download](https://sourceforge.net/projects/capd/files/) zone). The CAPD library is also available as debian deb package. There are also prebuilt versions for Debian, Ubuntu and OSX.

### Compilation instruction:
In principle it is possible to compile the CAPD library under MS Windows but we strongly recommend to compile and run the programs on linux-like systems. The following commands should be executed from the terminal (or msys environment under MS Windows)

- unzip the archive: ```unzip Roessler_periodic.zip```
- enter directory of the program: ```cd Roessler_periodic```
- build the programs: ```make CAPDBINDIR=relative_or_absolute_path_to_capd_bin_directory```

    for example: ```make CAPDBINDIR=../../capd-build/bin/```
    
    Warning: Do not forget the last slash character.
- If the CAPD library is installed in standard directories that are on system path (like /usr/bin/) then you can compile the program just by invoking ```make```
    
    The command generates executable files in the current directory.

### Programs options:

The programs can be run from the terminal, for example by ```./01-Roessler_a525```

There are two programs with computer-assisted proofs in the package:

- 01-Roessler_a525: the main menu asks to choose one of two procedures:
1. Proof of a 3-periodic orbit's existence by INM;
2. Proof of the chain horizontal covering relations:
        N0 =Pc⇒ N1 =Pc⇒ N1 =Pc⇒ N0. 
- 02-Roessler_a47: the main menu asks to choose one of six procedures:
1. Proof of a 5-periodic orbit's existence by INM
2. Proof of the chain horizontal covering relations:
        N0 =Pc⇒ N2 =Pc⇒ N2 =Pc^3⇒ N0.
3. Proof of a 2-periodic orbit's existence by INM
4. Proof of a 4-periodic orbit's existence by INM
5. Proof of forward invariance of A
6. Proof of non-existence of 3-periodic orbits by INM 

### Output:

All procedures execute within seconds on a laptop type computer with Intel Core i7 2.7GHz x 2 processor.

The programs do not create output files. They print to the screen what conditions are being checked and the result of this check (true or false). 
