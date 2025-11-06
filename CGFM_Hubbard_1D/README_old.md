# 1D Hubbard model solver based on the cumulant Green's functions method (CGFM)

Solver for the 1D Hubbard model based on the cumulant Green's functions method (CGFM) (https://doi.org/10.1088/1361-648X/acc628). The method is based on the iterative construction and exact diagonalization of a finite chain containing N correlated sites. The eigenvectors and eigenvalues of the finite chain are used to obtain the atomic Green's functions for the chain using the Lehmann representation. Using the cumulant expansion, the atomic cumulants are obtained from the atomic Green's functions. The atomic cumulants are then used as approximations to the full cumulants to calculate the Green's functions for the infinite chain limit.

### Pre-requisites and how to install them

- gfortran (https://fortran-lang.org/)

  `sudo apt-get install gfortran`
- LAPACK   (https://netlib.org/lapack/)

  `sudo apt-get install liblapack-dev`
- BLAS     (https://netlib.org/blas/)

  `sudo apt-get install libblas-dev `

### How to compile and run the code

- Compiling

  `gfortran -std=legacy -mcmodel=medium -o CGFM_Hubbard_1D.exe CGFM_Hubbard_1D.f90 -llapack -lblas`

  After compiling, the module files and the executable will be created.
  
- Running

  `nohup ./CGFM_Hubbard_1D.exe > CGFM_Hubbard_1D.out &`

  Using this command, the code runs in the background, freeing up the terminal. To check the execution status, use `top`.

  While running, the temporary files "fort.*" are created. They can be removed after the code is done running. After running, the output files containing the density of states and the occupations numbers are created.

### More details

Check the "README.md" inside the "CGFM_Hubbard_1D" directory for a detailed description of the code parameters, variables, subroutines, functions, and output files.

### Description of the important parameters and variables (in order of appearance)

#### In main program

- LFIN:              defines the total number of sites in the chain (LFIN-1 conduction sites connected to 1 impurity);
- IZI:               defines on which site the Green's functions will be calculated (IZI=LFIN means on the last site, which is the impurity site);
- IDOWN:             defines if the down spin Green's functions should be calculated (IDOWN=0 means up spin only; IDOWN=1 means both spins);
- ACOPtI:            electron hopping term between sites of the chain;
- ACOPtR:            electron hopping term in the Green's functions for the lattice (set equal to the chain hopping);
- D:                 bandwidth for the zero order density of states;
- ACOPH:             external magnetic field;
- ACt:               individual hopping terms between sites;
- AMUI:              chemical potential MU;
- BETAI:             BETA=1/T;
- C1,C2,C3,C4:       parameters for the exact double occupation curve;
- AMI:               starting chemical potential amu_i for the chemical potential loop;
- AMF:               ending chemical potential amu_i for the chemical potential loop;
- MT:                number of points for the chemical potential loop;
- UI:                starting electronic correlation U_i for the correlation loop;
- UF:                ending electronic correlation U_f for the correlation loop;
- NT:                number of points for the electron correlation loop;
- ACOPnI:            local energy E_n of the electrons;
- COMPP:             completeness number calculated from up spin Green's functions;
- SVACP:             vacuum occupation number calculated from up spin Green's functions;
- SFUP:              up spin occupation number calculated from up spin Green's functions;
- SFDP:              down spin occupation number calculated from up spin Green's functions;
- SD2P:              double spin occupation number calculated from up spin Green's functions;
- COMPN:             completeness number calculated from down spin Green's functions;
- SVACN:             vacuum occupation number calculated from down spin Green's functions;
- SFUN:              up spin occupation number calculated from down spin Green's functions;
- SFDN:              down spin occupation number calculated from down spin Green's functions;
- SD2N:              double spin occupation number calculated from down spin Green's functions;
- D2:                exact double occupation curve;

#### In subroutine GROUND

- N1:               number of points for the ground-state integration;

#### In subroutine OCUPT

- A(1):             starting frequency for the occupation numbers;
- B(1):             ending frequency for the occupation numbers;

#### In subroutine OCUP

- NT:                number of points for the calculation;
- ETTA:              imaginary part of the complex frequency;
- A:                 starting frequency for the occupation numbers;
- B:                 ending frequency for the occupation numbers;
  
#### In subroutine DENSI2

- NT:                number of points for the calculation;
- ETTA:              imaginary part of the complex frequency;
- OI:                starting frequency for the DOS;
- OF:                ending frequency for the DOS.

### Subroutines and functions (in order of appearance)

- PROGRAM AMHUBBARD: main routine, defines the parameters of the model and calls other subroutines;
- EGGAPBA:           calculates the ground-state energy and the gap of the density of states as functions of the correlation U for the Bethe ansatz solution;
- f(x):              defines the function that represents the gap of the density of states as a function of the correlation U for the Bethe ansatz;
- g(x):              defines the function that represents the ground-state energy as a function of the correlation U for the Bethe ansatz;
- GROUND:            calculates the ground-state energy for the infinite chain;
- OCUPT:             prepares the parameters to calculate the occupation numbers for the infinite chain;
- OCUP:              calculates the occupation numbers for the infinite chain;
- DENSI2:            calculates the density of states for the infinite chain;
- GKONDO:            calculates the frequency dependent up spin and down spin Green's functions for the finite chain using the residues (numerators) from MEOCUP and MEOCDOWN, and afterwards calculates the up spin and down spin Green's functions for the infinite chain;
- RESIDUOS:          builds the finite chain by iteratively adding site by site;
- SIMP:              Simpson integration;
- FERM:              real Fermi function;
- ZFER:              complex Fermi function;
- MEOCUP:            calculates the residues (numerator) of the up spin Green's functions for the finite chain using the Lehmann representation;
- MEOCDOWN:          calculates the residues (numerator) of the down spin Green's functions for the finite chain using the Lehmann representation;
- LEE:               reads fort.* files that contain information about the Fock space from the previous iteration to be used in the current iteration;
- TRANSFNEW:         transfers information about the Fock space from the current iteration to be used in the next iteration;
- DEALLOC_IGSIS:     deallocates memory used to store charge and spin of the Hamiltonian blocks of the current iteration;
- TCUPES1:           calculates the up spin creation operator of the current iteration;
- TCDOWNES1:         calculates the down spin creation operator of the current iteration;
- TCUPSISQ:          calculates the up spin creation operator of the first site;
- TCUPSIS:           calculates the up spin creation operator of the last site;
- TCDOWNSISQ:        calculates the down spin creation operator of the first site;
- TCDOWNSIS:         calculates the down spin creation operator of the last site;
- THSIS:             calculates the non-local elements of the Hamiltonian matrix, adds them to the local elements and diagonalizes the Hamiltonian matrix;
- CONECTIONS:        defines which Hamiltonian blocks connect to each other and to creation operators ;
- ESTRUCTURA:        defines the size, charge, and spin of Hamiltonian blocks;
- DEALLOC_ESP1:      deallocates memory used to store information about the Fock space of the current iteration;
- MULTI:             calculates the local elements of the Hamiltonian matrix;
- VALIN:             defines the operators and Hamiltonian for a single site.

### Output files

- CGFM_Hubbard_1D.out

    Main output file. The file contains information about the execution of the program, such as: input parameters; Hamiltonian blocks with their respective charge, spin, and eigenvalues; lowest eigenvalue; partition function; and others.
  
- YYYYMMDD_aocup_LFIN=X.dat
  
    Output file for the occupation numbers. The name contains the date the program was run (in YYYYMMDD format) and the number of sites in the chain. The file is structured in the following way:

  1. the chemical potential in column 1;
  2. the total vacuum occupation number in column 2;
  3. the total up spin occupation number in column 3;
  4. the total down spin occupation number in column 4;
  5. the total double spin occupation number in column 5;
  6. the total completeness number in column 6;
  7. the total occupation number in column 7.

- YYYYMMDD_adensity_LFIN=X_ACOPU=Y.dat

    Output file for the density of states. The name contains the date the program was run (in YYYYMMDD format), the number of sites in the chain, and the electronic correlation. The file is structured in the following way:

  1. the frequency in column 1;
  2. the up spin density of states in column 2;
  3. the down spin density of states in column 3.
