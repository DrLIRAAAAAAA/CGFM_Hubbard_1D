Co-authored-by: M. S. Figueira <marcosfigueira@id.uff.br>
Co-authored-by: J. Silva-Valencia <jsilvav@unal.edu.co>





# Anderson impurity model (AIM) solver based on the cumulant Green's functions method (CGFM)

Solver for the Anderson impurity model (AIM) based on the cumulant Green's functions method (CGFM) (https://arxiv.org/abs/2409.16881). The method is based on the iterative construction and exact diagonalization of a finite chain containing a single correlated impurity site connected to a chain of uncorrelated sites. The eigenvectors and eigenvalues of the finite chain are used to obtain the atomic Green's functions for the chain using the Lehmann representation. Using the cumulant expansion, the atomic cumulants are obtained from the atomic Green's functions. The atomic cumulants are then used as approximations to the full cumulants to calculate the Green's functions for the infinite chain limit.

### Pre-requisites and how to install them

- gfortran (https://fortran-lang.org/)

  `sudo apt-get install gfortran`
  
- LAPACK   (https://netlib.org/lapack/)

  `sudo apt-get install liblapack-dev`
  
- BLAS     (https://netlib.org/blas/)

  `sudo apt-get install libblas-dev `

### How to compile and run the code

- Compiling

  `gfortran -std=legacy -mcmodel=medium -o CGFM_Anderson_impurity_1D.exe CGFM_Anderson_impurity_1D.f90 -llapack -lblas`

  After compiling, the module files and the executable will be created.
  
- Running

  `nohup ./CGFM_Anderson_impurity_1D.exe > CGFM_Anderson_impurity_1D.out &`

  Using this command, the code runs in the background, freeing up the terminal. To check the execution status, use `top`.

  While running, the temporary files "fort.*" are created. They can be removed after the code is done running. After running, the output files containing the density of states and the occupations numbers are created.

### More details

Check the "README.md" inside the "CGFM_Anderson_impurity_1D" directory for a detailed description of the code parameters, variables, subroutines, functions, and output files.

### Description of the important parameters and variables (in order of appearance)

#### In main program

- LFIN:              defines the total number of sites in the chain (LFIN-1 conduction sites connected to 1 impurity);
- IZI:               defines on which site the Green's functions will be calculated (IZI=LFIN means on the last site, which is the impurity site);
- IDOWN:             defines if the down spin Green's functions should be calculated (IDOWN=0 means up spin only; IDOWN=1 means both spins);
- PRINTTRANSITIONS:  defines if the atomic transitions should be printed on the output file (PRINTTRANSITIONS=0 means no; PRINTTRANSITIONS=1 means yes);
- ALLBLOCKS:         defines if all blocks of the Hamiltonian matrix should be included in the calculations (ALLBLOCKS=0 only considers blocks where the number of electrons is equal to LFIN, LFIN-1 or LFIN+1; ALLBLOCKS=1 considers all blocks);
- ACOPtI:            electron hopping term between sites of the chain;
- D:                 bandwidth for the zero order density of states;
- ACOPparam:         Anderson parameter DELTA (the Anderson parameter defines the energy scale);
- ACOPhib:           hybridization V (the hybridization between the impurity site and the first conduction site);
- BETAI:             BETA=1/T;
- AMUI:              chemical potential MU;
- ACOPUI:            starting electronic correlation U_i for the correlation loop;
- ACOPUF:            ending electronic correlation U_f for the correlation loop;
- ACOPnI:            local energy E_q of the conduction electrons;
- EFI:               starting impurity gate energy;
- EFF:               ending impurity gate energy;
- MT:                number of points for the impurity gate energy loop;
- LT:                number of points for the temperature loop;
- TI:                starting temperature;
- TF:                ending temperature;
- NT:                number of points for the electron correlation loop;
- UI:                starting electron correlation;
- UF:                ending electron correlation;
- TKHaldane:         Kondo temperature T_K as calculated by Haldane;
- A:                 starting frequency for the occupation numbers;
- B:                 ending frequency for the occupation numbers;
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
- FRIED:             Friedel sum rule.
  
#### In subroutine OCUP

- NT:                number of points for the calculation;
- ETTA:              imaginary part of the complex frequency;
  
#### In subroutine DENSI2

- NT:                number of points for the calculation;
- ETTA:              imaginary part of the complex frequency;
- OI:                starting frequency for the DOS;
- OF:                ending frequency for the DOS.

### Subroutines and functions (in order of appearance)

- PROGRAM AMHUBBARD: main routine, defines the parameters of the model and calls other subroutines;
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
- TCUPSIS:           calculates the up spin creation operator of the whole system;
- TCDOWNSIS:         calculates the down spin creation operator of the whole system;
- THSIS:             calculates the non-local elements of the Hamiltonian matrix, adds them to the local elements and diagonalizes the Hamiltonian matrix;
- CONECTIONS:        defines which Hamiltonian blocks connect to each other and to creation operators ;
- ESTRUCTURA:        defines the size, charge, and spin of Hamiltonian blocks;
- DEALLOC_ESP1:      deallocates memory used to store information about the Fock space of the current iteration;
- MULTI:             calculates the local elements of the Hamiltonian matrix;
- VALIN:             defines the operators and Hamiltonian for a single site.

### Output files

- CGFM_Anderson_impurity_1D.out

    Main output file. The file contains information about the execution of the program, such as: input parameters; Hamiltonian blocks with their respective charge, spin, and eigenvalues; lowest eigenvalue; partition function; and others.
  
- YYYYMMDD_aocup_LFIN=X.dat
  
    Output file for the occupation numbers. The name contains the date the program was run (in YYYYMMDD format) and the number of sites in the chain. The file is structured in the following way:

  1. the temperature in column 1;
  2. the impurity gate energy divided by the Anderson parameter in column 2;
  3. the chemical potential in column 3;
  4. the electron correlation in column 4;  
  5. the Friedel sum rule in column 5;  
  6. the total vacuum occupation number in column 6;  
  7. the total up spin occupation number in column 7;  
  8. the total down spin occupation number in column 8;  
  9. the total double spin occupation number in column 9;  
  10. the total completeness number in column 10;
  11. the total occupation number in column 11.

- YYYYMMDD_adensity_LFIN=X.dat

    Output file for the density of states. The name contains the date the program was run (in YYYYMMDD format) and the number of sites in the chain. The file is structured in the following way:

  1. the frequency in column 1;
  2. the up spin density of states in column 2;
  3. the down spin density of states in column 3.
