# Coupled Channels Solver

This suite of MATLAB code is designed for solving coupled channels equations that arise when considering either the scattering or bound state properties of ultracold alkali metal atoms in magnetic fields.  Code for handling scattering problems can be found in the folder 'scattering-solver', and code for handling bound state problems can be found in the folder 'bound-solver'.  Both of these use files from the folder 'common'.

Both the scattering and bound state solvers integrate the set of coupled radial Schrodinger equations for the set of _channels_ that are coupled together by the interaction potential.  For alkali metal atoms, and when neglecting the magnetic dipole-dipole interaction, these channels are defined by states with the same total orbital angular momentum L, orbital projection mL, and total spin projection mF which is the sum of each atoms' electronic spin projection mS_i and nuclear spin projection mI_i.  When the magnetic dipole-dipole interaction is included the coupled channels are those with the same mL+mF and whose spatial states have the same parity -- i.e., the values of L differ by 2.  So this means that s states (L=0) couple to d states (L=2) and also to g states (L=4) and so on.

Both solvers use the same method for integrating the coupled radial Schrodinger equations: namely, they use the improved log-derivative method of D. E. Manolopoulos to integrate the equivalent matrix Riccati equation using a form of invariant embedding.  A grid is formed where the step size is the shortest WKB wavelength for the system multiplied by a user-settable scale factor.  Integrate proceeds within user-settable ranges rmin to rmax.  For scattering problems integration proceeds only in the forward direction (rmin to rmax), whereas for bound state problems the integration has to proceed in both the forward and backward directions in order to satisfy both the r=0 and r=Inf boundary conditions.  Although the integration method is the same, the scattering and bound state integrators, labelled 'manolopoulos_scatt' and 'manolopolous_bound', are different functions because the bound state version needs to compute more quantities in order to be useful for bracketing energy ranges.

Currently, the coupled channels solver handles all isotopes and combinations of Rb and K.

In order to use either the scattering or bound state solvers, you will need to first add the 'common' folder to your MATLAB path.

## Basis states

The first step to using either the scattering or bound state solver is defining the basis states and the transformations between different bases.  These basis vectors (bv) and basis transformations (bt) are handled by the 'atompairbasis' class.  Different bases are better than others at defining various terms in the total Hamiltonian.  The Zeeman Hamiltonian is most conveniently written in terms of a fully uncoupled basis where all spins are kept separate.  This is 'basis 1' or 'bv1', and it is characterized by a set of 6 numbers |L mL mS1 mI1 mS2 mI2>.

The main interaction due to the Coulomb interaction of the various clouds of electrons gives rise to what are known as the Born-Oppenheimer (BO) potentials.  For alkali metal atoms, there are two potentials which depend on the total electron spin of the two valence electrons: either spin 0 (singlet) or spin 1 (triplet).  The singlet potential always has a deeper potential well because the anti-symmetric spin state allows for a symmetric spatial state that leads to a stronger covalent bond.  So when defining the BO contribution to the total Hamiltonian, it makes sense to break it into operators that differentiate based on the total electronic spin.  This basis, where S1+S2 = S but the nuclear spins are kept uncoupled, is known as 'basis 2' or 'bv2' and is characterized by a set of 6 numbers |L mL S mS mI1 mI2>

The hyperfine Hamiltonian is most conveniently written in terms of the coupled electronic and nuclear spins F_i = S_i+I_i for each atom individually.  This is known as 'basis 3' or 'bv3' and is characterized by a set of 6 numbers |L mL F1 mF1 F2 mF2>.

Finally, we need a set of basis vectors that are eigenstates of the single-atom internal Hamiltonians (hyperfine + Zeeman) which is denoted H_int.  The total potential approaches H_int as r -> Inf.  We call this basis the _internal_ basis, and for scattering problems the asymptotic states are always internal states.  These states are labelled by whole numbers starting at one in order of ascending energy.  So for Rb-87, the state |F=1,mF=1> is labelled |1> as it is the lowest energy state for a single atom, and the state |F=1,mF=0> is labelled |2>, and so on.  A pair of Rb-87 atoms in the |F1=1,mF1=1>|F2=1,mF2=0> state is labelled as |1,2>.  Note that since these are identical particles the states |1,2> and |2,1> are the same.  For KRb atom pairs, the K label is stated first so that |1,2> = |FK=9/2,mFK=-9/2>|FRb=1,mFRb=0>.  The complete basis vector is labelled by 6 numbers |L mL mF1 mF2 internal1 internal2>.  

Transformation matrices that convert vectors or matrices in one basis to another are written as btab where 'b' is the starting basis and 'a' is the final basis.  For instance, bt21 takes a vector written in basis 1 and transforms it into a vector written in basis 2.  The basis transformations bt21, bt31, and bt32 are independent of magnetic field, whereas the basis transformation bt2int is dependent on the magnetic field.

Before running either solver, you must create the basis states.  This can be done directly by using 

    basis = atompairbasis(AtomPair,LVec);

where 'AtomPair' is a string denoting which pair of atoms to consider, and 'LVec' is a vector of L values to use.  Valid 'AtomPair' values are:

* 'KRb' or 'K40Rb87'
* 'RbRb' or 'Rb87Rb87'
* 'KK' or 'K40K40'
* 'K41K41'
* 'K40K41'
* 'K39K39'
* 'K39K40'
* 'K39K41'
* 'K41Rb87'
* 'K39Rb87'
* 'K39Rb85'
* 'K40Rb85'
* 'K41Rb85'
* 'Rb85Rb85'
* 'Rb85Rb87'

Other combinations of the various isotopes can be added easily in the class defintion for 'atompairbasis'

You can also use the wrapper file 'BasisSetGeneration' which will automatically save the created basis as a .mat file in your current directory.  The usage is

    basis = BasisSetGeneration(AtomPair,LVec,filename);

where 'filename' is optional.  If it is omitted a filename is automatically generated.

## Scattering solver

Scattering problems are those where E - E_i > 0 for at least one channel in the system for a given energy E and channel threshold E_i which itself depends on the applied magnetic field B.  Atoms start in well-defined initial states in the internal basis.  In a scattering problem one needs to specify the entrance channel, and this is done in the internal basis.  You can calculate the scattering matrix for a particular energy E and magnetic field B using the function

    S = MultiChannelScattering(initLabel,E,B,outputFile,basis,opt);

Here, E is specified in uK and B in Gauss - these are converted internally into units appropriate for the solver.  The variable 'initLabel' is a 4 element vector '[L mL internal1 internal2]' (an abridged version of bvint) and this specifies the entrance channel.  'outputFile' is the file to save all the variables to at the end of the calculation - if it is set to an empty string then nothing is saved.  'basis' is either an instance of the class 'atompairbasis' or it is a path pointing towards a previously saved basis file.  Finally, 'opt' is an instance of the class 'scattoptions' which provides a number of options related to the integration of the coupled channels equations.  These are:

* rmin: minimum integration distance
* rmax: maximum integration distance
* drscale: the scaling factor between the shortest WKB wavelength and the step size
* drmax: maximum step size
* drmin: minimum step size
* parallel: set to true to use a parfor loop, false to use a normal loop
* getwf: set to true to calculate wavefunctions in addition to scattering matrix
* dipole: set to true to use dipole-dipole interaction - this takes a lot of memory and can be very slow
* blocksize: size of blocks (in Angstroms) for which the step size is constant.  The potential operator is calculated and store in blocks of this size, so if it is too large you will run out of memory.  You should rarely need to change this.

The options class can be instantiated using

    opt = scattoptions(property1,value1,property2,value2,...);

where 'property' is a character string indicating the property to set and 'value' is the value to set.  So

    opt = scattoptions('parallel',false,'getwf',true);

creates a 'scattoptions' object with default properties except that 'parallel' is set to false and 'getwf' is set to true.  If you are calculating the wavefunctions you can get them from the MultiChannelScattering function as the second output argument.

### Scattering Matrix

The scattering matrix S that is the output of MultiChannelScattering(...) is an instance of the class ScatteringMatrix.  It contains the full S and T matrices for the system, along with the appropriate basis vector labels.  Note that when considering identical particles the internal basis states are reduced to unique labels only -- |1,2> is the same as |2,1>, so only |1,2> is retained.  The scattering matrix is transformed appropriately to take into account the necessary symmetries.  The ScatteringMatrix class mostly exists to provide a convenient method by which one can access elements of the S or T matrices without having to muck around with an Nchannels x Nchannels x Nruns array.  The special aspect of the class is that its 'subsref' function is overloaded so that different elements of the S or T matrices can be accessed in a more convenient fashion.

General usage for an instance S of ScatteringMatrix is S(subscripts).(Property).  If (Property) is something other than S or T, the default behaviour applies, so one can access bv or Sfull or Tfull in a sensible manner.

If one uses S.S or S.T, it returns the elastic element of the S or T matrix, which is entrance channel to entrance channel.
            
If one uses S(idx).S or S(idx).T then the element that is accessed is actually S(entrance channel index,idx)

If one uses S(L,mL,int1,int2).S then the element that is accessed is entrance channel to the channel specified by the 4-element vector '[L,mL,int1,int2]'.

If one uses S(idx1,idx2).S then the element that is accessed is the one going from idx1 to idx2 where these are indices in the basis vector bv

If one uses S(L1,mL1,int1_1,int2_1,L2,mL2,int1_2,int2_2).S then the element that is accessed is the one starting from the 4-element vector '[L1,mL1,int1_1,int2_1]' and ending with '[L2,mL2,int1_2,int2_2]'

### Usage
Below is a set of MATLAB commands that can be used to calculate and plot the T-matrix for two Rb-87 atoms each in the |F=1,mF=1> state near a Feshbach resonance at 1007 G.  First, navigate to the folder 'scattering-solver' and run:

    addpath('../common');
    opt = scattoptions('parallel',true);
    basis = atompairbasis('RbRb',0:2);
    E = 1;  %In uK
    B = linspace(1006,1008,1e2);    %In G
    S = MultiChannelScattering([0,0,1,1],E,B,'',basis,opt);
    plot(S.B,abs(S.T).^2,'.-');

The resulting plot should show a Beutler-Fano profile close to 1007 G.

## Bound state solver

Bound state problems are those where E < E_i for all relevant channels so that only bound states can exist.  From a mathematical point of view, bound states occur at energies where the radial wavefunction is normalizable, which means that it has the boundary conditions of being zero at the origin and decay to zero as r -> Inf sufficiently fast to retain a finite integral.  Unlike scattering problems, there is no question of an _entrance channel_, and we do not really care as much about the internal basis as they become relevant only at infinite separations.  Instead, we care more about basis 2 which is where the electronic spin is coupled together and is the basis in which the BO potentials are diagonal.

Since solutions to bound state problems occur only at specific energies, rather than having a solution at _every_ energy, finding bound states is more computationally intensive than solving for a scattering matrix.  The procedure that is used here is:

1. Define the basis states and coupled channels via '[L, mL, mF]'.
2. Choose an energy range (E <= 0) in which to search for bound states.  The larger the energy range, the longer this takes.
3. Use node counting to break the initial energy range into subranges containing exactly one bound state.
4. Within those subranges, integrate the equations in the forward and backward directions, subtract the resulting log-derivative matrices at a matching point, and find the energies at which that difference matrix has a zero eigenvalue.  These are bound state energies.
5. Optionally calculate the bound state wavefunctions

### Basis states and coupled channels

The bound state solver uses the same class, 'atompairbasis', as the scattering solver for defining the basis vectors and transformations.  A major difference is encountered in specifying the 'initial state'.  In scattering calculations, this is the entrance channel to which the energy is referenced and it therefore makes sense to specify it in terms of the internal states.  In bound state calculations, where atoms do not fly out to infinite separations, it makes more sense to specify the 'initial state' in terms of its angular momementum and total spin projection mF.  In this sense it is no longer an 'initial state' but rather a set of conserved (or almost conserved) quantum numbers that determine the channels which are coupled.  Bound state energies are referenced to the lowest channel threshold among the coupled channels.

### Node counting

In order to maintain orthogonality, each bound state for a given set of coupled channels has a different number of nodes than the others.  Considering that the ground state for any quasi-harmonic oscillator has exactly zero nodes, it follows that each additional bound state has one additional node in the wavefunction.  We can use this feature to determine the number of bound states in any given energy range and then find subranges where there is exactly one bound state in each energy range.  These subranges are determined by recursive bisection of the ranges using the function 'findranges'.  Once a subrange with one bound state has been identified the subrange is narrowed by further iterations of bisection to improve convergence during integration - the number of iterations can be set using the 'boundoptions' class.

In the one dimensional case node counting is fairly simple.  One integrates the radial Schrodinger equation from rmin to rmax (where rmax is far, far into the classically forbidden region) and counts the number of nodes that one encounters.  When using log-derivative solution methods, one counts the number of poles that are encountered.  Multichannel node counting is somewhat harder, but in the end what one does is counts the number of times the solution matrix has a zero eigenvalue (for wavefunction methods).  A review paper by Jeremy Hutson (Computer Physics Communications 84 (1994) 1-18) makes a number of comments about counting the number of negative eigenvalues for log-derivative methods, but this is difficult to implement in practice when using MATLAB's eigenvalue function 'eig'.  Instead, I have found that it is necessary to adiabatically follow eigenvalues of the log-derivative matrix, by comparing old eigenvectors to current ones and matching eigenvalues based on the inner products, and then counting the number of poles encountered for each adiabatically followed eigenvalue.  This is slightly more computationally expensive but tends to be more robust than the method proposed by Hutson.


### Integration and convergence

For finding a bound state of the single channel Schrodinger equation, one needs to integrate the equation in the forward direction from rmin to a matching point somewhere in the classically allowed region and then integrate backwards from rmax to the matching point.  The condition for a bound state to exist at that energy is that the log-derivatives at the matching point must be the same.  For multi-channel problems, the condition is that the difference in the two log-derivative matrices must have a zero eigenvalue.  Hutson has some comments on how to pick this matching point, which apparently come from Manolopoulos's PhD thesis, based on finding the first node in the log-derivative when integrating in the forwards direction, but these don't work for alkali metal atoms (and probably not for other systems) because the matching point has to be in the classically allowed region for all channels that contribute.  For alkali metal atoms, the singlet and triplet potentials have very different repulsive cores and the classically allowed region for the singlet potential is much larger than for the triplet potential.  Using this method means that any solution that has any contribution from the singlet potential (i.e. all of them) will be excluded.

For bound states near threshold, one can simply use the equilibrium position of the triplet potential (about 6 Angstroms) as the matching point.  This will become a problem for very deeply bound states in the singlet potential because it will be in the classically forbidden region, but these situations are rare enough to be identifiable by the user.

The method for converging on an energy is the secant method with a slight modification to ensure that each iteration of the energy is within the starting energy range.  If the solver fails to converge within a settable number of iterations, then MATLAB's 'fsolve()' function is used, but be warned that it does not respect limits.  You may find that it gives unphysical (positive) energies or misses bound states altogether.

### Usage

Below is a set of MATLAB commands that can be used to calculate and plot the bound state energies for two Rb-87 atoms with total spin projection of 2 near a Feshbach resonance at 1007 G.  First, navigate to the folder 'bound-solver' and run:

    addpath('../common');
    opt = boundoptions();
    basis = atompairbasis('RbRb',0:2);
    E = [0,-1500];   %In uK
    B = linspace(1006,1008,25);    %In G
    [Eout,results] = MultiChannelBound([0,0,2],E,B,'',basis,opt);
    plot(B,Eout,'o-');

This will solve for and plot the bound state energies as a function of magnetic field.  For B < 1007 G there should be two bound states and for B > 1007 G there should be one bound state.  The disappearance of a bound state as it moves through the channel threshold creates a Feshbach resonance.

The class 'boundoptions' contains a number of properties related to the integration of the coupled channels equations.  Most of them will not need to be changed by the user, but ones that are of particular interest are:

* iter: number of iterations for converging on a bound state energy
* rangeiter: number of iterations for narrowing energy ranges using node counting
* output: create output wavefunctions
* debug: display debugging information
* stopR: r value to stop at when matching forwards and backwards integration
* parallel: set to true to use a parfor loop and false to use a normal for loop
* dipole: set to true to include the dipole-dipole interaction.  This will be very, very slow

The outputs from MultiChannelBound are a matrix of bound state energies (in uK) which has dimensions numel(B) x (maximum number of bound states), since the number of bound states in a given energy range is dependent on the magnetic field.  Where there are fewer than the maximum number of bound states, the values replaced by NaNs.

'results' is an instance of the class 'boundresults' which simply provides a convenient container for bound state energies at a given magnetic field as well as the wavefunctions (if calculated).