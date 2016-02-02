!    BROMOC  -  CG-GCMC-BD
!    Electrodiffusion, Gran Canonical Monte Carlo, Brownian,Dynamics 
!    and Coarse Grain Model DNA Simulation Program.
!    Copyright (C) 2014 Pablo M. De Biase (pablodebiase@gmail.com)
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.

! New variables related to DNA insertion in GCMC/BD code
! ======================================================
!
! Strands, nucleotides and interaction sites
! ------------------------------------------
!
! istrs                 -> Number of strands
! inuc                  -> Number of nucleotides in each strand
! cgnuc                 -> Nucleotide charge
! diffnuc               -> Nucleotide diffusion constant
! cdnuc                 -> Effective dielectric constant for DNA
! epsnuc                -> Nucleotide Lennard-Jones potential parameter 
! Qnucl                 -> Logical variable which indicates if nucleotide order
!                          have been called
! Qdie                  -> Logical variable which indicates if effective dielectric 
!                          constant for DNA have been called
! Qsolv                 -> Logical variable which indicates if solvent-induced 
!                          contributions have been called
! Qpar                  -> Logical variable which indicates if particles order
!                          have been called
! Qsystem               -> Logical variable which indicates if system order
!                          have been called
! Qbuf                  -> Logical variable which indicates if buffers order
!                          have been called
! Qtraj                 -> Logical variable which indicates if DNA and/or ions 
!                          trajectories have to be written in a outputfile
! Qfmemb                -> Logical variable which indicates if traslocation duration 
!                          for cylindrical pores has to be calculated
! Qstfx                 -> Logical variable which indicates if DNA fixed sites are 
!                          turned on
! eps(P,S,Ab,Tb,Cb,Gb)  -> Lennard-Jones parameter for interaction between ions 
!                          and sites (combination rules)
! sg(P,S,Ab,Tb,Cb,Gb)   -> Lennard-Jones parameter for interaction between ions
!                          and sites (combination rules)
! epsLJ(dtype,dtype)    -> Lennard-Jones parameter for interaction between ions
!                          and sites
! sgLJ(dtype,dtype)     -> Lennard-Jones parameter for interaction between ions
!                          and sites
! epsolv                -> Energy scale for solvnet-induced contribution
! maxsite               -> Maximum number of interaction sites
! mmxsites              -> mmxsites = maxsite*(maxsite-1)/2 + maxsite - 1
! nsites                -> Number of interaction sites
! nstfx                 -> Number of DNA fixed sites
! strand(maxsite)       -> This integer vector indicates the strand for the 
!                          interaction site
! typenuc(maxsite)       -> This integer vector indicates the nucleotide type 
!                          for the interaction site
! namnucl(maxsite)      -> This character vector store the name nucleotides
!                          (A,T,C,G) 
! namsite(maxsite)      -> This character vector stores the name interaction 
!                          sites P=Phosphate, S=Sugar, (Ab,Tb,Cb,Gb)=Bases
! stfx(maxsite)         -> This integer vector stores DNA fixed sites
! stfree(maxsite)       -> This logical*1 vector indicates if an site is not fixed
! (x,y,z,r,phi)nat(maxsite) -> Interaction site coordinates for the native 
!                              structure
! Qtras                 -> Logical variable which indicates if there is a traslation of DNA 
!                          sites positions
! Qrot                  -> Logical variable which indicates if there is a rotation of DNA
!                          sites positions
! typat(datom), typtyp(datom), nwtype(dtype)              
!                        -> These integer vectors are used in wrtrraj routine
! nold                  -> number of sites types
! ndna                  -> number of dna fragment types 

module nucleotmod
implicit none
integer     istrs, inuc, maxsite, mmxsites, nsites, nstfx,extraP,nsites1st
!parameter   (maxsite  = 1000)
!parameter   (mmxsites = maxsite*(maxsite-1)/2 + maxsite - 1)
integer,allocatable     ::  strand(:), typenuc(:), stfx(:)
logical*1,allocatable   ::  stfree(:)
character*1,allocatable :: namnucl(:)
character*2,allocatable :: namsite(:)
real,allocatable ::   xnat(:), ynat(:), znat(:), rnat(:), phinat(:)
real        cgnuc, diffnuc, epsnuc, cdnuc, epsolv, fctn, scalepairing
real        notrx,notry,notrz,insites  ! NOTRAN
integer     ndna,nion,nttyp,ndnaxnion,setframes
!CONTRA
real        xcon,ycon,zcon
real,allocatable :: kx(:),ky(:),kz(:),contrx(:),contry(:),contrz(:) 
integer     ctn
integer,allocatable ::  csn(:)
logical*1   Qcontrans,Qcontprint,Qunsplit
logical*1   Qnucl, Qpar, Qsystem, Qbuf, Qtraj, Qtrajcont, Qdie, Qsolv, Qfmemb, Qstfx
logical*1   Qtras, Qrot, Qnotrans, Qnotrx, Qnotry, Qnotrz, Qdnafree, Qinvstr, QfirstP
real cylall(6,3),din,ain

! Bonded terms
! ------------
!
! maxbond               -> Maximum number of bonds
! nbond                 -> Number of bonds
! sitebond(maxbond,2)   -> This integer matrix indicates the bonded 
!                          interaction sites
! distbond(maxbond)     -> This real vector stores the equilibrium 
!                          bond lenghts
! bond(mmxsites)        -> Logical supervector which indicates if two sites 
!                          are forming a bond
! typbond(maxbond)      -> This integer vector indicates if the bond 
!                          is either intranucleotide (0) or 
!                          internucleotide (1)
!
! maxang                -> Maximum number of bond angles
! nangle                -> Number of bond angles
! siteangle(maxang,3)   -> This integer matrix indicates the interaction 
!                          sites which are forming a bond angle
! valangle(maxang)      -> This real vector stores the equilibrium bond
!                          angles
! angle(mmxsites)       -> Logical supervector which indicates if two sites 
!                          are forming abond angle
!
! maxdihe               -> Maximum number of dihedral angles
! ndihe                 -> Number of dihedral angles
! sitedihe(maxdihe,4)   -> This integer matrix indicates the interaction
!                          sites which are forming a dihedral angle
! valdihe(maxdihe)      -> This real vector stores the equilibrium dihedral
!                          angles
! (dSAb,dSTb,dSCb,dSGb,dPS5,dPS3) -> Natural bond lenghts
! (phPSAb,phPSTb,phPSCb,phPSGb,phPSAb2,phPSTb2,phPSCb2,phPSGb2,phSPS,phPSP) ->
!                                    Natural bond angles
! (dhAbSPS,dhTbSPS,dhCbSPS,dhGbSPS,dhSPSAb,dhSPSTb,dhSPSCb,dhSPSGb,dhSPSP,dhPSPS) -> 
!                                    Natural dihedral angles

integer   maxbond, maxang, maxdihe  
!parameter (maxbond = maxsite-1)
!parameter (maxang  = maxsite-1)
!parameter (maxdihe = maxsite-2)
integer   nbond, nangle, ndihe
integer,allocatable ::   sitebond(:,:),siteangle(:,:),sitedihe(:,:)
real,allocatable ::      distbond(:),valangle(:),valdihe(:)
logical*1,allocatable :: bond(:),angle(:)
integer,allocatable ::   typbond(:)
real      dSAb, dSTb, dSCb, dSGb, dPS5, dPS3
real      phPSAb, phPSTb, phPSCb, phPSGb, phPSAb2, phPSTb2, phPSCb2, phPSGb2, phSPS, phPSP
real      dhAbSPS, dhTbSPS, dhCbSPS, dhGbSPS, dhSPSAb, dhSPSTb, dhSPSCb, dhSPSGb, dhSPSP, dhPSPS 

! Non-bonded terms
! ----------------
!
! maxstack              -> Maximum number of native contacts 
! maxbp                 -> Maximum number of hydrogen bondings
! maxex                 -> Maximum number of excluded volume interactions
! maxqq                 -> Maximum number of Coulomb interactions
! maxsolv               -> Maximum number of solvent-induced contributions
! nstack                -> Number of native contacts
! nbp                   -> Number of hydrogen bondings
! nex                   -> Number of excluded volume interactions
! nqq                   -> Number of Coulomb interactions
! nsolv                 -> Number of solvent-induced contributions
! sitestack(maxstack,2) -> This integer matrix indicates the interaction
!                          sites which are forming a native contact
! sitebp(maxbp,2)       -> This integer matrix indicates the interaction
!                          sites which are forming a hydrogen bonding
! siteex(maxex,2)       -> This integer matrix indicates the interaction
!                          sites which are forming a excluded volume interaction
! siteqq(maxqq,2)       -> This integer matrix indicates the interaction
!                          sites which are forming a Coulomb interaction
! siteslv(maxsolv,2)    -> This integer matrix indicates the interaction
!                          sites which are forming a solvent-induced contributions
! sgstack(maxstack)     -> sigma parameter in intra-strand native contacts term
! sgbp(maxbp)           -> sigma parameter in hydrogen bonding terms
! sgex(maxex)           -> sigma parameter in excluded volume

integer           maxstack, maxbp, maxex, maxqq, maxsolv
!parameter         (maxstack = 1000000)
!parameter         (maxbp    = 1000000)
!parameter         (maxex    = 1000000)
!parameter         (maxqq    = 1000000)
!parameter         (maxsolv  = 1000000)
integer           nstack, nbp, nex, nqq, nsolv
integer,allocatable ::           sitestack(:,:),sitebp(:,:),siteex(:,:),siteqq(:,:),siteslv(:,:)
real,allocatable ::              sgstack(:),sgbp(:),sgex(:)

! Potential energy contributions
! ------------------------------
! ebond    -> strech energy
! eang     -> bending energy
! edihe    -> torsional energy
! estack   -> native contacts
! ebp      -> hydrogen bonding
! eex      -> excluded volume
! eqq      -> coulomb interaction
! esolv    -> solvent-induced contribution
! eqqmx    -> coulomb interaction between ions and sites
! evdwmx   -> LJ potential between ions and sites
! Qdeby    -> Logical variable which indicates if Coulomb
!             interactions are taken into account using the
!             Debye-Hückel approximation 
! Qionsite -> Logical variable which indicates if there are 
!             ions-sites interactions using combination rules 
!             for obtaining LJ parameters
! Qljpar   -> Logical variable which indicates if there are
!             ions-sites interactions using LJ parameters
!             directly (i.e, don't using combination rules)
! ionstr   -> Ionic strength [Mol/L]

real ebond, eang, edihe
real estack, ebp, eex, eqq, esolv
real eqqmx, evdwmx
logical*1 Qdeby, Qljsin, Qljpar, Qninfo, Qdebyhyb
real  ionstr

! Fraction of denatured bases
! ---------------------------

! Qfbases -> Logical variable which indicates if fraction of denatured bases will be
!            calculated after simulation of DNA
! fbases  -> Fraction of denatured bases

real fbases
logical*1 Qfbases

end module 
