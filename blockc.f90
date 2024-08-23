      Module  Blockc 
        
!!        Use Lattices_v3


        Real (Kind=8),save :: BETA , RJ, RHUB, RHUB1, RHUB2, RT1, DTAU, PI,  TwistX
        Integer,      save :: LTROT, LTROT_quench, LTROT_ST, LTROT_EN, LTROT_ME, NWRAP, N_SUN, NBIN, NSWEEP,  NORB, NORB1, NDIM, LTAU,  &
             &                LFAM, NFAM, Nbond, nspin, LQ, NLX, NLY, NE, Itwist, CH
        real(kind=8), DIMENSION(:), save :: a1_p(2), a2_p(2), b1_p(2), b2_p(2)

        Integer,  Dimension(:,:), Allocatable, Save :: list(:,:), invlist(:,:), nlist(:,:), invnlist(:,:,:,:)
        Integer,  dimension(:,:), allocatable, save :: lattimj(:,:)
	    INTEGER,  Dimension(:,:), Allocatable, Save :: L_bonds, NSIGL_U
        INTEGER,  Dimension(:,:,:), Allocatable, Save :: NSIGL_K
        INTEGER,  Save ::  NFLIPL(-1:1)
        Complex (kind=8), Dimension(:,:), Allocatable, Save :: GRUP(:,:), GRUPC(:,:)

        Complex(Kind=8), save :: phase
        
        COMPLEX (Kind=8), Dimension(:,:),   Allocatable, Save ::  PROJ, ZKRON
        COMPLEX (Kind=8), Dimension(:,:)  , Allocatable, Save ::  URT_tot, URTM1_tot

        COMPLEX (Kind=8),  Save ::  DELTA_xL(-1:1), DELTA_xR(-1:1), XSIGMA_xL(-1:1), XSIGMA_xR(-1:1), &
             &                       UR_K(2,2), URT_K(2,2), obs(30)

        REAL (Kind=8), Save ::   XSIGP2(-2:2),XSIGM2(-2:2),ETAL(-2:2), &
             &                   GAML(-2:2),   DELLP2(-2:2,3), DELLM2(-2:2,3), DGAML(-2:2,3), &
             &                   FD(4), ZERO, IsingS(-1:1)

        Logical :: L_Trot_hop

        Contains

        Subroutine Allocate_Blockc

          Implicit none

          Zero    =  1.0D-10
          PI      =  acos(-1.d0)

          a1_p(1) = 1.d0
          a1_p(2) = 0.d0
          a2_p(1) = 0.d0
          a2_p(2) = 1.d0
          b1_p(1) = 2.d0*Pi
          b1_p(2) = 0.d0
          b2_p(1) = 0.d0
          b2_p(2) = 2.d0*Pi
          
          NORB = 1 
          NORB1 = 1
          LQ   = NLX * NLY  
          NFAM = 4                 ! no use in single Dirac cone model
          LFAM = LQ
          Nbond = NFam
          nspin = 2                ! **IMPORTANT CHANGE** for up-down hopping
          NDIM = Norb*LQ * nspin   ! **IMPORTANT CHANGE** for spin degree of freedom
          DTAU   = BETA/DBLE(LTROT)
          LTROT_ME = LTROT/2
          LTROT_ST = LTROT_ME - LTROT_quench; LTROT_EN = LTROT_ME + LTROT_quench + 1
          
          
          Allocate( L_Bonds(LQ,0:Nbond), list(LQ,2), invlist(NLX,NLY), nlist(Ndim, 4), invnlist(NLX,NLY,norb,nspin), &
               &    NSIGL_K(LQ,Nfam,LTROT), NSIGL_U(LQ,LTROT), Lattimj(LQ,LQ)  )
          Allocate ( PROJ(NDIM,NDIM),   ZKRON(NDIM,NDIM) ) 
          Allocate ( URT_tot(NDIM,NDIM), URTM1_tot(NDIM,NDIM) )
          Allocate ( GRUP(Ndim,Ndim),GRUPC(Ndim,Ndim))

        end Subroutine Allocate_Blockc

      end Module Blockc
