SUBROUTINE SALPH
! Calculate DTAU and ALPHA
! Here so-called ALPHA is just lambda

   Use Blockc
   Implicit Real (KIND=8) (A-G,O-Z)
   Implicit Integer (H-N)

   REAL (Kind=8) :: lambda

   ! set lambda, for Hubbard U
   lambda = acosh(EXP((DTAU*RHUB)/2))

   ! auxiliary Ising field S (field value)
   IsingS(-1) = -1.d0
   IsingS(+1) = +1.d0

   ! flip the index of auxiliary Ising field S (index value)
   NFLIPL(-1) = +1
   NFLIPL(+1) = -1

   ! x-channel rotation between i_up and i_do (on-site)
   DO M = 1,2
      DO N = 1,2
         UR_K(M,N) = DCMPLX(0.D0,0.D0)
      ENDDO
   ENDDO
   UR_K(1,1) =   DCMPLX( 1.D0/SQRT(2.D0),0.D0 )
   UR_K(1,2) =   DCMPLX( 1.D0/SQRT(2.D0),0.D0 )
   UR_K(2,1) =   DCMPLX( 1.D0/SQRT(2.D0),0.D0 )
   UR_K(2,2) =   DCMPLX(-1.D0/SQRT(2.D0),0.D0 )
   DO M = 1,2
      DO N = 1,2
         URT_K(M,N) = DCONJG(UR_K(N,M))
      ENDDO
   ENDDO

   ! U-term propagator in x eigen-basis
   XSIGMA_xL(-1) = EXP( lambda * DCMPLX( +1.d0*IsingS(-1), 0.D0) )
   XSIGMA_xL(+1) = EXP( lambda * DCMPLX( +1.d0*IsingS(+1), 0.D0) )
   XSIGMA_xR(-1) = EXP( lambda * DCMPLX( -1.d0*IsingS(-1), 0.D0) )
   XSIGMA_xR(+1) = EXP( lambda * DCMPLX( -1.d0*IsingS(+1), 0.D0) )

   ! upgrading delta
   DELTA_xL(-1)  = EXP( lambda * DCMPLX( -2.d0*IsingS(-1), 0.D0) ) - DCMPLX(1.D0,0.D0)
   DELTA_xL(+1)  = EXP( lambda * DCMPLX( -2.d0*IsingS(+1), 0.D0) ) - DCMPLX(1.D0,0.D0)
   DELTA_xR(-1)  = EXP( lambda * DCMPLX( +2.d0*IsingS(-1), 0.D0) ) - DCMPLX(1.D0,0.D0)
   DELTA_xR(+1)  = EXP( lambda * DCMPLX( +2.d0*IsingS(+1), 0.D0) ) - DCMPLX(1.D0,0.D0)
   
   RETURN
END SUBROUTINE SALPH
