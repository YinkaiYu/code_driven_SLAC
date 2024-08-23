SUBROUTINE MMUUR(A, NF, NTAU, NFLAG)

   !	In A Out U(NF) * A                if NFLAG = 1
   !	In A Out EXP(D(NF)) * UT(NF) * A  if NFLAG = 2
   !	In A Out U(zx) * A                if NFLAG = 3
   !	In A Out EXP(D(xx)) * UT(xz) * A  if NFLAG = 4

   Use Blockc
   Implicit Real (KIND=8) (A-G,O-Z)
   Implicit Integer (H-N)

   COMPLEX (Kind=8), Dimension(:,:) :: A 
   INTEGER :: NF, NT , NFLAG

   !	Local
   COMPLEX (Kind=8),  Dimension(:), Allocatable ::  V1, V2
   COMPLEX (Kind=8) ::  UT(2,2), U(2,2)

   N = SIZE(A,2) ! particle number

   ALLOCATE (V1(N), V2(N))
      
   ! rotation between i up-down and i left-right
   NF1 = NF
   DO I = 1,2
      DO J = 1,2
         U(I,J)  = UR_K (I,J)
         UT(I,J) = URT_K(I,J)
      ENDDO   
   ENDDO

   IF (NFLAG.EQ.5.and.RHUB.GT.zero) THEN
      DO J = 1,N
         DO nx = 1,NLX
            DO ny = 1,NLY

               i  = invlist(nx,ny)
               i1 = INVNLIST(nx,ny,1,1) ! i up
               i2 = INVNLIST(nx,ny,1,2) ! i down

               ! U-term propagator in z eigen-basis
               A(I1,J) = XSIGMA_xL(NSIGL_U(I,NTAU)) * A(I1,J) ! i up
               A(I2,J) = XSIGMA_xR(NSIGL_U(I,NTAU)) * A(I2,J) ! i down

            ENDDO
         ENDDO
      ENDDO
   ENDIF

   IF (NFLAG.EQ.4.and.RHUB.GT.zero) THEN
      DO J = 1,N
         DO nx = 1,NLX
            DO ny = 1,NLY

               i  = invlist(nx,ny)
               i1 = INVNLIST(nx,ny,1,1) ! i up
               i2 = INVNLIST(nx,ny,1,2) ! i down

               ! rotation from i up-down to i left-right
               V1(J) =  UT(1,1) * A(I1,J) + UT(1,2) * A(I2,J) ! i left
               V2(J) =  UT(2,1) * A(I1,J) + UT(2,2) * A(I2,J) ! i right

               ! U-term propagator in x eigen-basis
               A(I1,J) = XSIGMA_xL(NSIGL_U(I,NTAU)) * V1(J) ! i left
               A(I2,J) = XSIGMA_xR(NSIGL_U(I,NTAU)) * V2(J) ! i right

            ENDDO
         ENDDO
      ENDDO
   ENDIF

   IF (NFLAG.EQ.3.and.RHUB.GT.zero) THEN
      DO J = 1,N
         DO nx = 1,NLX
            DO ny = 1,NLY

               i  = invlist(nx,ny)
               i1 = INVNLIST(nx,ny,1,1) ! i up
               i2 = INVNLIST(nx,ny,1,2) ! i down

               ! rotation back from i left-right to i up-down
               V1(J) =  U(1,1) * A(I1,J) + U(1,2) * A(I2,J) ! i up
               V2(J) =  U(2,1) * A(I1,J) + U(2,2) * A(I2,J) ! i down
               A(I1,J) = V1(J) ! i up
               A(I2,J) = V2(J) ! i down
               
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   IF (NFLAG.EQ.2 .AND.  RJ .GT. ZERO ) THEN
      DO I = 1,LFAM
         I1 = L_bonds  (I,0  )
         I2 = L_bonds  (I,nf1)
         DO J = 1,N
            V1(J)   =  UT(1,1) * A(I1,J) + UT(1,2) * A(I2,J)
            V2(J)   =  UT(2,1) * A(I1,J) + UT(2,2) * A(I2,J) 
         ENDDO
            
         !Kinetic
         DO J = 1,N
            A(I1,J) = DCMPLX(XSIGP2(NSIGL_K(I,nf1,NTAU)),0.D0)*V1(J)
            A(I2,J) = DCMPLX(XSIGM2(NSIGL_K(I,nf1,NTAU)),0.D0)*V2(J)   
         ENDDO
      ENDDO
   ENDIF

   IF (NFLAG.EQ.1 .AND. RJ .GT. ZERO) THEN
      DO I = 1,LFAM
         I1 = L_bonds  (I,0  )
         I2 = L_bonds  (I,nf1)
         DO J = 1,N
            V1(J)   =  U(1,1) * A(I1,J) +  U(1,2) * A(I2,J) 
            V2(J)   =  U(2,1) * A(I1,J) +  U(2,2) * A(I2,J) 
         ENDDO
         DO J = 1,N
            A(I1,J) = V1(J)
            A(I2,J) = V2(J)
         ENDDO
      ENDDO
   ENDIF
   
   DEALLOCATE (V1, V2)

   RETURN
END SUBROUTINE MMUUR



