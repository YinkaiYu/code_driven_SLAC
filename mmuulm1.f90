SUBROUTINE MMUULM1(A,NF,NTAU,NFLAG)

   Use Blockc
   Implicit Real (KIND=8) (A-G,O-Z)
   Implicit Integer (H-N)

   ! Arguments:
   COMPLEX (Kind=8), Dimension(:,:) ::  A
   Integer :: NF, NTAU,  NFLAG

   ! Local
   COMPLEX (Kind=8), Dimension(:), Allocatable ::   V1, V2
   COMPLEX (Kind=8)  :: UT(2,2), U(2,2)

   N = Size(A,1)
   Allocate (V1(N),V2(N))

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
               A(J,I1) = A(J,I1) / XSIGMA_xL(NSIGL_U(I,NTAU)) ! i up
               A(J,I2) = A(J,I1) / XSIGMA_xR(NSIGL_U(I,NTAU)) ! i down

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
               V1(J) = A(J,I1) * U(1,1) + A(J,I2) * U(2,1)  ! i left
               V2(J) = A(J,I1) * U(1,2) + A(J,I2) * U(2,2)  ! i right

               ! U-term propagator in x eigen-basis
               A(J,I1) = V1(J) / XSIGMA_xL(NSIGL_U(I,NTAU)) ! i left
               A(J,I2) = V2(J) / XSIGMA_xR(NSIGL_U(I,NTAU)) ! i right

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
               V1(J) = A(J,I1) * UT(1,1) +  A(J,I2) * UT(2,1) ! i up
               V2(J) = A(J,I1) * UT(1,2) +  A(J,I2) * UT(2,2) ! i down
               A(J,I1) = V1(J) ! i up
               A(J,I2) = V2(J) ! i down
               
            ENDDO
         ENDDO
      ENDDO
   ENDIF
   
   IF (NFLAG.EQ.2 .AND. RJ.GT.ZERO) THEN
      DO I = 1,LFAM
         I1 = L_Bonds(I,0  )
         I2 = L_Bonds(I,nf1)
         DO J = 1,N
            V1(J)   =  A(J,I1) * U(1,1) + A(J,I2) * U(2,1) 
            V2(J)   =  A(J,I1) * U(1,2) + A(J,I2) * U(2,2) 
         ENDDO
         !          Kenitic
         DO J = 1,N
            A(J,I1) = V1(J) / DCMPLX(XSIGP2(NSIGL_K(I,Nf1,NTAU)),0.D0)
            A(J,I2) = V2(J) / DCMPLX(XSIGM2(NSIGL_K(I,Nf1,NTAU)),0.D0)
         ENDDO
      ENDDO
   ENDIF
   

   IF (NFLAG.EQ.1 .AND. RJ.GT. ZERO) THEN
      DO I = 1,LFAM
         I1 = L_Bonds(I,0   ) 
         I2 = L_Bonds(I,Nf1 ) 
         DO J = 1,N
            V1(J) = A(J,I1) * UT(1,1) +  A(J,I2) * UT(2,1)
            V2(J) = A(J,I1) * UT(1,2) +  A(J,I2) * UT(2,2)
         ENDDO
         DO J = 1,N
            A(J,I1) = V1(J)
            A(J,I2) = V2(J)
         ENDDO
      ENDDO
   ENDIF
   
   Deallocate (V1,V2)

RETURN
END SUBROUTINE MMUULM1


