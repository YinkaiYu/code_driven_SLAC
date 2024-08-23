SUBROUTINE MMUURM1(A, NF, NTAU,  NFLAG)

   Use Blockc
   Implicit Real (KIND=8) (A-G,O-Z)
   Implicit Integer (H-N)
   COMPLEX (KIND=8), Dimension(:,:) :: A
   Integer :: NF, NTAU,  NFLAG

   !Local
   COMPLEX (KIND=8), Dimension(:), Allocatable ::  V1, V2
   COMPLEX (KIND=8) :: UT(2,2), U(2,2)

   N = Size(A,2)

   Allocate (V1(N), V2(N))
      
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
               A(I1,J) = A(I1,J) / XSIGMA_xL(NSIGL_U(I,NTAU)) ! i up
               A(I2,J) = A(I2,J) / XSIGMA_xR(NSIGL_U(I,NTAU)) ! i down
               
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
               A(I1,J) = V1(J) ! i left
               A(I2,J) = V2(J) ! i right

            ENDDO
         ENDDO
      ENDDO
   ENDIF

   IF (NFLAG.EQ.3.and.RHUB.GT.zero) THEN
      DO J = 1,N
         DO nx = 1,NLX
            DO ny = 1,NLY

               i  = invlist(nx,ny)
               i1 = INVNLIST(nx,ny,1,1) ! i left
               i2 = INVNLIST(nx,ny,1,2) ! i right

               ! U-term propagator in x eigen-basis
               A(I1,J) = A(I1,J) / XSIGMA_xL(NSIGL_U(I,NTAU)) ! i left
               A(I2,J) = A(I2,J) / XSIGMA_xR(NSIGL_U(I,NTAU)) ! i right

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
         I1 = L_Bonds(I,0  )
         I2 = L_Bonds(I,nf1)
         DO J = 1,N
            V1(J)   =  UT(1,1) * A(I1,J) + UT(1,2) * A(I2,J)
            V2(J)   =  UT(2,1) * A(I1,J) + UT(2,2) * A(I2,J) 
         ENDDO
         DO J = 1,N
            A(I1,J) = V1(J)
            A(I2,J) = V2(J)
         ENDDO
         
      ENDDO
   ENDIF
   
   IF (NFLAG.EQ.1 .AND. RJ .GT. ZERO) THEN
      DO I = 1,LFAM
         I1 = L_Bonds(I,0  )
         I2 = L_Bonds(I,nf1)
         ! Kinetic 
         DO J = 1,N
            A(I1,J) =  A(I1,J) / DCMPLX(XSIGP2(NSIGL_K(I,Nf1,NTAU)),0.D0)   
            A(I2,J) =  A(I2,J) / DCMPLX(XSIGM2(NSIGL_K(I,Nf1,NTAU)),0.D0)
         ENDDO
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
   
   Deallocate (V1, V2)
   
   Return
END SUBROUTINE MMUURM1

