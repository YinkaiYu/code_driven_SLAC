SUBROUTINE UPGRADEU(NTAU,ISEED, UL,UR, ULRINV)

    Use Blockc
    Use Block_obs
    Implicit Real (KIND=8) (A-G,O-Z)  
    Implicit Integer (H-N)
    !Arguments.
    COMPLEX (Kind=8), Dimension(:,:) :: UL, UR, ULRINV
    INTEGER :: NTAU, ISEED
    
    !Space for tests
    COMPLEX (Kind=8), Dimension(:), Allocatable ::  V1, U1, VHLP1, UHLP1
    COMPLEX (Kind=8), Dimension(:), Allocatable ::  V2, U2, VHLP2, UHLP2
    COMPLEX (Kind=8)  RATIOUP, RATIOTOT, DENOM
    COMPLEX (Kind=8)  G44UP, G55UP, G45UP, G54UP, DEL44, DEL55, Z1, Z2, Z3, Z4
    
    ! Local
    Real (Kind=8) :: ACCM, Ratio_abs, Random, Weight
    Real (Kind=8), external :: Ranf
    Integer, external :: nranf
    Integer :: I4, NL, NL1, NL2

    ALLOCATE(V1(NE), U1(NE), VHLP1(NE), UHLP1(NE))
    ALLOCATE(V2(NE), U2(NE), VHLP2(NE), UHLP2(NE))

	ACCM  = 0.D0
	DO nx = 1, NLX
        DO ny = 1, NLY

            I  = INVLIST(NX,NY)
            I4 = INVNLIST(NX,NY,1,1)    ! spin left  in site I
            I5 = INVNLIST(NX,NY,1,2)    ! spin right in site I

            DEL44 = DELTA_xL( NSIGL_U(I,NTAU) )  
            DEL55 = DELTA_xR( NSIGL_U(I,NTAU) )

            G44UP = CMPLX(0.D0,0.D0)
            G45UP = CMPLX(0.D0,0.D0)
            G54UP = CMPLX(0.D0,0.D0)
            G55UP = CMPLX(0.D0,0.D0)
            ! ZGEMV
            DO NL = 1,NE
                VHLP1(NL) = DEL44*UR(I4,NL)
                VHLP2(NL) = DEL55*UR(I5,NL)
                UHLP1(NL) = UL(NL,I4)
                UHLP2(NL) = UL(NL,I5)
            ENDDO
            DO NL = 1,NE
                V1(NL) = CMPLX(0.D0,0.D0)
                V2(NL) = CMPLX(0.D0,0.D0)
                U1(NL) = CMPLX(0.D0,0.D0)
                U2(NL) = CMPLX(0.D0,0.D0)
            ENDDO
            DO NL = 1,NE
                DO NL1 = 1,NE
                V1(NL) = V1(NL) + VHLP1(NL1)*ULRINV(NL1,NL)
                V2(NL) = V2(NL) + VHLP2(NL1)*ULRINV(NL1,NL)
                ENDDO
            ENDDO
            DO NL = 1,NE
                G44UP = G44UP + V1(NL)*UHLP1(NL)
                G54UP = G54UP + V2(NL)*UHLP1(NL)
                G45UP = G45UP + V1(NL)*UHLP2(NL)
                G55UP = G55UP + V2(NL)*UHLP2(NL)
            ENDDO
            RATIOUP = (DCMPLX(1.D0,0.D0) + G44UP) * &
                        &    (DCMPLX(1.D0,0.D0) + G55UP) - &
                        &     G45UP*G54UP

            RATIOTOT  = RATIOUP ! C/C=1
            RATIO_ABS = abs(RATIOTOT)

            Random = RANF(ISEED)
            IF ( RATIO_ABS.GT.Random ) THEN

                ACCM  = ACCM + 1.D0

                ! Upgrade Phase.
                WEIGHT = SQRT(DBLE(RATIOTOT*DCONJG(RATIOTOT)))
                phase = phase * ratiotot/dcmplx(weight,0.d0)
                phasetot = phasetot + real(phase)
                ncount = ncount + 1

                DO NL  = 1,NE
                    DO NL1 = 1,NE
                        U1(NL) = U1(NL) + ULRINV(NL,NL1)*UHLP1(NL1)
                        U2(NL) = U2(NL) + ULRINV(NL,NL1)*UHLP2(NL1)
                    ENDDO
                ENDDO

                Z1 =  CMPLX(1.D0,0.D0)/(CMPLX(1.D0,0.D0) + G55UP)
                Z2 =  G54UP*Z1
                Z3 =  G45UP*Z1
                Z4 =  DCMPLX(1.D0,0.D0) + G44UP - G45UP*G54UP*Z1
                Z4 =  DCMPLX(1.D0,0.D0)/Z4
                DO NL = 1,NE
                    UHLP1(NL) = U2(NL)
                    VHLP1(NL) = V2(NL)*Z1
                    UHLP2(NL) = Z4*( U1(NL) - U2(NL)*Z2 )
                    VHLP2(NL) =      V1(NL) - V2(NL)*Z3
                ENDDO
                DO NL1 = 1,NE
                    DO NL2 = 1,NE
                        ULRINV(NL1,NL2) = ULRINV(NL1,NL2) &
                        &        -  UHLP1(NL1)*VHLP1(NL2)&
                        &        -  UHLP2(NL1)*VHLP2(NL2)
                    ENDDO
                ENDDO

                ! Upgrade  UR
                DO NL = 1,NE
                    UR(I4,NL) = UR(I4,NL) + &    
                    &     DEL44*UR(I4,NL)
                    UR(I5,NL) = UR(I5,NL) +  &      
                    &     DEL55*UR(I5,NL)   
                ENDDO

                ! Flip:              
                NSIGL_U(I,NTAU) =  NFLIPL( NSIGL_U(I,NTAU) )

            ENDIF

        ENDDO
	ENDDO
	OBS(27) = OBS(27) + DCMPLX(ACCM/DBLE(Ndim),0.D0)
	OBS(28) = OBS(28) + DCMPLX(1.D0,0.D0)

    DEALLOCATE(V1, U1, VHLP1, UHLP1)
    DEALLOCATE(V2, U2, VHLP2, UHLP2)
	 
END SUBROUTINE UPGRADEU
