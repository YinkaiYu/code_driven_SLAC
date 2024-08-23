    SUBROUTINE PREQ(NOBS)
        Use Blockc
        Use Block_obs
        Implicit Real (KIND=8) (A-G,O-Z)
        Implicit Integer (H-N)
        COMPLEX (KIND =8) :: Znorm
        COMPLEX (KIND=8), Dimension(:,:,:), Allocatable :: Collect2
        REAL( kind=8) :: collect4
        Character(16) :: filek,filek1, filek2
        INTEGER :: NOBS
        real(kind=8) :: momx, momy
        Interface
            Subroutine Fourier_Trans(gr,filek)
                Complex (Kind=8), dimension(:,:,:) :: gr
                Character (16) :: filek
            end Subroutine Fourier_Trans
            Subroutine orderparameter(gr1, file)
                real(kind=8) :: gr1
                Character (16) :: file
            end subroutine orderparameter
            Subroutine orderparameter_complex(gr1, file)
                complex(kind=8) :: gr1
                Character (16) :: file
            end subroutine orderparameter_complex
            Subroutine structurefactor(gr,filek,mom_x,mom_y,no1,no2)
                Complex (Kind=8), dimension(:,:,:) :: gr
                Character (16) :: filek
                real(kind=8) :: mom_x,mom_y
                integer :: no1, no2
            end subroutine structurefactor
            Subroutine realspace(gr,filek)
                Complex (Kind=8), dimension(:,:,:) :: gr
                Character (16) :: filek
            end subroutine realspace
        end Interface
!#define DEC
        INCLUDE 'mpif.h'
        INTEGER STATUS(MPI_STATUS_SIZE)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
        ZNORM = CMPLX(1.D0,0.D0)/ DCMPLX( DBLE(NOBS), 0.D0 )
        spin = ZNORM * spin
        density = density / DBLE(NOBS)
        ferromagnetism = ferromagnetism / DBLE(NOBS)
        fermicor11 = fermicor11 / DBLE(NOBS)
        fermicor12 = fermicor12 / DBLE(NOBS)
        fermicor21 = fermicor21 / DBLE(NOBS)
        fermicor22 = fermicor22 / DBLE(NOBS)
        fermicor11_quarter = fermicor11_quarter / DBLE(NOBS)
        fermicor12_quarter = fermicor12_quarter / DBLE(NOBS)
        fermicor21_quarter = fermicor21_quarter / DBLE(NOBS)
        fermicor22_quarter = fermicor22_quarter / DBLE(NOBS)
        fermicor11_deltaq = fermicor11_deltaq / DBLE(NOBS)
        fermicor12_deltaq = fermicor12_deltaq / DBLE(NOBS)
        fermicor21_deltaq = fermicor21_deltaq / DBLE(NOBS)
        fermicor22_deltaq = fermicor22_deltaq / DBLE(NOBS)
        phasetot = phasetot/dble(Ncount)
        Allocate(Collect2(LQ,norb1,norb1))
        N = LQ*norb1*norb1 !
        Collect2 = CMPLX(0.D0,0.D0)
        CALL MPI_REDUCE(spin,Collect2,N,MPI_COMPLEX16,MPI_SUM,&
            & 0,MPI_COMM_WORLD,IERR)
        spin = Collect2/CMPLX(DBLE(ISIZE),0.D0)
        N=1
        Collect3 = 0.d0
        CALL MPI_REDUCE(phasetot,Collect3,N,MPI_REAL8,MPI_SUM,&
            & 0,MPI_COMM_WORLD,IERR)
        phasetot = Collect3/CMPLX(DBLE(ISIZE),0.D0) 
        collect4 = 0.0d0
        CALL MPI_REDUCE(density, Collect4, 1, MPI_REAL8, MPI_SUM, &
            &   0, MPI_COMM_WORLD, IERR)
        density = collect4/CMPLX(DBLE(ISIZE),0.D0)
        collect4 = 0.0d0
        CALL MPI_REDUCE(ferromagnetism, Collect4, 1, MPI_REAL8, MPI_SUM, &
            &   0, MPI_COMM_WORLD, IERR)
        ferromagnetism = collect4/CMPLX(DBLE(ISIZE),0.D0)
        collect4 = 0.0d0
        CALL MPI_REDUCE(fermicor11, Collect4, 1, MPI_REAL8, MPI_SUM, &
            &   0, MPI_COMM_WORLD, IERR)
        fermicor11 = collect4/CMPLX(DBLE(ISIZE),0.D0)
        collect4 = 0.0d0
        CALL MPI_REDUCE(fermicor12, Collect4, 1, MPI_REAL8, MPI_SUM, &
            &   0, MPI_COMM_WORLD, IERR)
        fermicor12 = collect4/CMPLX(DBLE(ISIZE),0.D0)
        collect4 = 0.0d0
        CALL MPI_REDUCE(fermicor21, Collect4, 1, MPI_REAL8, MPI_SUM, &
            &   0, MPI_COMM_WORLD, IERR)
        fermicor21 = collect4/CMPLX(DBLE(ISIZE),0.D0)
        collect4 = 0.0d0
        CALL MPI_REDUCE(fermicor22, Collect4, 1, MPI_REAL8, MPI_SUM, &
            &   0, MPI_COMM_WORLD, IERR)
        fermicor22 = collect4/CMPLX(DBLE(ISIZE),0.D0)
        collect4 = 0.0d0
        CALL MPI_REDUCE(fermicor11_quarter, Collect4, 1, MPI_REAL8, MPI_SUM, &
            &   0, MPI_COMM_WORLD, IERR)
        fermicor11_quarter = collect4/CMPLX(DBLE(ISIZE),0.D0)
        collect4 = 0.0d0
        CALL MPI_REDUCE(fermicor12_quarter, Collect4, 1, MPI_REAL8, MPI_SUM, &
            &   0, MPI_COMM_WORLD, IERR)
        fermicor12_quarter = collect4/CMPLX(DBLE(ISIZE),0.D0)
        collect4 = 0.0d0
        CALL MPI_REDUCE(fermicor21_quarter, Collect4, 1, MPI_REAL8, MPI_SUM, &
            &   0, MPI_COMM_WORLD, IERR)
        fermicor21_quarter = collect4/CMPLX(DBLE(ISIZE),0.D0)
        collect4 = 0.0d0
        CALL MPI_REDUCE(fermicor22_quarter, Collect4, 1, MPI_REAL8, MPI_SUM, &
            &   0, MPI_COMM_WORLD, IERR)
        fermicor22_quarter = collect4/CMPLX(DBLE(ISIZE),0.D0)
        collect4 = 0.0d0
        CALL MPI_REDUCE(fermicor11_deltaq, Collect4, 1, MPI_REAL8, MPI_SUM, &
            &   0, MPI_COMM_WORLD, IERR)
        fermicor11_deltaq = collect4/CMPLX(DBLE(ISIZE),0.D0)
        collect4 = 0.0d0
        CALL MPI_REDUCE(fermicor12_deltaq, Collect4, 1, MPI_REAL8, MPI_SUM, &
            &   0, MPI_COMM_WORLD, IERR)
        fermicor12_deltaq = collect4/CMPLX(DBLE(ISIZE),0.D0)
        collect4 = 0.0d0
        CALL MPI_REDUCE(fermicor21_deltaq, Collect4, 1, MPI_REAL8, MPI_SUM, &
            &   0, MPI_COMM_WORLD, IERR)
        fermicor21_deltaq = collect4/CMPLX(DBLE(ISIZE),0.D0)
        collect4 = 0.0d0
        CALL MPI_REDUCE(fermicor22_deltaq, Collect4, 1, MPI_REAL8, MPI_SUM, &
            &   0, MPI_COMM_WORLD, IERR)
        fermicor22_deltaq = collect4/CMPLX(DBLE(ISIZE),0.D0)
     
        IF (IRANK.EQ.0) THEN    

            filek = "density"
            Call orderparameter(density, filek )
            filek = "ferromagnetism"
            Call orderparameter(ferromagnetism, filek )
            filek = "fermicor11"
            Call orderparameter(fermicor11, filek )
            filek = "fermicor12"
            Call orderparameter(fermicor12, filek )
            filek = "fermicor21"
            Call orderparameter(fermicor21, filek )
            filek = "fermicor22"
            Call orderparameter(fermicor22, filek )
            filek = "fermicor11_quarter"
            Call orderparameter(fermicor11_quarter, filek )
            filek = "fermicor12_quarter"
            Call orderparameter(fermicor12_quarter, filek )
            filek = "fermicor21_quarter"
            Call orderparameter(fermicor21_quarter, filek )
            filek = "fermicor22_quarter"
            Call orderparameter(fermicor22_quarter, filek )
            filek = "fermicor11_deltaq"
            Call orderparameter(fermicor11_deltaq, filek )
            filek = "fermicor12_deltaq"
            Call orderparameter(fermicor12_deltaq, filek )
            filek = "fermicor21_deltaq"
            Call orderparameter(fermicor21_deltaq, filek )
            filek = "fermicor22_deltaq"
            Call orderparameter(fermicor22_deltaq, filek )
            filek = "phasetot"
            call orderparameter(phasetot,filek)
            filek = "FMzz"
            momx = 0.d0
            momy = 0.d0
            call structurefactor(spin,filek,momx,momy,1,1)
            filek = "FMzz_deltaq"
            momx = 1.0d0/dble(NLx)
            momy = 1.0d0/dble(NLy)
            call structurefactor(spin,filek,momx,momy,1,1)

        ENDIF
    END SUBROUTINE PREQ
    
    Subroutine Fourier_Trans(gr,filek)
        Use Blockc
        Use Block_obs
        Implicit Real (KIND=8) (A-G,O-Z)
        Implicit Integer (H-N)
        Complex (Kind=8), dimension(:) :: gr
        Integer :: lp
        Character (16) :: filek
        Real (Kind=8) :: xk_p(2), aimj_p(2)
        Complex (Kind=8), allocatable , dimension(:) :: gk
        allocate (gk(LQ))
        gk = cmplx(0.d0,0.d0)
        do imj = 1,LQ
            lx = list(imj,1)
            ly = list(imj,2)
            aimj_p(1) = dble(lx)* a1_p(1) + dble(ly)* a2_p(1)
            aimj_p(2) = dble(lx)* a1_p(2) + dble(ly)* a2_p(2)
            do nk = 1,LQ
                kx = list(nk,1)
                ky = list(nk,2)
                xk_p(1) = dble(kx - 1) * b1_p(1)/dble(NLX) + dble(ky - 1) * b2_p(1)/dble(NLX)
                xk_p(2) = dble(kx - 1) * b1_p(2)/dble(NLX) + dble(ky - 1) * b2_p(2)/dble(NLY)
                gk(nk) = gk(nk) + &
                    & exp( cmplx( 0.d0, xk_p(1)*aimj_p(1) + xk_p(2)*aimj_p(2) ) ) * gr(imj)                
            enddo
        enddo
        gk = gk/cmplx(LQ,0.d0)
        OPEN (UNIT=20,FILE=filek,STATUS='UNKNOWN', action="write", position="append")
        do nk = 1,LQ
            kx = list(nk,1)
            ky = list(nk,2)
            !  write(20,*)  kx,ky
            xk_p(1) = dble(kx - 1) * b1_p(1)/dble(NLX) + dble(ky - 1) * b2_p(1)/dble(NLY)
            xk_p(2) = dble(kx - 1) * b1_p(2)/dble(NLX) + dble(ky - 1) * b2_p(2)/dble(NLY)
!! to convert the FFA convention to correct convention
            write(20,*) xk_p(1), xk_p(2)
            write(20,*) gk(nk)
        enddo
        close(20)
        deallocate (gk)
    end Subroutine Fourier_Trans
    
    Subroutine orderparameter(gr1, file)
        Use Blockc
        Use Block_obs
        real(kind=8) :: gr1
        Character (16) :: file
        OPEN (UNIT=20,FILE=file,STATUS='UNKNOWN', action="write", position="append")
        write(20,*) gr1      
        close(20)
    end Subroutine orderparameter

    Subroutine orderparameter_complex(gr1, file)
        Use Blockc
        Use Block_obs
        complex(kind=8) :: gr1
        Character (16) :: file
        OPEN (UNIT=20,FILE=file,STATUS='UNKNOWN', action="write", position="append")
        write(20,*) gr1      
        close(20)
    end Subroutine orderparameter_complex
    
    Subroutine structurefactor(gr,filek,mom_x,mom_y,no1,no2)
         Use Blockc
         Use Block_obs
         Implicit Real (KIND=8) (A-G,O-Z)
         Implicit Integer (H-N)
         Complex (Kind=8), dimension(:,:,:) :: gr
         Integer :: lp
         Character (16) :: filek
         Real (Kind=8) :: mom_x,mom_y,xk_p(2), aimj_p(2)
         gk = cmplx(0.d0,0.d0)
         xk_p(1) = mom_x * b1_p(1) + mom_y * b2_p(1)
         xk_p(2) = mom_x * b1_p(2) + mom_y * b2_p(2)
         do imj = 1,LQ
             lx = list(imj,1)
             ly = list(imj,2)
            aimj_p(1) = dble(lx)* a1_p(1) + dble(ly)* a2_p(1)
            aimj_p(2) = dble(lx)* a1_p(2) + dble(ly)* a2_p(2)
            gk = gk + exp( cmplx( 0.d0, xk_p(1)*aimj_p(1)+xk_p(2)*aimj_p(2) ) ) * gr(imj,no1,no2)
         enddo
         gk = gk/cmplx(LQ,0.d0)
         OPEN (UNIT=20,FILE=filek,STATUS='UNKNOWN', action="write", position="append")
         write(20,*) real(gk)
         close(20)
    end subroutine structurefactor
    
    
        Subroutine realspace(gr,filek)
         Use Blockc
         Use Block_obs
         Implicit Real (KIND=8) (A-G,O-Z)
         Implicit Integer (H-N)
         Complex (Kind=8), dimension(:,:,:) :: gr
         Character (16) :: filek
    
         OPEN (UNIT=20,FILE=filek,STATUS='UNKNOWN', action="write", position="append")
         do imj = 1,LQ
             lx = list(imj,1)
             ly = list(imj,2)
             write(20,*) lx, ly
             do no1 = 1,norb
                 do no2 = 1, norb
                     write(20,*) gr(imj,no1,no2)
                 enddo
             enddo
         enddo
         close(20)
    end subroutine realspace