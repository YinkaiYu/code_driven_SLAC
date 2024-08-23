SUBROUTINE OBSER
    Use Blockc
    Use Block_obs
    Use MyMats
    Implicit Real (KIND=8) (A-G,O-Z)
    Implicit Integer (H-N)
    
    complex(kind=8) :: Z1
    Real (Kind=8) :: mom_x,mom_y,xk_p(2), aimj_p(2)

    
    do ix = 1, nlx
        do iy = 1, nly

            iup = invnlist(ix,iy,1,1)
            ido = invnlist(ix,iy,1,2)
            density = density + real(phase) * real( GRUPC(iup,iup) + GRUPC(ido,ido) )/dble(LQ)
            ferromagnetism = ferromagnetism + real(phase) * real( GRUPC(iup,iup) - GRUPC(ido,ido) )/dble(LQ)

            ! fermion correlation function (sciadv.aau1463)
            jx = npbcx( ix + (nlx-1)/2 )
            jy = npbcx( iy + (nly-1)/2 )
            jup = invnlist(jx,jy,1,1)
            jdo = invnlist(jx,jy,1,2)
            fermicor11 = fermicor11 + real(phase) * real( GRUPC(iup,jup) + GRUPC(jup,iup) )/dble(LQ)
            fermicor12 = fermicor12 + real(phase) * real( GRUPC(iup,jdo) + GRUPC(jdo,iup) )/dble(LQ)
            fermicor21 = fermicor21 + real(phase) * real( GRUPC(ido,jup) + GRUPC(jup,ido) )/dble(LQ)
            fermicor22 = fermicor22 + real(phase) * real( GRUPC(ido,jdo) + GRUPC(jdo,ido) )/dble(LQ)

            
            jx = npbcx( ix + (nlx-1)/4 )
            jy = npbcx( iy + (nly-1)/4 )
            jup = invnlist(jx,jy,1,1)
            jdo = invnlist(jx,jy,1,2)
            fermicor11_quarter = fermicor11_quarter + real(phase) * real( GRUPC(iup,jup) + GRUPC(jup,iup) )/dble(LQ)
            fermicor12_quarter = fermicor12_quarter + real(phase) * real( GRUPC(iup,jdo) + GRUPC(jdo,iup) )/dble(LQ)
            fermicor21_quarter = fermicor21_quarter + real(phase) * real( GRUPC(ido,jup) + GRUPC(jup,ido) )/dble(LQ)
            fermicor22_quarter = fermicor22_quarter + real(phase) * real( GRUPC(ido,jdo) + GRUPC(jdo,ido) )/dble(LQ)

        enddo
    enddo

    mom_x = 1.0d0*PI/dble(NLx)
    mom_y = 1.0d0*PI/dble(NLy)
    xk_p(1) = mom_x * b1_p(1) + mom_y * b2_p(1)
    xk_p(2) = mom_x * b1_p(2) + mom_y * b2_p(2)

    do ix = 1, nlx
        do  iy = 1, nly
            do jx = 1, nlx
                do jy = 1, nly  
                    i     = invlist(ix,iy)                
                    j     = invlist(jx,jy)
                    iup   = invnlist(ix,iy,1,1)
                    ido   = invnlist(ix,iy,1,2)
                    jup   = invnlist(jx,jy,1,1)
                    jdo   = invnlist(jx,jy,1,2)
                    imj_x = ix-jx ! no npbcx here
                    imj_y = iy-jy ! no npbcy here
                    imj   = invlist(npbcx(imj_x),npbcx(imj_y))
                    spin(imj,1,1) = spin(imj,1,1) + phase * ( &
                        &   GRUPC(iup,jup) * GRUP(iup,jup) + GRUPC(iup,iup) * GRUPC(jup,jup) &
                        & + GRUPC(ido,jdo) * GRUP(ido,jdo) + GRUPC(ido,ido) * GRUPC(jdo,jdo) &
                        & - GRUPC(iup,jdo) * GRUP(iup,jdo) - GRUPC(iup,iup) * GRUPC(jdo,jdo) &
                        & - GRUPC(ido,jup) * GRUP(ido,jup) - GRUPC(ido,ido) * GRUPC(jup,jup) &
                    & )/dble(LQ)
                    ! fermion correlation function (PRL 128, 225701 (2022))
                    aimj_p(1) = dble(imj_x)* a1_p(1) + dble(imj_y)* a2_p(1)
                    aimj_p(2) = dble(imj_x)* a1_p(2) + dble(imj_y)* a2_p(2)
                    fermicor11_deltaq = fermicor11_deltaq + real( exp( cmplx( 0.d0, xk_p(1)*aimj_p(1)+xk_p(2)*aimj_p(2) ) ) * GRUPC(iup,jup) * phase ) / dble(LQ)
                    fermicor12_deltaq = fermicor12_deltaq + real( exp( cmplx( 0.d0, xk_p(1)*aimj_p(1)+xk_p(2)*aimj_p(2) ) ) * GRUPC(iup,jdo) * phase ) / dble(LQ)
                    fermicor21_deltaq = fermicor21_deltaq + real( exp( cmplx( 0.d0, xk_p(1)*aimj_p(1)+xk_p(2)*aimj_p(2) ) ) * GRUPC(ido,jup) * phase ) / dble(LQ)
                    fermicor22_deltaq = fermicor22_deltaq + real( exp( cmplx( 0.d0, xk_p(1)*aimj_p(1)+xk_p(2)*aimj_p(2) ) ) * GRUPC(ido,jdo) * phase ) / dble(LQ)
                enddo
            enddo    
        enddo
    enddo

    ! write(6,*) 'matrix of GRUPC'
    ! WRITE(6, '(A7)',advance='no') "Col/Row" 
    ! DO j = 1, Ndim
    !     WRITE(6,'(7X,I5,8X)', advance='no') j
    ! ENDDO
    ! WRITE(6,*) ! \n
    ! DO i = 1, Ndim
    !     WRITE(6, '(I5,2X)',advance='no') i 
    !     DO j = 1, Ndim
    !         WRITE(6,'( "(", F8.4, ",", F8.4, ") " )', advance='no') REAL(GRUPC(i,j)), AIMAG(GRUPC(i,j))
    !     ENDDO
    !     WRITE(6,*) ! \n
    ! ENDDO

END SUBROUTINE OBSER