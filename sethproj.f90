Subroutine  SetHproj(HLP2)

   Use Blockc
   Use MyMats
   Implicit Real (KIND=8) (A-G,O-Z)
   Implicit Integer (H-N)
   

   Complex (Kind=8), Dimension(:,:) :: HLP2 
   Complex (Kind=8) :: Z1

   IseedHop = 3958195
   HLP2 = CMPLX(0.D0,0.D0)

   !! single dirac cone 
   do ix = 1, nlx
      do iy = 1, nly

         i  = invlist(ix,iy)
         ii = invnlist(ix,iy,1,1)      ! i_up

         do NRx = 1, nlx-1
            jx = npbcx(ix+NRx)
            jy = iy
            jj = invnlist(jx,jy,1,2)   ! j_down
            Z1 = CMPLX( 0.d0, (-1.d0)**dble(NRx) ) / ( (NLx/Pi) * sin(Pi*NRx/NLx) )
            if(Itwist==2) then
               random = ranf(IseedHop)
               Z1 = Z1 * CMPLX( 1.d0 + TwistX*(random - 0.5d0), 0.d0 )
            endif            
            if(Itwist==3) then
               random = ranf(IseedHop)
               Z1 = Z1 * exp( CMPLX( 0.d0, TwistX*(random - 0.5d0) ) )
            endif
            HLP2(ii,jj)  =  Z1         ! from i_up   to j_down
            HLP2(jj,ii)  =  conjg(Z1)  ! from j_down to i_up
         enddo

         do NRy = 1, nly-1
            jx = ix
            jy = npbcy(iy+NRy)
            jj = invnlist(jx,jy,1,2)   ! j_down
            Z1 = CMPLX( (-1.d0)**dble(NRy), 0.d0 ) / ( (NLy/Pi) * sin(Pi*NRy/NLy) ) 
            if(Itwist==2) then
               random = ranf(IseedHop)
               Z1 = Z1 * CMPLX( 1.d0 + TwistX*(random - 0.5d0), 0.d0 )
            endif            
            if(Itwist==3) then
               random = ranf(IseedHop)
               Z1 = Z1 * exp( CMPLX( 0.d0, TwistX*(random - 0.5d0) ) )
            endif
            HLP2(ii,jj)  =  Z1         ! from i_up   to j_down
            HLP2(jj,ii)  =  conjg(Z1)  ! from j_down to i_up
         enddo

      enddo
   enddo

   if(Itwist==1) then
      do ix = 1, nlx
         do iy = 1, nly

            random = ranf(IseedHop)

            i  = invlist(ix,iy)
            i_up = invnlist(ix,iy,1,1)
            i_do = invnlist(ix,iy,1,2)      

            ! HLP2(i_up,i_up) = HLP2(i_up,i_up) + TwistX*(random - 0.5d0)
            ! HLP2(i_do,i_do) = HLP2(i_do,i_do) + TwistX*(random - 0.5d0)
            HLP2(i_up,i_up) = HLP2(i_up,i_up) + CMPLX( TwistX*(random - 0.5d0), 0.d0 )
            HLP2(i_do,i_do) = HLP2(i_do,i_do) + CMPLX( TwistX*(random - 0.5d0), 0.d0 )

         enddo
      enddo
   endif
        

   !! FM initial
   if(Itwist==0) then
      HLP2 = CMPLX(0.D0,0.D0)
      do ix = 1, nlx
            do iy = 1, nly

               random = ranf(IseedHop)

               do ns = 1,nspin

                  i  = invlist(ix,iy)
                  ii = invnlist(ix,iy,1,ns)
                     
                  HLP2(ii,ii) = RHUB1*(-1.d0)**dble(ns) + TwistX*(random - 0.5d0)

               enddo
            enddo
      enddo
   endif

   !! RS initial: random up or down
   if(Itwist==-1) then
      HLP2 = CMPLX(0.D0,0.D0)
      do ix = 1, nlx
            do iy = 1, nly
               
               random = ranf(IseedHop)

               do ns = 1,nspin

                  i  = invlist(ix,iy)
                  ii = invnlist(ix,iy,1,ns)

                  HLP2(ii,ii) = RT1 * (-1.d0)**dble(ns) * (random - 0.5d0)

               enddo
            enddo
      enddo
   endif

   !! RS initial: up and down ramdom separately
   if(Itwist==-2) then
      HLP2 = CMPLX(0.D0,0.D0)
      do ix = 1, nlx
            do iy = 1, nly
               do ns = 1,nspin

                  i  = invlist(ix,iy)
                  ii = invnlist(ix,iy,1,ns)

                  random = ranf(IseedHop)
                  HLP2(ii,ii) = RT1 *  (random - 0.5d0)

               enddo
            enddo
      enddo
   endif


   ! write(6,*) 'matrix of HLP2'
   ! WRITE(6, '(A7)',advance='no') "Col/Row" 
   ! DO j = 1, Ndim
   !    WRITE(6,'(7X,I5,8X)', advance='no') j
   ! ENDDO
   ! WRITE(6,*) ! \n
   ! DO i = 1, Ndim
   !    WRITE(6, '(I5,2X)',advance='no') i 
   !    DO j = 1, Ndim
   !       WRITE(6,'( "(", F8.4, ",", F8.4, ") " )', advance='no') REAL(HLP2(i,j)), AIMAG(HLP2(i,j))
   !    ENDDO
   !    WRITE(6,*) ! \n
   ! ENDDO

   ! write(6,*) 'diag of HLP2'
   ! DO i = 1, Ndim
   !    WRITE(6,*) i, HLP2(i,i)
   ! ENDDO
          

end Subroutine SetHproj
