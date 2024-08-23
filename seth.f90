Subroutine  SetH(HLP2)

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
        

end Subroutine SetH
