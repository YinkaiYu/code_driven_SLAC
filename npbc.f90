      INTEGER FUNCTION NPBCX(NR)
        use blockc
        NPBCX = NR
        IF (NR.GT.NLX) NPBCX = NR - NLX
        IF (NR.LT.1) NPBCX = NR + NLX
        RETURN
    END FUNCTION NPBCX
    
    INTEGER FUNCTION NPBCY(NR)
        use blockc
        NPBCY = NR
        IF (NR.GT.NLY) NPBCY = NR - NLY
        IF (NR.LT.1) NPBCY = NR + NLY
        RETURN
    END FUNCTION NPBCY
    
     INTEGER FUNCTION NPBC_tau(NR)
        use blockc
       
       
        NPBC_tau = NR
        IF (NR.GT.LTROT) NPBC_tau = NR - LTROT
        IF (NR.LT.1) NPBC_tau = NR + LTROT
        
        RETURN
      END FUNCTION NPBC_tau

    INTEGER FUNCTION NPBC_SPIN(NR)
      use blockc
      NPBC_SPIN = NR
      IF (NR.GT.nspin) NPBC_SPIN = NR - nspin
      IF (NR.LT.1) NPBC_SPIN = NR + nspin
      RETURN
    END FUNCTION NPBC_SPIN   

    INTEGER FUNCTION NPBC_RX(NR)
      use blockc
      NPBC_RX = NR
      IF (NR.EQ.(NLX-1)) NPBC_RX = 1
      RETURN
    END FUNCTION

    INTEGER FUNCTION NPBC_RY(NR)
      use blockc
      NPBC_RY = NR
      IF (NR.EQ.(NLY-1)) NPBC_RY = 1
      RETURN
    END FUNCTION
