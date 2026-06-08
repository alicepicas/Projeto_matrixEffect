PROGRAM  tables_mie

!!# MG JIMENEZ, 2019. This version writes segment pair combinations for mdp file
!!# MG JIMENEZ, 2018. CG segment names with up to 5 characteres
!!# MG JIMENEZ, 2016, Based on CHErdes' version 2012

  IMPLICIT NONE

  INTEGER :: i, j, k,  nlines, nr

  DOUBLE PRECISION :: r, r2, rlr, rla, rlrb, rlab, lr, la
  DOUBLE PRECISION :: f1, f2, f3, factor, sij, sij3, s13, s23, epsij
  DOUBLE PRECISION :: sijla, sijlr, cgmx, agmx
  DOUBLE PRECISION :: cutoff, f, fprime, g, gprime, h, hprime
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: sig, eps, llr, lla, mas


  CHARACTER(LEN=18) :: file_name
  CHARACTER(LEN=5), DIMENSION(:), ALLOCATABLE :: seg


  OPEN( UNIT=1, FILE='data.inp', STATUS='old' )
  OPEN( UNIT=2, FILE='param_list.dat', STATUS='new' )
  OPEN( UNIT=3, FILE='gmx_self.dat', STATUS='new' )
  OPEN( UNIT=4, FILE='gmx_cross.dat', STATUS='new' )
  OPEN( UNIT=6, FILE='mdp_list.dat', STATUS='new' )

  ! ... reading lambda_r and lambda_a

  READ(1,*)
  READ(1,*) cutoff
  nr= INT(cutoff/0.002d0) + 1 

  READ(1,*)
  READ(1,*) nlines      
  ALLOCATE ( seg(nlines) )
  ALLOCATE ( sig(nlines) )
  ALLOCATE ( eps(nlines) )
  ALLOCATE ( llr(nlines) )
  ALLOCATE ( lla(nlines) )
  ALLOCATE ( mas(nlines) )

  READ(1,*)
  DO i=1, nlines
    READ(1,*) seg(i), sig(i), eps(i), llr(i), lla(i), mas(i)
  END DO

  WRITE( 2,* )"#  i    j    sig_ij*10    eps_ij    lr_ij    la_ij "
  WRITE( 3,*)"# Self interactions for GROMACS:            'C'            'A'"
  WRITE( 4,*)"# Cross interactions for GROMACS: 'C'     'A'"


  DO i=1, nlines
    DO j=i, nlines
      lr = 3 + ( ( llr(i)-3 )*( llr(j)-3) )**0.5d0
      la = 3 + ( ( lla(i)-3 )*( lla(j)-3) )**0.5d0
      f1 = lr/la
      f2 = 1.0d0/(lr-la)
      f3 = f1**(la*f2)
      factor = lr*f2*f3
      sij = ( sig(i) + sig(j) )/2.d0
      s13 = sig(i)**3
      s23 = sig(j)**3
      sij3 = sij**3
      sijla = sij**la
      sijlr = sij**lr
      epsij = SQRT( s13*s23*eps(i)*eps(j) )/sij3
      cgmx = factor*epsij*8.31451d0*sijla/1000.d0
      agmx = factor*epsij*8.31451d0*sijlr/1000.d0

      WRITE( 2,101 ) seg(i), seg(j), sij, epsij, lr, la
      WRITE( 6,107 ) seg(i), seg(j)

      IF ( i .EQ. j ) THEN
        WRITE( 3,102) seg(i), mas(i), 0.0d0, cgmx, agmx
      ELSE
        WRITE( 4,108) seg(i), seg(j), 1, cgmx, agmx
      END IF

      WRITE( file_name, *)"table_",trim(adjustl(seg(i))),"_",trim(adjustl(seg(j))),".xvg"
      OPEN( UNIT=5, FILE=trim(adjustl(file_name)), STATUS='new')
      WRITE(5,103)
      WRITE(5,104) lr, la
      WRITE(5,105) factor, lr, factor, la
    

      DO k=1, nr
        r = 0.002d0*(k-1)
        r2 = r*r
        rlr = r**lr
        rlrb = r**( lr+1 )
        rla = r**la
        rlab = r**( la+1 )
        f = 1/r
        fprime = 1/r2
        g = -1/rla
        gprime = -la/rlab
        h = 1/rlr
        hprime = lr/rlrb


        IF ( hprime .GT. 1e27) THEN
          WRITE(5,106) r, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        ELSE
          WRITE(5,106) r, f, fprime, g, gprime, h, hprime
        END IF
      END DO
     
      CLOSE(UNIT=5)

    END DO
  END DO

  CLOSE(UNIT=1)
  CLOSE(UNIT=2)
  CLOSE(UNIT=3)
  CLOSE(UNIT=4)
  CLOSE(UNIT=6)

101 FORMAT( 2x, 2(A5, 2x), F6.3, 4x, F7.2, 4x, 2(F6.2, 4x) )
102 FORMAT( 2x,A5,1x,"C",10x,F8.5,2x,F5.3,4x,"A",7x,2(E11.5, 5x) )
103 FORMAT('# MGJimenez, 2016. Based on CHerdes, 2012. Imperial College London' )
104 FORMAT('# Tabulated Mie Potential, lr=', F9.6, 'and la=', F9.6 )
105 FORMAT('# A = ', F9.6, 'x Epsilon x Sigma^', F9.6, ' and C = ', F9.6, 'x Epsilon x Sigma^', F9.6 )
106 FORMAT( 7(E17.10,3x) )
107 FORMAT( A5,1x,A5 )
108 FORMAT( 1x,A5,1x,A5,4x,I1,7x, 2(E11.5, 5x) )

  DEALLOCATE( seg )
  DEALLOCATE( sig )
  DEALLOCATE( eps )
  DEALLOCATE( llr )
  DEALLOCATE( lla )

END PROGRAM tables_mie
