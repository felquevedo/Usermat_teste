*deck,usermat      USERDISTRIB  parallel                                gal
      subroutine usermat(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7, var8)
c*************************************************************************
c     *** primary function ***
c
c           user defined material constitutive model
c
c      Attention:
c           User must define material constitutive law properly
c           according to the stress state such as 3D, plane strain
c           and axisymmetry, plane stress and 3D/1D beam.
c
c           A 3D material constitutive model can be used for
c           plane strain and axisymmetry cases.
c
c           When using shell elements, a plane stress algorithm
c           must be used.
c
c                                             gal July, 1999
c
c       The following demonstrates a USERMAT subroutine for 
c       a plasticity model, which is the same as TB, BISO,
c       for different stress states. 
c       See "ANSYS user material subroutine USERMAT" for detailed
c       description of how to write a USERMAT routine.
c
c       This routine calls four routines,
c       usermat3d.F, usermatps.F usermatbm.F and usermat1d.F, w.r.t.
c       the corresponding stress states.
c       Each routine can be also a usermat routine for the specific 
c       element.
c
c*************************************************************************
c
c     input arguments
c     ===============
c      matId     (int,sc,i)               material #
c      elemId    (int,sc,i)               element #
c      kDomIntPt (int,sc,i)               "k"th domain integration point
c      kLayer    (int,sc,i)               "k"th layer
c      kSectPt   (int,sc,i)               "k"th Section point
c      ldstep    (int,sc,i)               load step number
c      isubst    (int,sc,i)               substep number
c      nDirect   (int,sc,in)              # of direct components
c      nShear    (int,sc,in)              # of shear components
c      ncomp     (int,sc,in)              nDirect + nShear
c      nstatev   (int,sc,i)               Number of state variables
c      nProp     (int,sc,i)               Number of material constants
c
c      Temp      (dp,sc,in)               temperature at beginning of
c                                         time increment
c      dTemp     (dp,sc,in)               temperature increment 
c      Time      (dp,sc,in)               time at beginning of increment (t)
c      dTime     (dp,sc,in)               current time increment (dt)
c
c      Strain   (dp,ar(ncomp),i)          Strain at beginning of time increment
c      dStrain  (dp,ar(ncomp),i)          Strain increment
c      prop     (dp,ar(nprop),i)          Material constants defined by TB,USER
c      coords   (dp,ar(3),i)              current coordinates
c      defGrad_t(dp,ar(3,3),i)            Deformation gradient at time t
c      defGrad  (dp,ar(3,3),i)            Deformation gradient at time t+dt
c
c     input output arguments              
c     ======================             
c      stress   (dp,ar(ncomp),io)         stress
c      ustatev   (dp,ar(nstatev),io)      user state variables
c      sedEl    (dp,sc,io)                elastic work
c      sedPl    (dp,sc,io)                plastic work
c      epseq    (dp,sc,io)                equivalent plastic strain
c      epsPl   (dp,ar(ncomp),io)          plastic strain
c      var?     (dp,sc,io)                not used, they are reserved arguments 
c                                         for further development
c
c     output arguments
c     ================
c      keycut   (int,sc,o)                loading bisect/cut control
c                                         0 - no bisect/cut
c                                         1 - bisect/cut 
c                                         (factor will be determined by solution control)
c      dsdePl   (dp,ar(ncomp,ncomp),o)    material jacobian matrix
c      tsstif   (dp,ar(2),o)              transverse shear stiffness
c                                         tsstif(1) - Gxz
c                                         tsstif(2) - Gyz
c                                         tsstif(1) is also used to calculate hourglass
c                                         stiffness, this value must be defined when low
c                                         order element, such as 181, 182, 185 with uniform 
c                                         integration is used.
c      epsZZ    (dp,sc,o)                 strain epsZZ for plane stress,
c                                         define it when accounting for thickness change
c                                         in shell and plane stress states
c
c*************************************************************************
c
c      ncomp   6   for 3D  (nshear=3)
c      ncomp   4   for plane strain or axisymmetric (nShear = 1)
c      ncomp   3   for plane stress (nShear = 1)
c      ncomp   3   for 3d beam      (nShear = 2)
c      ncomp   1   for 1D (nShear = 0)
c
c      stresses and strains, plastic strain vectors
c          11, 22, 33, 12, 23, 13    for 3D
c          11, 22, 33, 12            for plane strain or axisymmetry
c          11, 22, 12                for plane stress
c          11, 13, 12                for 3d beam
c          11                        for 1D
c
c      material jacobian matrix
c        3D
c           dsdePl    |  1111   1122   1133   1112   1123   1113 |
c           dsdePl    |  2211   2222   2233   2212   2223   2213 |
c           dsdePl    |  3311   3322   3333   3312   3323   3313 |
c           dsdePl    |  1211   1222   1233   1212   1223   1213 |
c           dsdePl    |  2311   2322   2333   2312   2323   2313 |
c           dsdePl    |  1311   1322   1333   1312   1323   1313 |
c        plane strain or axisymmetric (11, 22, 33, 12)
c           dsdePl    |  1111   1122   1133   1112 |
c           dsdePl    |  2211   2222   2233   2212 |
c           dsdePl    |  3311   3322   3333   3312 |
c           dsdePl    |  1211   1222   1233   1212 |
c        plane stress (11, 22, 12)
c           dsdePl    |  1111   1122   1112 |
c           dsdePl    |  2211   2222   2212 |
c           dsdePl    |  1211   1222   1212 |
c        3d beam (11, 13, 12)
c           dsdePl    |  1111   1113   1112 |
c           dsdePl    |  1311   1313   1312 |
c           dsdePl    |  1211   1213   1212 |
c        1d
c           dsdePl    |  1111 |
c
c*************************************************************************
#include "impcom.inc"
c
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),       
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
c
      EXTERNAL         usermat3d, usermatps, usermatbm, usermat1d

      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7, var8
c
c*************************************************************************
c
      IF(ncomp .GE. 4) THEN
c ***    3d, plane strain and axisymmetric example
         call usermat3d (
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords,
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7, var8)

      ELSE IF(nDirect.eq. 2 .and. ncomp .EQ. 3) THEN
c ***    plane stress example
         call usermatps (
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords,
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7, var8)

      ELSE IF(ncomp .EQ. 3) THEN
c ***    3d beam example
         call usermatbm (
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords,
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7, var8)

      ELSE IF(ncomp .EQ. 1) THEN
c ***    1d beam example
         call usermat1d (
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords,
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7, var8)

      END IF
      return
      end
*deck,usermat1d    USERDISTRIB  parallel                                gal
      subroutine usermat1d(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7, var8)
c*************************************************************************
c     *** primary function ***
c
c           user defined material constitutive model
c
c      Attention:
c           User must define material constitutive law properly
c           according to the stress state such as 3D, plane strain
c           and axisymmetry, plane stress and beam.
c
c           a 3D material constitutive model can use for
c           plane strain and axisymmetry cases.
c
c           When using shell elements, a plane stress algorithm
c           must be use.
c
c                                             gal July, 1999
c
c       The following demonstrates a USERMAT subroutine for
c       a plasticity model of 1D truss element (LINK180). 
c       The plasticity model is the same as TB, BISO.
c
c       See "ANSYS user material subroutine USERMAT" for detailed
c       description of how to write a USERMAT routine.
c
c*************************************************************************
c
c     input arguments
c     ===============
c      matId     (int,sc,i)               material #
c      elemId    (int,sc,i)               element #
c      kDomIntPt (int,sc,i)               "k"th domain integration point
c      kLayer    (int,sc,i)               "k"th layer
c      kSectPt   (int,sc,i)               "k"th Section point
c      ldstep    (int,sc,i)               load step number
c      isubst    (int,sc,i)               substep number
c      nDirect   (int,sc,in)              # of direct components
c      nShear    (int,sc,in)              # of shear components
c      ncomp     (int,sc,in)              nDirect + nShear
c      nstatev   (int,sc,l)               Number of state variables
c      nProp     (int,sc,l)               Number of material ocnstants
c
c      Temp      (dp,sc,in)               temperature at beginning of
c                                         time increment
c      dTemp     (dp,sc,in)               temperature increment 
c      Time      (dp,sc,in)               time at beginning of increment (t)
c      dTime     (dp,sc,in)               current time increment (dt)
c
c      Strain   (dp,ar(ncomp),i)          Strain at beginning of time increment
c      dStrain  (dp,ar(ncomp),i)          Strain increment
c      prop     (dp,ar(nprop),i)          Material constants defined by TB,USER
c      coords   (dp,ar(3),i)              current coordinates
c      defGrad_t(dp,ar(3,3),i)            Deformation gradient at time t
c      defGrad  (dp,ar(3,3),i)            Deformation gradient at time t+dt
c
c     input output arguments              
c     ======================             
c      stress   (dp,ar(nTesn),io)         stress
c      ustatev   (dp,ar(nstatev),io)      user state variables
c            ustatev(1)                     - equivalent plastic strain
c            ustatev(2) - ustatev(1+ncomp)   - plastic strain vector
c            ustatev(nStatev)               - von-Mises stress
c      sedEl    (dp,sc,io)                elastic work
c      sedPl    (dp,sc,io)                plastic work
c      epseq    (dp,sc,io)                equivalent plastic strain
c      tsstif   (dp,ar(2),io)             transverse shear stiffness
c                                         tsstif(1) - Gxz
c                                         tsstif(2) - Gyz
c                                         tsstif(1) is also used to calculate hourglass
c                                         stiffness, this value must be defined when low
c                                         order element, such as 181, 182, 185 with uniform 
c                                         integration is used.
c      var?     (dp,sc,io)                not used, they are reserved arguments 
c                                         for further development
c
c     output arguments
c     ================
c      keycut   (int,sc,io)               loading bisect/cut control
c                                         0 - no bisect/cut
c                                         1 - bisect/cut 
c                                         (factor will be determined by ANSYS solution control)
c      dsdePl   (dp,ar(ncomp,ncomp),io)   material jacobian matrix
c      epsZZ    (dp,sc,o)                 strain epsZZ for plane stress,
c                                         define it when accounting for thickness change 
c                                         in shell and plane stress states
c
c*************************************************************************
c
c      ncomp   6   for 3D  (nshear=3)
c      ncomp   4   for plane strain or axisymmetric (nShear = 1)
c      ncomp   3   for plane stress (nShear = 1)
c      ncomp   3   for 3d beam      (nShear = 2)
c      ncomp   1   for 1D (nShear = 0)
c
c      stresss and strains, plastic strain vectors
c          11, 22, 33, 12, 23, 13    for 3D
c          11, 22, 33, 12            for plane strain or axisymmetry
c          11, 22, 12                for plane stress
c          11, 13, 12                for 3d beam
c          11                        for 1D
c
c      material jacobian matrix
c        3D
c           dsdePl    |  1111   1122   1133   1112   1123   1113 |
c           dsdePl    |  2211   2222   2233   2212   2223   2213 |
c           dsdePl    |  3311   3322   3333   3312   3323   3313 |
c           dsdePl    |  1211   1222   1233   1212   1223   1213 |
c           dsdePl    |  2311   2322   2333   2312   2323   2313 |
c           dsdePl    |  1311   1322   1333   1312   1323   1313 |
c        plane strain or axisymmetric (11, 22, 33, 12)
c           dsdePl    |  1111   1122   1133   1112 |
c           dsdePl    |  2211   2222   2233   2212 |
c           dsdePl    |  3311   3322   3333   3312 |
c           dsdePl    |  1211   1222   1233   1212 |
c        plane stress (11, 22, 12)
c           dsdePl    |  1111   1122   1112 |
c           dsdePl    |  2211   2222   2212 |
c           dsdePl    |  1211   1222   1212 |
c        3d beam (11, 13, 12)
c           dsdePl    |  1111   1113   1112 |
c           dsdePl    |  1311   1313   1312 |
c           dsdePl    |  1211   1213   1212 |
c        1d
c           dsdePl    |  1111 |
c
c*************************************************************************
#include "impcom.inc"
c
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
c
c***************** User defined part *************************************
c
c --- parameters
c
      INTEGER          mcomp
      DOUBLE PRECISION ZERO, HALF, ONE, TWO, SMALL
      PARAMETER       (ZERO       = 0.d0,
     &                 HALF       = 0.5d0,
     &                 ONE        = 1.d0,
     &                 TWO        = 2.d0,
     &                 SMALL      = 1.d-08,
     &                 mcomp      = 1
     &                 )
c
c --- local variables
c
c      sigElp   (dp,ar(6  ),l)            trial stress
c      dsdeEl   (dp,ar(6,6),l)            elastic moduli
c      sigDev   (dp,ar(6  ),l)            deviatoric stress tensor
c      dfds     (dp,ar(6  ),l)            derivative of the yield function 
c      JM       (dp,ar(6,6),l)            2D matrix for a 4 order tensor
c      pEl      (dp,sc     ,l)            hydrostatic pressure stress
c      qEl      (dp,sc     ,l)            von-mises stress
c      pleq_t   (dp,sc     ,l)            equivalent plastic strain at beginnig of time increment
c      pleq     (dp,sc     ,l)            equivalent plastic strain at end of time increment
c      dpleq    (dp,sc     ,l)            incremental equivalent plastic strain
c      sigy_t   (dp,sc     ,l)            yield stress at beginnig of time increments
c      sigy     (dp,sc     ,l)            yield stress at end of time increment
c      young    (dp,sc     ,l)            Young's modulus
c      posn     (dp,sc     ,l)            Poiss's ratio
c      sigy0    (dp,sc     ,l)            initial yield stress
c      dsigdep  (dp,sc     ,l)            plastic slop
c      twoG     (dp,sc     ,l)            two time of shear moduli
c
c
      DOUBLE PRECISION sigElp(mcomp), dsdeEl(mcomp,mcomp)

      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7, var8

      DOUBLE PRECISION qEl,   pleq_t,  sigy_t , sigy,
     &                 dpleq, pleq,    signTens,
     &                 young, posn,    sigy0,   dsigdep, 
     &                 twoG,  fratio
c*************************************************************************
c
      keycut   = 0
      dsigdep  = ZERO 
      pleq_t   = ustatev(1)
      pleq     = pleq_t
c *** get Young's modulus and Poisson's ratio, initial yield stress and others
      young    = prop(1)
      posn     = prop(2)
      sigy0    = prop(3)
c *** calculate plastic slope
      dsigdep  = young*prop(4)/(young-prop(4))
      twoG     = young / (ONE+posn)
c *** define tsstif(1) since it is used for calculation of hourglass stiffness
      tsstif(1) = HALF * twoG
c
c *** calculate elastic stiffness matrix 
c
      dsdeEl(1,1)= young
c
c *** calculate the trial stress and 
c     copy elastic moduli dsdeEl to material Jacobian matrix
      sigElp(1)   = stress(1)
      dsdePl(1,1) = dsdeEl(1,1)
      sigElp(1)   = sigElp(1) + dsdeEl(1,1) * dStrain(1)
c *** sign of predicted stress
      signTens = sign (ONE, sigElp(1))
c *** compute von-mises equivalent stress
      qEl = abs(sigElp(1))
c *** compute current yield stress
      sigy    = sigy0 + dsigdep * pleq
c
      fratio = qEl / sigy - ONE
c *** check for yielding
      IF (sigy .LE. ZERO.or.fratio .LE. -SMALL) GO TO 500
c
      sigy_t   = sigy
c *** initial guess of incremental equivalent plastic strain   
      dpleq    = (qEl - sigy) / young
      pleq     = pleq_t + dpleq
      sigy     = sigy0 + dsigdep * pleq
c
c ***  update plastic strains, stresses
      epsPl(1) = epsPl(1) + dpleq * signTens
      stress(1) =  signTens * sigy
c
c ***  update plastic strains
      epseq  = pleq
c *** Update state variables
      ustatev(1) = pleq
      ustatev(2) = epsPl(1)
c *** Update plastic work
      sedPl = sedPl + HALF * (sigy_t + sigy) * dpleq
c
c *** Material Jcobian matrix
c
      dsdePl(1,1) = dsdeEl(1,1) * dsigdep /(dsdeEl(1,1) + dsigdep)
c *** Allow a small number for Jcobian matrix if it is an ideal plasticity
      if(dsdePl(1,1).LE.ZERO) dsdePl(1,1) = SMALL*dsdeEl(1,1)
c
      goto 600
  500 continue

c *** Update stress in case of elastic/unloading
      stress(1) = sigElp(1)

  600 continue
c *** elastic strain energy
      sedEl = HALF * stress(1) * (Strain(1)+dStrain(1)-epsPl(1))
c *** update state variables
      ustatev(nStatev) = sigy
c
      return
      end
*deck,usermat3d    USERDISTRIB  parallel                                gal
      subroutine usermat3d(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ, cutFactor, 
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7)
c*************************************************************************
c     *** primary function ***
c
c           user defined material constitutive model
c
c      Attention:
c           User must define material constitutive law properly
c           according to the stress state such as 3D, plane strain
c           and axisymmetry, plane stress and beam.
c
c           a 3D material constitutive model can use for
c           plane strain and axisymmetry cases.
c
c           When using shell elements, a plane stress algorithm
c           must be use.
c
c                                             gal July, 1999
c
c       The following demonstrates a USERMAT subroutine for
c       a plasticity model of 3D solid elements or plane elements
c       in plane strain or axisymmetric stress state. The plasticity
c       model is the same as TB, BISO.
c       See "ANSYS user material subroutine USERMAT" for detailed
c       description of how to write a USERMAT routine.
c
c*************************************************************************
c
c     input arguments
c     ===============
c      matId     (int,sc,i)               material #
c      elemId    (int,sc,i)               element #
c      kDomIntPt (int,sc,i)               "k"th domain integration point
c      kLayer    (int,sc,i)               "k"th layer
c      kSectPt   (int,sc,i)               "k"th Section point
c      ldstep    (int,sc,i)               load step number
c      isubst    (int,sc,i)               substep number
c      nDirect   (int,sc,in)              # of direct components
c      nShear    (int,sc,in)              # of shear components
c      ncomp     (int,sc,in)              nDirect + nShear
c      nstatev   (int,sc,l)               Number of state variables
c      nProp     (int,sc,l)               Number of material ocnstants
c
c      Temp      (dp,sc,in)               temperature at beginning of
c                                         time increment
c      dTemp     (dp,sc,in)               temperature increment 
c      Time      (dp,sc,in)               time at beginning of increment (t)
c      dTime     (dp,sc,in)               current time increment (dt)
c
c      Strain   (dp,ar(ncomp),i)          Strain at beginning of time increment
c      dStrain  (dp,ar(ncomp),i)          Strain increment
c      prop     (dp,ar(nprop),i)          Material constants defined by TB,USER
c      coords   (dp,ar(3),i)              current coordinates
c      defGrad_t(dp,ar(3,3),i)            Deformation gradient at time t
c      defGrad  (dp,ar(3,3),i)            Deformation gradient at time t+dt
c
c     input output arguments              
c     ======================             
c      stress   (dp,ar(nTesn),io)         stress
c      ustatev   (dp,ar(nstatev),io)      user state variable
c            ustatev(1)                     - equivalent plastic strain
c            ustatev(2) - statev(1+ncomp)   - plastic strain vector
c            ustatev(nStatev)               - von-Mises stress
c      sedEl    (dp,sc,io)                elastic work
c      sedPl    (dp,sc,io)                plastic work
c      epseq    (dp,sc,io)                equivalent plastic strain
c      tsstif   (dp,ar(2),io)             transverse shear stiffness
c                                         tsstif(1) - Gxz
c                                         tsstif(2) - Gyz
c                                         tsstif(1) is also used to calculate
c                                          hourglass
c                                         stiffness, this value must be defined when 
c                                          low
c                                         order element, such as 181, 182, 185 with 
c                                          uniform 
c                                         integration is used.
c      var?     (dp,sc,io)                not used, they are reserved arguments 
c                                         for further development
c
c     output arguments
c     ================
c      keycut   (int,sc,io)               loading bisect/cut control
c                                         0 - no bisect/cut
c                                         1 - bisect/cut 
c                                         (factor will be determined by ANSYS solution 
c                                          control)
c      dsdePl   (dp,ar(ncomp,ncomp),io)   material jacobian matrix
c      epsZZ    (dp,sc,o)                 strain epsZZ for plane stress,
c                                         define it when accounting for thickness 
c                                          change 
c                                         in shell and plane stress states
c      cutFactor(dp,sc,o)                 time step size cut-back factor 
c                                         define it if a smaller step size is wished
c                                         recommended value is 0~1
c
c*************************************************************************
c
c      ncomp   6   for 3D  (nshear=3)
c      ncomp   4   for plane strain or axisymmetric (nShear = 1)
c      ncomp   3   for plane stress (nShear = 1)
c      ncomp   3   for 3d beam      (nShear = 2)
c      ncomp   1   for 1D (nShear = 0)
c
c      stresss and strains, plastic strain vectors
c          11, 22, 33, 12, 23, 13    for 3D
c          11, 22, 33, 12            for plane strain or axisymmetry
c
c      material jacobian matrix
c        3D
c           dsdePl    |  1111   1122   1133   1112   1123   1113 |
c           dsdePl    |  2211   2222   2233   2212   2223   2213 |
c           dsdePl    |  3311   3322   3333   3312   3323   3313 |
c           dsdePl    |  1211   1222   1233   1212   1223   1213 |
c           dsdePl    |  2311   2322   2333   2312   2323   2313 |
c           dsdePl    |  1311   1322   1333   1312   1323   1313 |
c        plane strain or axisymmetric (11, 22, 33, 12)
c           dsdePl    |  1111   1122   1133   1112 |
c           dsdePl    |  2211   2222   2233   2212 |
c           dsdePl    |  3311   3322   3333   3312 |
c           dsdePl    |  1211   1222   1233   1212 |
c
c*************************************************************************
#include "impcom.inc"
      implicit none
c
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ,   cutFactor
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
c
c***************** User defined part *************************************
c
c --- parameters
c
      INTEGER          mcomp
      DOUBLE PRECISION HALF, THIRD, ONE, TWO, SMALL, ONEHALF,
     &                 ZERO, TWOTHIRD, ONEDM02, ONEDM05, sqTiny
      PARAMETER       (ZERO       = 0.d0,
     &                 HALF       = 0.5d0,
     &                 THIRD      = 1.d0/3.d0,
     &                 ONE        = 1.d0,
     &                 TWO        = 2.d0,
     &                 SMALL      = 1.d-08,
     &                 sqTiny     = 1.d-20,
     &                 ONEDM02    = 1.d-02,
     &                 ONEDM05    = 1.d-05,
     &                 ONEHALF    = 1.5d0,
     &                 TWOTHIRD   = 2.0d0/3.0d0,
     &                 mcomp      = 6
     &                 )
c
c --- local variables
c
c      sigElp   (dp,ar(6  ),l)            trial stress
c      dsdeEl   (dp,ar(6,6),l)            elastic moduli
c      sigDev   (dp,ar(6  ),l)            deviatoric stress tensor
c      dfds     (dp,ar(6  ),l)            derivative of the yield function 
c      JM       (dp,ar(6,6),l)            2D matrix for a 4 order tensor
c      pEl      (dp,sc     ,l)            hydrostatic pressure stress
c      qEl      (dp,sc     ,l)            von-Mises stress
c      pleq_t   (dp,sc     ,l)            equivalent plastic strain at beginnig of time increment
c      pleq     (dp,sc     ,l)            equivalent plastic strain at end of time increment
c      dpleq    (dp,sc     ,l)            incremental equivalent plastic strain
c      sigy_t   (dp,sc     ,l)            yield stress at beginnig of time increments
c      sigy     (dp,sc     ,l)            yield stress at end of time increment
c      young    (dp,sc     ,l)            Young's modulus
c      posn     (dp,sc     ,l)            Poiss's ratio
c      sigy0    (dp,sc     ,l)            initial yield stress
c      dsigdep  (dp,sc     ,l)            plastic slop
c      twoG     (dp,sc     ,l)            two time of shear moduli
c      threeG   (dp,sc     ,l)            three time of shear moduli
c
c --- temperary variables for solution purpose
c      i, j
c      threeOv2qEl, oneOv3G, qElOv3G, con1, con2, fratio
c
      EXTERNAL         vzero, vmove, get_ElmData
      DOUBLE PRECISION sigElp(mcomp), dsdeEl(mcomp,mcomp), G(mcomp),
     &                 sigDev(mcomp), JM    (mcomp,mcomp), dfds(mcomp),
     &                 sigi  (mcomp), strainEl(mcomp)
c
      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7
c
      DATA G/1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,0.0D0/
c
      INTEGER          i, j, k
      DOUBLE PRECISION pEl,   qEl,     pleq_t,  sigy_t , sigy,
     &                 dpleq, pleq, 
     &                 young, posn,    sigy0,   dsigdep, 
     &                 elast1,elast2,
     &                 twoG,  threeG,  oneOv3G, qElOv3G, threeOv2qEl, 
     &                 fratio,  con1,    con2, dperr(3)
c
c --- variaveis locais do usuario necessárias ao cálculo da fluência e retração do 
c concreto
c    
      ! Variaveis de entrada do usuario TBDATA: relacionadas aos dados do concreto
      !
      ! fck [kN/cm²] resistencia característica do concreto 
      ! s [adm] coef. que depende do tipo de concreto 
      ! ti [dias] tempo da analise apartir da qual o concreto iniciou 
      ! rh [%] umidade relativa do ambiente 
      ! hf [cm] espessura fictítica 
      ! ts [dias] idade do concreto no inicio da secagem 
      ! betasc [adm] coeficiente que depende do tipo de cimento 
      ! alpha [adm] efeito do tipo de cimento durante a cura 
      ! timei [dias] Tempo inicial para o ajuste da cadeia de Kelvin 
      ! timef [dias] Tempo final para o ajuste da cadeia de Kelvin
      ! ndec [pontos por decada] No. de pontos para discretização
      ! ligafluencia 0 - desliga a fluência, 1 - liga fluência
      ! ligaretracao 0 - desliga a retração, 1 - liga a retração
      !
      ! Variaveis calculadas na USERMAT
      !
      ! tmat  [dias] idade do material durante a análise (time-ti)
      ! fcm28 [kN/cm²] resistência média do concreto aos 28 dias 
      ! Ec28 [kN/cm²] módulo de elasticidade do concreto aos 28 dias
      ! mkelvin [un] no. máximo de elementos de Kelvin
      ! ndt [un] No. de pontos para calcular os parâmetros da cadeia de Kelvin 
      ! nkelvin [un] No. de camadas da cadeia de Kelvin (DIAS,2013: N da eq3.3, p44) 
      ! tau [cm²/(kN*dia)] tempos de retardação dos elementos da cadeia de Kelvin 
      ! Emu [kN/cm²] módulos de elasticidade dos elementos da cadeia de Kelvin
      ! tmeio [dias] idade do material somada meio incremento de tempo (tmat+dtime/2)
      ! tmatadj [dias] idade do material durante a análise corrigida com a temperatura
      ! fi0 [adm] coeficiente de fluência que depende do tempo t0
      ! V [adm] 
      ! dsdeElinv [cm²/kN] matriz de rigidez inversa
      ! epsloncf1 [adm] primeira parcela deform. fluência de uma unidade de Kelvin 
      ! epsloncf2 [adm] segunda parcela deform. por fluência de uma unidade de Kelvin
      ! epsloncf3 [adm] deformação por fluência da unidade i 
      ! deltay [adm]
      ! lbd [adm]
      !
      ! Outras variaveis
      !
      ! iott, wrinqr variaveis que guardam o numero da output window do Ansys
      !
      ! Variaveis de estado
      !
      ! ustatev(1)    primeirarodada    ! 9999 – indica primeira rodada
      ! ustatev(2)    fcm28                                       
      ! ustatev(3)    Ec28                                        
      ! ustatev(4)    nkelvin                                     
      ! ustatev(5)    timei                                       
      ! ustatev(11..) Emu                                         
      ! ustatev(21..) tau                                         
      ! ustatev(31..) epsloncsant                                
      ! ustatev(41..) epsloncs                                   
      ! ustatev(51..) depsloncs                                   
      ! ustatev(101..) epsloncfantmu                              
      !        (111..) epsloncfantmu
      !        (121..) epsloncfantmu
      !        (130..) epsloncfantmu
      !        (140..) epsloncfantmu
      !        (150..) epsloncfantmu
      !        (160..) epsloncfantmu
      ! ustatev(201..) epsloncfmu                                
      !        (211..) epsloncfmu
      !        (221..) epsloncfmu
      !        (230..) epsloncfmu
      !        (240..) epsloncfmu
      !        (250..) epsloncfmu
      !        (260..) epsloncfmu
      ! ustatev(301..) depsloncf                                  
      ! ustatev(401..) dstress
      !
      DOUBLE PRECISION fck,s,ti,rh,hf,ts,betasc,alpha,timei,timef                                                                                                             
      INTEGER          ndec,ligafluencia,ligaretracao                                                                  
      DOUBLE PRECISION tmat,fcm28,Ec28
      INTEGER          mkelvin,ndt,nkelvin
      PARAMETER        (mkelvin=7)
      DOUBLE PRECISION tau(mkelvin),EMU(mkelvin),tmeio,tmatadj,fi0,V,
     &                 dsdeElinv(mcomp,mcomp),epsloncf1(mcomp),
     &                 epsloncf2(mcomp),epsloncf3(mcomp),
     &                 deltay(mkelvin),lbd(mkelvin)
      EXTERNAL         wrinqr
      INTEGER          iott,wrinqr                                  
      DOUBLE PRECISION wr_output
      DOUBLE PRECISION CEB90_fcm
      DOUBLE PRECISION CEB90_Ec
      DOUBLE PRECISION CEB90_epsloncs
      DOUBLE PRECISION CEB90_tadj
      DOUBLE PRECISION CEB90_tT
      DOUBLE PRECISION CEB90_fi0
c
c*************************************************************************c
      !
      !*************************************************************************!
      ! PASSO1) Entrada de dados do script                                      !
      !*************************************************************************!
      fck             = prop(1)                 
      s               = prop(2)
      ti              = prop(3)
      posn            = prop(4)
      rh              = prop(5)
      hf              = prop(6)
      ts              = prop(7)
      betasc          = prop(8)
      alpha           = prop(9)
      timei           = prop(10)
      timef           = prop(11)
      ndec            = prop(12)
      ligafluencia    = prop(13)
      ligaretracao    = prop(14)
      iott            = wrinqr(WR_OUTPUT)
      !
      !*************************************************************************!
      ! PASSO2) Obtendo parametros do modelo de Kelvin                          !
      !*************************************************************************!
      !
      ! Calcula os parâmetros do modelo de Kelvin apenas uma vez por ponto de integração
      if (ustatev(1).NE.9999) then
          !
          ! Calculos iniciais
          fcm28 = CEB90_fcm(fck,28.0d0,20.0d0,s)
          Ec28 = CEB90_Ec(fck,s,28.0d0,20.0d0)
          !
          ! Calculando numero de elementos de Kelvin e o numero de pontos do ajuste
          call KELPAR(timei,timef,ndec,nkelvin,ndt) 
          !
          ! Calculando os tempos de retardação e módulos do ajuste
          call KELVIN_betactt0(nkelvin,timei,timef,ndec,ndt,
     & tau,Emu,temp,rh,hf)
          !
          ! Guardando os parâmetros
          ustatev(2) = fcm28
          ustatev(3) = Ec28
          ustatev(4) = nkelvin
          ustatev(5) = timei
          !
          do i=1,nkelvin
              ustatev(10+i)=Emu(i)
              ustatev(20+i)=tau(i)
          enddo
          !
          ! Definindo variavel que indica o fim da primeira rodada
          ustatev(1) = 9999
          !
      else
          ! Recupera os parâmetros
          fcm28 = ustatev(2)
          Ec28 = ustatev(3)
          nkelvin = ustatev(4)
          timei = ustatev(5)
          do i=1,nkelvin
              Emu(i) = ustatev(10+i)
              tau(i) = ustatev(20+i)
          enddo           
      endif
      !    
      !*************************************************************************!
      ! PASSO3) Calculo da idade do material                                    !
      !*************************************************************************!
      !
      if (time.GT.ti) then
          tmat = time-ti+timei
      else
          tmat = 0.0d0
      endif
      ustatev(6) = tmat
      !
      !*************************************************************************!
      ! PASSO4) Calculo do incremento de deformacao por retracao                !
      !*************************************************************************!
      !
      if (time.GE.ti) then
          !
          ! Calcula o coeficiente de retração para o próximo passo de tempo (tmat+dtime)
          do i=1,3
              ustatev(40+i) = 
     & CEB90_epsloncs(tmat+dtime,ts,fcm28,betasc,rh,hf,temp)
          enddo
          !
          ! Calcula o incremento de deformação por retração (apenas nas três direções principais)
          do i=1,3
              ustatev(50+i) = 
     & ustatev(40+i)-ustatev(30+i)
          enddo
      endif
      !
      !*************************************************************************!
      ! PASSO5) Calculo do incremento de deformacao por fluencia                !
      !*************************************************************************!
      !
      if (time.GE.ti) then  
          !
          ! Guardando deformação por fluência do ultimo passo que convergiu
          do i=1,nkelvin
              do j=1,ncomp
                  ustatev(100+(i-1)*10+j) = ustatev(200+(i-1)*10+j)
              enddo
          enddo             
          !
          ! Calculando o fator do envelhecimento
          tmatadj = CEB90_tadj(tmat,alpha,temp,tmat)
          tmeio = tmatadj+dtime/2.0d0
          fi0 = CEB90_fi0(tmeio,rh,hf,fcm28,temp)
          V = Ec28/fi0
          !
          ! Inicializando o incremento de deformação por fluência
          do i =1,ncomp
              ustatev(300+i) = 0.0d0
          enddo
          !
          ! Calculo do incremento de deformação por fluência
          do i=1,nkelvin
              !
              ! Calculo da Matriz Dmu
              deltay(i) = dtime/tau(i)
              lbd(i) = (1.0d0-dexp(-deltay(i)))/deltay(i)
              young = (Emu(i)*V)/(lbd(i))
              call MatrizD(young,posn,dsdeEl)                   
              call matinv(6,dsdeEl,dsdeElinv)
              !
              ! Calculo da primeira parcela da deformação por fluência de uma unidade
              epsloncf1=0.0d0
              do j=1,ncomp
                  do k=1,ncomp
                      epsloncf1(j) = epsloncf1(j)+
     & dsdeElinv(j,k)*ustatev(400+k)
                  enddo
              enddo
              !
              ! Calculo da segunda parcela da deformação por fluência de uma unidade 
              do j =1,ncomp
                  epsloncf2(j) = 
     & ustatev(100+(i-1)*10+j)*dexp(-deltay(i))
              enddo
              !
              ! Calculo da deformação por fluência da unidade i
              do j=1,ncomp
                  ustatev(200+(i-1)*10+j) =
     & epsloncf1(j)+epsloncf2(j)
              enddo
              !                
              ! Calculo do incremento de deformação por fluência de todas unidades
              do j =1,ncomp
                  ustatev(300+j) =
     & ustatev(300+j)+
     & ustatev(200+(i-1)*10+j)*(1-dexp(-deltay(i)))
              enddo
              !                                    
          enddo
          !
      endif
      !
      !*************************************************************************!
      ! PASSO6) Calculo do incremento de deformacao elastica                    !
      !*************************************************************************!
      !
      ! calculo do deltay
      deltay = 0.0d0
      do i = 1,nkelvin
          deltay(i) = dtime/tau(i)
      enddo
      !
      ! calculo do lbd
      lbd = 0.0d0
      do i = 1,nkelvin
          lbd(i) = (1.0d0-dexp(-deltay(i)))/deltay(i)
      enddo
      !
      ! calculo do v
      tmatadj = CEB90_tadj(tmat,alpha,temp,tmat)
      tmeio = tmatadj+dtime/2.0d0
      fi0 = CEB90_fi0(tmeio,rh,hf,fcm28,temp)
      V = Ec28/fi0
      !
      ! calculo do modulo de Young corrigido
      young = 0.0d0
      do i=1,nkelvin
          young = young + (1.0d0-lbd(i))/(Emu(i)*V)*ligafluencia  
      enddo
      young = young + 1.0d0/CEB90_Ec(fck,s,tmeio,temp)
      young = 1.0d0/young
      !
      ! Calculando a matriz de rigidez elástica
      call MatrizD(young,posn,dsdeEl)
      !
      ! Descontando o incremento de deformação por fluência e retração do incremento de deformação total (Equação 5.60)
      do i=1,ncomp
          dstrain(i) = dstrain(i)
     & -ustatev(300+i)*ligafluencia 
     & -ustatev(50+i)*ligaretracao    
      end do
      !
      !*************************************************************************!
      ! PASSO7) Calculo do incremento de tensão                                 !
      !*************************************************************************!
      !     
      sigElp = 0.0d0
      do i=1,ncomp
          do j=1,ncomp
              sigElp(i) = sigElp(i)+dsdeEl(i,j)*(dstrain(j))
          end do
      end do 
      !
      !
      !*************************************************************************!
      ! PASSO8) Atualizando a tensão                                            !
      !*************************************************************************!
      !      
      do i=1,ncomp
         stress(i) = stress(i)+sigElp(i)
      end do 
      !
      !*************************************************************************!
      ! PASSO9) Atualizando a matriz Jacobiana                                  !
      !*************************************************************************!
      !
      do i=1,ncomp
          do j=1,ncomp
              dsdePl(j,i) = dsdeEl(j,i)               
          end do
      end do 
      !
      !*************************************************************************!
      ! PASSO10) Guardando incremento de tensao e deformação por retração       !
      !*************************************************************************!
      !     
      do i=1,ncomp
          ustatev(400+i) = sigElp(i)
      enddo
      !
      do i=1,ncomp
          ustatev(30+i) = ustatev(40+i)     
      enddo
      !
      return
  100 FORMAT(1X,100(E10.4,'   '))
   90 FORMAT(E14.8)
   80 FORMAT(f10.2,f10.2,f8.2,f8.2,f8.2)
      End

      subroutine matrizD(E,Poisson,D)
      !**************************************************************************!
      !** Função: matrizD                                                      **!
      !**                                                                      **!
      !** Objetivo: calcula a matriz consitutiva do material isotrópico        **!
      !**                                                                      **!
      !** Situação: (28-09-2016) OK                                            **!
      !**                                                                      **!
      !**************************************************************************!
      implicit none
      DOUBLE PRECISION E         ! módulo de elasticidade                  
      DOUBLE PRECISION Poisson    ! coeficiente de Poisson
      DOUBLE PRECISION D(6,6)     ! matriz constitutiva elástica isotrópica
      !
	D=0.0d0
	D(1,1)=(E*(1.0d0-Poisson))/((1.0d0+Poisson)*(1.0d0-2.0d0*Poisson))
	D(1,2)=(E*Poisson)/((1.0d0+Poisson)*(1.0d0-2.0d0*Poisson))
	D(1,3)=D(1,2)
	D(2,1)=D(1,2)
	D(2,2)=D(1,1)
	D(2,3)=D(1,2)
	D(3,1)=D(1,3)
	D(3,2)=D(2,3)
	D(3,3)=D(1,1)
	D(4,4)=(E)/((1.0d0+Poisson)*2.0d0)
	D(5,5)=D(4,4)
	D(6,6)=D(4,4)
      !
      end subroutine MatrizD
 
      function CEB90_tT(temp,deltat)
      !**************************************************************************!
      !** Função: CEB90_tT                                                     **!
      !**                                                                      **!
      !** Objetivo: considera o efeito da temperatura na maturidade do         **!
      !**           concreto                                                   **!
      !**                                                                      **!
      !** Referencia: CEB-MC90                                                 **!
      !** Situação: (28-09-2016) OK                                            **!
      !**                                                                      **!
      !**************************************************************************!
      implicit none
      DOUBLE PRECISION temp   ! [oC] temperatura       
      DOUBLE PRECISION deltat ! [dias] intervalo em dias em que a temp. prevalece     
      DOUBLE PRECISION CEB90_tT
      !
      CEB90_tT = deltat*dexp(13.65d0-4000.0d0/(273.0d0+temp))
      !
      end
      
      function CEB90_tadj(t,alpha,temp,deltat)
      !**************************************************************************!
      !** Função: CEB90_tadj                                                   **!
      !**                                                                      **!
      !** Objetivo: considera o efeito do tipo de cimento e temperatura de     **!
      !**           cura na idade de carregamento.                             **!
      !**                                                                      **!
      !** Referencia: CEB-MC90                                                 **!
      !** Situação: (28-09-2016) OK                                            **!
      !**                                                                      **!
      !**************************************************************************!
      implicit none
      DOUBLE PRECISION t      ! [dias] idade         
      DOUBLE PRECISION alpha  ! [adm] coeficiente que depende do tipo de cimento    
      DOUBLE PRECISION temp   ! [oC] temperatura       
      DOUBLE PRECISION deltat ! [dias] intervalo em dias em que a temp. prevalece 
      DOUBLE PRECISION tT     ! [dias] idade ajustada     
      DOUBLE PRECISION CEB90_tadj
      DOUBLE PRECISION CEB90_tT
      !
      tT = CEB90_tT(temp,deltat)
      CEB90_tadj = tT * (9.0d0 / (2.0d0 + (tT) ** 1.2d0) + 1.0d0) 
     & ** alpha
      If (CEB90_tadj.LE.0.5d0) Then 
          CEB90_tadj = 0.5d0
      endif
      !
      end
      function CEB90_fcm(fck,t,temp,s)
      !**************************************************************************!
      !** Função: CEB90_fcm                                                    **!
      !**                                                                      **!
      !** Objetivo: calcular a resistência média do concreto considerando o    **!
      !**           o tempo t e a temperatura constante temp.                  **!
      !**                                                                      **!
      !** Referencia: CEB-MC90                                                 **!
      !** Situação: (28-09-2016) OK                                            **!
      !**                                                                      **!
      !**************************************************************************!
      implicit none
      DOUBLE PRECISION fck     ! [kN/cm²] resistência característica a compressão    
      DOUBLE PRECISION t       ! [dias] idade
      DOUBLE PRECISION s       ! [adm] coeficiente que depende do tipo de cimento
      DOUBLE PRECISION temp    ! [oC] temperatura
      DOUBLE PRECISION Bcc     ! [adm] coeficiente
      DOUBLE PRECISION fcm     ! [kN/cm²] resistência média a compressão
      DOUBLE PRECISION CEB90_fcm
      !
      fcm = fck + 0.8d0 
      Bcc = dexp(s*(1.0d0-dsqrt(28.0d0/t)))  
      fcm = Bcc*fcm
      fcm = fcm*(1.06d0-0.003d0*temp)
      CEB90_fcm = fcm
      !
      End
           
      function CEB90_Ec(fck,s,t,temp)  
      !**************************************************************************!
      !** Função: CEB90_Eci                                                    **!
      !**                                                                      **!
      !** Objetivo: calcular o módulo de elasticidade do concreto considerando **!
      !**           o tempo t e a temperatura constante temp.                  **!
      !**                                                                      **!
      !** Referencia: CEB-MC90                                                 **!
      !** Situação: (28-09-2016) OK                                            **!
      !**                                                                      **!
      !**************************************************************************!
      implicit none
      DOUBLE PRECISION fck        ! [kN/cm²] resistência característica a compressão    
      DOUBLE PRECISION s          ! [adm] coeficiente que depende do tipo de cimento       
      DOUBLE PRECISION t          ! [dias] idade      
      DOUBLE PRECISION temp       ! [oC] temperatura       
      DOUBLE PRECISION Ec0        ! [kN/cm²] módulo de elasticidade de referência  
      DOUBLE PRECISION fcm0       ! [kN/cm²] resistência de referência
      DOUBLE PRECISION Eci        ! [kN/cm²] módulo de elasticidade inicial
      DOUBLE PRECISION Bcc        ! [adm] coeficiente
      DOUBLE PRECISION Be         ! [adm] coeficiente
      DOUBLE PRECISION Ec         ! [kN/cm²] módulo de elasticidade do concreto
      DOUBLE PRECISION CEB90_Ec
      !
      Ec0 = 21500.0d0/10.0d0                                    
      fcm0 = 1.0d0            
      Eci = Ec0*((fck+0.8d0)/fcm0)**(1.0d0/3.0d0)                   
      Eci = Eci*(1.06d0-0.003d0*temp)
      Bcc = dexp(s*(1.0d0-dsqrt(28.0d0/t)))               
      Be = SQRT(Bcc)                             
      Ec = Be*Eci 
      !
      CEB90_Ec = Ec
      !
      end 

      function CEB90_fi0(t0,rh,hf,fcm,temp)
      !**************************************************************************!
      !** Função: CEB90_fi0                                                    **!
      !**                                                                      **!
      !** Objetivo: calcular o coeficiente de fluência que depende do tempo t0 **!
      !**                                                                      **!
      !** Referencia: CEB-MC90                                                 **!
      !** Situação: (28-09-2016) OK                                            **!
      !**                                                                      **!
      !**************************************************************************!
      implicit none
      DOUBLE PRECISION t0             ! [dias] idade no instante de aplicação da carga         
      DOUBLE PRECISION rh             ! [%] umidade relativa
      DOUBLE PRECISION hf             ! [cm] espessura ficticia
      DOUBLE PRECISION fcm            ! [kN/cm²] resistência média a compressão
      DOUBLE PRECISION temp           ! [oC] temperatura
      DOUBLE PRECISION fcm0           ! [kN/cm²^] resistência de referencia
      DOUBLE PRECISION firh           ! [adm] coeficiente
      DOUBLE PRECISION fit            ! [adm] coeficiente
      DOUBLE PRECISION betafcm        ! [adm] coeficiente
      DOUBLE PRECISION betat0         ! [adm] coeficiente
      DOUBLE PRECISION CEB90_fi0
      !
      fcm0 = 1.0d0                                                      
      firh = 1.0d0+(1.0d0-RH/100.0d0)/(0.46d0*(hf*10.0d0/100.0d0)**
     & (1.0d0/3.0d0))     
      fit = dexp(0.015d0*(temp-20.0d0))
      firh = fit+(firh-1.0d0)*fit**1.2d0
      betafcm = 5.3d0/((fcm/fcm0)**0.5d0)
      betat0 = 1.0d0/(0.1d0+(t0/1.0d0)**0.2d0)
      CEB90_fi0 = firh*betafcm*betat0
      !
      end
      
      function CEB90_betactt0(t,t0,rh,hf,temp)
      !**************************************************************************!
      !** Função: CEB90_betactt0                                               **!
      !**                                                                      **!
      !** Objetivo: calcular o coeficiente de fluência que depende do          **!
      !**           intervalo de tempo de aplicação da carga t-t0              **!      
      !**                                                                      **!
      !** Referencia: CEB-MC90                                                 **!
      !** Situação: (28-09-2016) OK                                            **!
      !**                                                                      **!
      !**************************************************************************!
      implicit none
      !
      DOUBLE PRECISION t                  ! [dias] idade      
      DOUBLE PRECISION t0                 ! [dias] idade na aplicação da carga
      DOUBLE PRECISION rh                 ! [%] umidade relativa
      DOUBLE PRECISION hf                 ! [cm] espessura fictícia
      DOUBLE PRECISION temp               ! [oC] temperatura
      DOUBLE PRECISION betah              ! [adm] coeficiente
      DOUBLE PRECISION betat              ! [adm] coeficiente
      DOUBLE PRECISION CEB90_betactt0
      !
      betah = 150.0d0*(1.0d0+(1.2d0*rh/100.0d0)**18.0d0)*(hf*10.0d0)
     & /100.0d0+250.0d0
      betat = exp(1500.0d0/(273.0d0+temp)-5.12d0)
      betah = betah*betat
      if (betah.GT.1500.0d0) then
          betah = 1500.0d0
      endif
      CEB90_betactt0 = ((t-t0)/(betah+(t-t0)))**0.3d0
      !
      end
      
      function CEB90_epsloncs(t,ts,fcm,betasc,rh,hf,temp)
      !**************************************************************************!
      !** Função: CEB90_epsloncs                                               **!
      !**                                                                      **!
      !** Objetivo: calcular o deformação de retração                          **!    
      !**                                                                      **!
      !** Referencia: CEB-MC90                                                 **!
      !** Situação: (28-09-2016) OK                                            **!
      !**                                                                      **!
      !**************************************************************************!
      implicit none
      DOUBLE PRECISION t           ! [dias] idade    
      DOUBLE PRECISION ts          ! [dias] idade no final da cura
      DOUBLE PRECISION fcm         ! [kN/cm²] resistência média a compressão
      DOUBLE PRECISION betasc      ! [adm] coeficien que depende do tipo de cimento
      DOUBLE PRECISION rh          ! [%] umidade relativa
      DOUBLE PRECISION hf          ! [cm] espessura ficticia
      DOUBLE PRECISION temp        ! [oC] temperatura
      DOUBLE PRECISION fcm0        ! [kN/cm²] resistência de referencia
      DOUBLE PRECISION betasrh     ! [adm] coeficiente
      DOUBLE PRECISION betarh      ! [adm] coeficiente
      DOUBLE PRECISION betast      ! [adm] coeficiente
      DOUBLE PRECISION betastts    ! [adm] coeficiente
      DOUBLE PRECISION epsloncsfcm  ! [adm] deformação específica
      DOUBLE PRECISION epsloncs0    ! [adm] deformação específica
      DOUBLE PRECISION CEB90_epsloncs
      !
      if (t.GT.ts) then
          fcm0 = 1.0d0    
          epsloncsfcm = (160.0d0+10.0d0*betasc*(9.0d0-fcm/fcm0))*10.0d0
     % **(-6)
          if (rh.GE.40.0d0.AND.rh.LE.99.0d0) THEN
              betasrh = 1.0d0-(rh/100.0d0)**3
              betarh = -1.55d0*betasrh    
          else
              betarh = 0.25d0
          endif
          betast = 1.0d0+(8.0d0/(103.0d0-rh))*((temp-20.0d0)/40.0d0)
          betarh = betarh*betast
          epsloncs0 = epsloncsfcm*betarh
          betastts = ((t-ts)/(350.0d0*(hf*10.0d0/100.0d0)**2*
     & exp(-0.06d0*(temp-20.0d0))+(t-ts)))**0.5d0
          CEB90_epsloncs = epsloncs0*betastts
      else
          CEB90_epsloncs = 0.0d0
      endif 
      !
      end
      
      subroutine KELPAR(timei,timef,ndec,nkelvin,ndt)
      !**************************************************************************!
      !** Subrotina: KELPAR                                                    **!
      !**                                                                      **!
      !** Objetivo: subrotina para determinar o numero de camadas de kelvin    **!
      !**           e o numero de intervalos de tempo                          **!
      !**                                                                      **!
      !** Referencia: DIAS,M.M. Análise numérica de vigas mistas aço-concreto  **!
      !**             pelo método dos elementos finitos: efeitos de longa      **!
      !**             duração. 2013. Dissertação (Mestrado em Engenharia)      **!
      !**             Programa de Pós Graduação em Engenharia Civil, UFRGS,    **!
      !**             Porto Alegre.                                            **!
      !**                                                                      **!
      !**             BAZANT,Z.P.; PRASANNAN, S. Solidification Theory for     **!
      !**             concrete creep. I: Formulation, Journal of Enginnering   **!
      !**             mechanics, v.115, n.8, p.1691-1703, 1989.                **!
      !**                                                                      **!      
      !**             BAZANT,Z.P.; PRASANNAN, S. Solidification Theory for     **!
      !**             concrete creep. II: verification and application.        **!
      !**             Journal of Enginnering Mechanics, v.115, n.8,            **!
      !**             p.1704-1725, 1989.                                       **!
      !**                                                                      **!
      !** Situação: (28-09-2016) OK                                            **!
      !**                                                                      **!
      !**************************************************************************!
      implicit none
      !
      DOUBLE PRECISION timei          ! [dias] idade inicial para o ajuste      
      DOUBLE PRECISION timef          ! [dias] idade final para o ajuste
      INTEGER          ndec           ! [un] numero de pontos por decada no ajuste
      INTEGER          nkelvin        ! [un] numero de elementos de Kelvin
      INTEGER          ndt            ! [un] total de pontos do ajuste     
      DOUBLE PRECISION d4             ! auxiliar  
      DOUBLE PRECISION d2             ! auxiliar
      !
      ! numero de camadas do modelo de Kelvin Generalizado
      nkelvin=0
	do
		nkelvin=nkelvin+1
		d4=(10.0d0**(nkelvin-1))*(0.01d0*timei)
		if(d4>=(timef/2.0d0)) then
			exit
		end if
	end do
      if(nkelvin>=7) then
		nkelvin=6
      end if
      !
      ! numero de tempos discretos
	ndt=0
	d2=timei*0.1d0
	do 
		ndt=ndt+1
		d2=(10.0d0**(1.0d0/ndec))*d2
		if(d2>=timef-timei) then
			exit
		end if
      end do
      !
      end
           


      subroutine KELVIN_betactt0(nkelvin,timei,timef,ndec,ndt,
     & tau,Emu,temp,rh,hf)
      !**************************************************************************!
      !** Subrotina: Kelvin_Generalizado                                       **!
      !**                                                                      **!
      !** Objetivo: calcula os parâmetros do modelo de Kelvin-Generalizado     **!
      !**                                                                      **!
      !** Referencia: DIAS,M.M. Análise numérica de vigas mistas aço-concreto  **!
      !**             pelo método dos elementos finitos: efeitos de longa      **!
      !**             duração. 2013. Dissertação (Mestrado em Engenharia)      **!
      !**             Programa de Pós Graduação em Engenharia Civil, UFRGS,    **!
      !**             Porto Alegre.                                            **!
      !**             Porto Alegre.                                            **!
      !**                                                                      **!
      !**             BAZANT,Z.P.; PRASANNAN, S. Solidification Theory for     **!
      !**             concrete creep. I: Formulation, Journal of Enginnering   **!
      !**             mechanics, v.115, n.8, p.1691-1703, 1989.                **!
      !**                                                                      **!      
      !**             BAZANT,Z.P.; PRASANNAN, S. Solidification Theory for     **!
      !**             concrete creep. II: verification and application.        **!
      !**             Journal of Enginnering Mechanics, v.115, n.8,            **!
      !**             p.1704-1725, 1989.                                       **!
      !**                                                                      **! 
      !**             CEB-FIP MODEL CODE 1990                                  **!
      !**                                                                      **!
      !** Situação: (29-09-2016) TESTE                                         **!
      !**                                                                      **!
      !**************************************************************************!
      implicit none
      !
      INTEGER          nkelvin        ! [un] numero de elementos de Kelvin           
      DOUBLE PRECISION timei          ! [dias] idade inicial para o ajuste               
      DOUBLE PRECISION timef          ! [dias] idade final para o ajuste             
      INTEGER          ndec           ! [un] numero de pontos por decada no ajuste          
      INTEGER          ndt            ! [un] total de pontos do ajuste             
      DOUBLE PRECISION tau(nkelvin)   ! [cm²/(kN*dia)] tempos de retardação da cadeia    
      DOUBLE PRECISION Emu(nkelvin)   ! [kN/cm²] modulos de elasticidade da cadeia
      DOUBLE PRECISION temp           ! [oC] temperatura
      DOUBLE PRECISION rh             ! [%] umidade relativa
      DOUBLE PRECISION hf             ! [cm] espessura ficticia           
      !
      DOUBLE PRECISION dti(ndt)       ! [dias] intervalos de tempo          
      DOUBLE PRECISION fit            ! [adm] guarda o coeficiente que será ajustado
      DOUBLE PRECISION C(nkelvin,nkelvin) ! matriz de coeficiente 
      DOUBLE PRECISION b(nkelvin)     ! vetor de termos independentes
      DOUBLE PRECISION x(nkelvin)     ! vetor de incognitas
      INTEGER          ntau           ! auxiliar
      INTEGER          ndti           ! auxiliar
      INTEGER          ndt0i          ! auxiliar
      INTEGER          i,j,k          ! contador
      DOUBLE PRECISION a1,a2,aa       ! auxiliar
      INTEGER          er             ! variavel que indica erro na solução do sistema
      DOUBLE PRECISION CEB90_betactt0
      !
      ! Inicializa a variavel
      tau = 0.0d0
      dti = 0.0d0
      !
      ! Zerando a matriz C e o vetor b
      C = 0.0d0
      b = 0.0d0
      x = 0.0d0
      !
      ! Calcula os intervalos de tempo de aplicação de carga 
      dti(1)=0.1d0*timei
      DO ndti=2,ndt
          dti(ndti)=(10.0d0**(1.0d0/ndec))*dti(ndti-1)
      ENDDO
      dti(ndt)=timef-timei
      !
      ! Calcula os tempos de retaração
      tau(1) = 0.01d0*timei          
      DO ntau=2,nkelvin
          tau(ntau)=(10.0d0**(ntau-1))*tau(1)
      ENDDO
      !
      ! Calculo da matriz C
      DO i=1,nkelvin
          DO j=1,nkelvin
              DO k=1,ndt
                  a1=1.0d0-DEXP(-dti(k)/tau(i))
                  a2=1.0d0-DEXP(-dti(k)/tau(j))
                  aa=a1*a2
                  c(i,j)=c(i,j)+aa
              ENDDO
          ENDDO
      ENDDO
      !
      DO i=1,nkelvin
          DO k=1,ndt
              a1=1.0d0-DEXP(-dti(k)/tau(i))
              fit = CEB90_betactt0(dti(k),0.0d0,rh,hf,temp)
              a2= fit
              aa=a1*a2
              b(i)=b(i)+aa
          ENDDO
      ENDDO
      !
      !
      ! Resolução do sistema
      call gauss_parcial(C,b,nkelvin,x,1d-20,er)
      !
      ! Módulos de elasticidade do modelo de Kelvin-Generalizado
      DO i=1,nkelvin
          Emu(i) = 1.0d0/x(i)      
      ENDDO
      !
      !
  100 FORMAT(1X,100(E10.4,'   '))
   90 FORMAT(E14.4)
      END     
                 
      subroutine gauss_parcial(a,b,n,x,tol,er)
      !**************************************************************************!
      !** Subrotina: gauss_parcial                                             **!
      !**                                                                      **!
      !** Objetivo: resolve sistema de equações pelo metodo de Gauss.          **!
      !**                                                                      **!
      !** Referencia: CHAPRA, S. C.; CANALE, R. P. Métodos numéricos para      **!
      !**             engenharia. 5.ed, São Paulo: McGraw-Hill, 2008.          **!      
      !**                                                                      **!
      !** Situação: (26-10-2016) OK                                            **!
      !**                                                                      **!
      !**************************************************************************!
      implicit none
      !      
      DOUBLE PRECISION a(n,n)     ! matriz dos coeficientes
      DOUBLE PRECISION b(n)       ! vetor de termos independentes
      INTEGER          n          ! dimensão do sistema
      DOUBLE PRECISION x(n)       ! vetor solução
      DOUBLE PRECISION tol        ! tolerância
      INTEGER          er         ! -1, indica que houve erro
      DOUBLE PRECISION s(n)       ! vetor que guarda o maior valor da linha
      INTEGER          i,j        ! contador
      !
      er = 0
      !
      ! Identifica o maior valor da linha
      DO i=1,n
          s(i) = abs(a(i,1))
          DO j=2,n
              IF (abs(a(i,j)).GT.s(i)) THEN
                  s(i) = abs(a(i,j))
              ENDIF
          ENDDO
      ENDDO
      !
      ! Inicia a eliminação de Gauss
      CALL eliminate(a,s,n,b,tol,er)
      !
      ! Faz a retrosubstituição
      IF (er.NE.-1) THEN
          call substitute(a,n,b,x)
      ENDIF
      !
      end
      
      subroutine eliminate(a,s,n,b,tol,er)
      !**************************************************************************!
      !** Subrotina: eliminate                                                 **!
      !**                                                                      **!
      !** Objetivo: faz a eliminação do sistema de Gauss.                      **!
      !**                                                                      **!
      !** Referencia: CHAPRA, S. C.; CANALE, R. P. Métodos numéricos para      **!
      !**             engenharia. 5.ed, São Paulo: McGraw-Hill, 2008.          **!      
      !**                                                                      **!
      !** Situação: (26-10-2016) OK                                            **!
      !**                                                                      **!
      !**************************************************************************!
      implicit none
      !
      DOUBLE PRECISION a(n,n) ! matriz dos coeficientes
      DOUBLE PRECISION s(n)   ! vetor que guarda o maior valor da linha
      INTEGER          n      ! dimensão do sistema
      DOUBLE PRECISION b(n)   ! vetor de termos independentes  
      DOUBLE PRECISION x(n)   ! vetor solução
      DOUBLE PRECISION tol    ! tolerância
      INTEGER          er     ! -1, indica que houve erro
      INTEGER          i,j,k  ! contador
      DOUBLE PRECISION fator
      !
      ! Iniciando a eliminação
      DO k=1,n-1
          !
          ! Verificação do pivot
          call pivot(a,b,s,n,k)
          !
          ! Testando pivot zero ou próximo de zero
          if (abs(a(k,k)/s(k)).LT.tol) then
              er = -1
              return
          endif
          !
          ! Divide a linha pelo fator
          do i=k+1,n
              fator = a(i,k)/a(k,k)
              do j=k+1,n
                  a(i,j) = a(i,j)-fator*a(k,j)
              enddo
              b(i) = b(i)-fator*b(k)
          enddo
      enddo
      !
      ! Testando pivot zero ou próximo de zero
      if (abs(a(k,k)/s(k)).LT.tol) then
          er=-1
      endif
      
      end
      
      subroutine pivot(a,b,s,n,k)
      !**************************************************************************!
      !** Subrotina: pivot                                                     **!
      !**                                                                      **!
      !** Objetivo: faz a troca de linhas para evitar pivot negativo ou zero   **!
      !**                                                                      **!
      !** Referencia: CHAPRA, S. C.; CANALE, R. P. Métodos numéricos para      **!
      !**             engenharia. 5.ed, São Paulo: McGraw-Hill, 2008.          **!      
      !**                                                                      **!
      !** Situação: (26-10-2016) OK                                            **!
      !**                                                                      **!
      !**************************************************************************!      
      implicit none
      !
      DOUBLE PRECISION a(n,n) ! matriz dos coeficientes
      DOUBLE PRECISION b(n)   ! vetor de termos independentes  
      DOUBLE PRECISION x(n)   ! vetor solução
      DOUBLE PRECISION s(n)   ! vetor que guarda o maior valor da linha 
      INTEGER          n      ! dimensão do sistema
      INTEGER          k,p    ! auxiliar
      INTEGER          ii,jj  ! contador
      DOUBLE PRECISION dummy  ! auxiliar
      DOUBLE PRECISION maior  ! auxiliar
      !
      ! Identificando o maior da linha
      p=k
      maior = abs(a(k,k)/s(k))
      DO ii=k+1,n
          dummy = abs(a(ii,k)/s(ii))
          if (dummy.GT.maior) then
              maior = dummy
              p=ii
          endif
      enddo
      !
      ! Trocando a linha para ter o maior pivot
      if (p.NE.k) then
          do jj=k,n
              dummy = a(p,jj)
              a(p,jj) = a(k,jj)
              a(k,jj) = dummy
          enddo
          dummy = b(p)
          b(p)=b(k)
          b(k) = dummy
          dummy = s(p)
          s(p) = s(k)
          s(k) = dummy
      endif
      !
      End
      
      subroutine substitute(a,n,b,x)
      !**************************************************************************!
      !** Subrotina: substitute                                                **!
      !**                                                                      **!
      !** Objetivo: faz a retrosubstituição do método de Gauss                 **!
      !**                                                                      **!
      !** Referencia: CHAPRA, S. C.; CANALE, R. P. Métodos numéricos para      **!
      !**             engenharia. 5.ed, São Paulo: McGraw-Hill, 2008.          **!      
      !**                                                                      **!
      !** Situação: (26-10-2016) OK                                            **!
      !**                                                                      **!
      !**************************************************************************!      
      implicit none
      !      
      DOUBLE PRECISION a(n,n) ! matriz dos coeficientes
      INTEGER          n      ! dimensão do sistema
      DOUBLE PRECISION b(n)   ! vetor de termos independentes  
      DOUBLE PRECISION x(n)   ! vetor solução
      INTEGER          i,j    ! contador
      DOUBLE PRECISION soma   ! auxiliar
      !
      ! Fazendo retrosubstituicao
      x(n) = b(n)/a(n,n)
      do i=n-1,1,-1
          soma = 0.0d0
          do j = i+1,n
              soma = soma+a(i,j)*x(j)
          enddo
          x(i) = (b(i)-soma)/a(i,i)
      enddo
      !
      end
      


      subroutine matinv(n,a,ainv)
      !**************************************************************************!
      !** Subrotina: matinv                                                    **!
      !**                                                                      **!
      !** Objetivo: inverte uma matriz pela técnica de pivotamento             **!
      !**                                                                      **!
      !** Situação: (26-10-2016) OK                                            **!
      !**                                                                      **!
      !**************************************************************************!      
      implicit none
      !      
      INTEGER          n          ! dimensão do sistema
      DOUBLE PRECISION a(n,n)     ! matriz dos coeficientes
      DOUBLE PRECISION ainv(n,n)  ! matriz inversa  
      DOUBLE PRECISION b(n,2*n)   ! matriz aumentada
      DOUBLE PRECISION pivot      ! pivô
      DOUBLE PRECISION xnum       ! auxiliar
      INTEGER          i,j,k      ! contador
      !
      ! Fazer matriz aumentada
      do i=1,n
          do j=1,n
              b(i,j) = 0.0d0
              b(i,j+n) = 0.0d0
              b(i,j)=a(i,j)
              if(i.eq.j)then
                  b(i,j+n)=1.0d0
              endif
          enddo
      enddo
      !
      do i=1,n
          ! Escolher o elemento não nulo mais a esquerda como pivot
          do j=1,n
              if (dabs(b(i,j)).gt.0.0d0)then
                  pivot=b(i,j)
                  exit
              endif
          enddo
          !
          ! Passo 1: alterar o pivo escolhido para 1 dividindo a linha do pivo pelo valor do pivo
          do j=1,2*n
              b(i,j)=b(i,j)/pivot
          enddo
          pivot=b(i,i)
          !
          ! Passo 2: mudando o restante da coluno do pivo para 0, adicionando a cada linha um multiplo adequado do pivot
          do k=1,n
              if(k.ne.i)then
                  xnum=b(k,i)/pivot
                  do j=1,2*n
                      b(k,j)=b(k,j)-xnum*b(i,j)
                  enddo
              endif
          enddo
      enddo
      !
      ! Prepara a matriz inversa final
      do i=1,n
          do j=1,n
              ainv(i,j)=b(i,j+n)
          enddo
      enddo
      return        
      end
    
     
      
      
      
*deck,usermatbm    USERDISTRIB  parallel                                gal
      subroutine usermatbm(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7, var8)
c*************************************************************************
c     *** primary function ***
c
c           user defined material constitutive model
c
c      Attention:
c           User must define material constitutive law properly
c           according to the stress state such as 3D, plane strain
c           and axisymmetry, plane stress and beam.
c
c           a 3D material constitutive model can use for
c           plane strain and axisymmetry cases.
c
c           When using shell elements, a plane stress algorithm
c           must be use.
c
c                                             gal July, 1999
c
c       The following demonstrates a USERMAT subroutine for 
c       a plasticity model in 3D beam(188, 189). The plasticity
c       model is the same as TB, BISO.
c       See "ANSYS user material subroutine USERMAT" for detailed
c       description of how to write a USERMAT routine.
c
c*************************************************************************
c
c     input arguments
c     ===============
c      matId     (int,sc,i)               material #
c      elemId    (int,sc,i)               element #
c      kDomIntPt (int,sc,i)               "k"th domain integration point
c      kLayer    (int,sc,i)               "k"th layer
c      kSectPt   (int,sc,i)               "k"th Section point
c      ldstep    (int,sc,i)               load step number
c      isubst    (int,sc,i)               substep number
c      nDirect   (int,sc,in)              # of direct components
c      nShear    (int,sc,in)              # of shear components
c      ncomp     (int,sc,in)              nDirect + nShear
c      nStatev   (int,sc,l)               Number of state variables
c      nProp     (int,sc,l)               Number of material ocnstants
c
c      Temp      (dp,sc,in)               temperature at beginning of
c                                         time increment
c      dTemp     (dp,sc,in)               temperature increment 
c      Time      (dp,sc,in)               time at beginning of increment (t)
c      dTime     (dp,sc,in)               current time increment (dt)
c
c      Strain   (dp,ar(ncomp),i)          Strain at beginning of time increment
c      dStrain  (dp,ar(ncomp),i)          Strain increment
c      prop     (dp,ar(nprop),i)          Material constants defined by TB,USER
c      coords   (dp,ar(3),i)              current coordinates
c      defGrad_t(dp,ar(3,3),i)            Deformation gradient at time t
c      defGrad  (dp,ar(3,3),i)            Deformation gradient at time t+dt
c
c     input output arguments              
c     ======================             
c      stress   (dp,ar(nTesn),io)         stress
c      ustatev   (dp,ar(nStatev),io)       statev
c           ustatev(1)                     - equivalent plastic strain
c           ustatev(2) - ustatev(1+ncomp)  - plastic strain vector
c           ustatev(nStatev)               - von-Mises stress
c      sedEl    (dp,sc,io)                elastic work
c      sedPl    (dp,sc,io)                plastic work
c      epseq    (dp,sc,io)                equivalent plastic strain
c      tsstif   (dp,ar(2),io)             transverse shear stiffness
c                                         tsstif(1) - Gxz
c                                         tsstif(2) - Gyz
c                                         tsstif(1) is also used to calculate hourglass
c                                         stiffness, this value must be defined when low
c                                         order element, such as 181, 182, 185 with uniform 
c                                         integration is used.
c      var?     (dp,sc,io)                not used, they are reserved arguments 
c                                         for further development
c
c     output arguments
c     ================
c      keycut   (int,sc,io)               loading bisect/cut control
c                                         0 - no bisect/cut
c                                         1 - bisect/cut 
c                                         (factor will be determined by ANSYS solution control)
c      dsdePl   (dp,ar(ncomp,ncomp),io)   material jacobian matrix
c      epsZZ    (dp,sc,o)                 strain epsZZ for plane stress,
c                                         define it when accounting for thickness change 
c                                         in shell and plane stress states
c
c*************************************************************************
c
c      ncomp   6   for 3D
c      ncomp   4   for plane strain, axisymmetric (nShear = 1)
c      ncomp   3   for plane stress (nShear = 1)
c      ncomp   3   for 3D beam (nShear = 2), beam188/189
c      ncomp   1   for 1D beam, link180
c
c      stresss and strains, plastic strain vectors
c          11, 22, 33, 12, 23, 13    for 3D
c          11, 22, 33, 12            for Plane strain and axisymmetry
c          11, 22, 12                for Plane stress
c          11, 13, 12                for 3d beam
c          11                        for 1D
c
c      material jacobian matrix
c        3D
c           dsdePl    |  1111   1122   1133   1112   1123   1113 |
c           dsdePl    |  2211   2222   2233   2212   2223   2213 |
c           dsdePl    |  3311   3322   3333   3312   3323   3313 |
c           dsdePl    |  1211   1222   1233   1212   1223   1213 |
c           dsdePl    |  2311   2322   2333   2312   2323   2313 |
c           dsdePl    |  1311   1322   1333   1312   1323   1313 |
c        plane strain, axisymmetric
c           dsdePl    |  1111   1122   1133   1112 |
c           dsdePl    |  2211   2222   2233   2212 |
c           dsdePl    |  3311   3322   3333   3312 |
c           dsdePl    |  1211   1222   1233   1212 |
c        plane stress
c           dsdePl    |  1111   1122   1112 |
c           dsdePl    |  2211   2222   2212 |
c           dsdePl    |  1211   1222   1212 |
c        3d beam plasticity
c           dsdePl    |  1111   1113   1112 |
c           dsdePl    |  1311   1313   1312 |
c           dsdePl    |  1211   1213   1212 |
c        1d
c           dsdePl    |  1111 |
c
c*************************************************************************
#include "impcom.inc"
c
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp), sigi(ncomp),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),       
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
c
c***************** User defined part *************************************
c
c --- parameters
c
      INTEGER          NEWTON, mcomp
      DOUBLE PRECISION HALF, ONE, TWO, SMALL, SQTWOTHIRD,
     &                 ZERO, TWOTHIRD, ONEDM02, ONEDM05, sqTiny
      PARAMETER       (ZERO       = 0.d0,
     &                 HALF       = 0.5d0,
     &                 ONE        = 1.d0,
     &                 TWO        = 2.d0,
     &                 SMALL      = 1.d-08,
     &                 sqTiny     = 1.d-20,
     &                 ONEDM02    = 1.d-02,
     &                 ONEDM05    = 1.d-05,
     &                 TWOTHIRD   = 2.0d0/3.0d0,
     &                 SQTWOTHIRD = 0.816496580927726030d0,
     &                 NEWTON     = 20,
     &                 mcomp      = 3
     &                 )
c
c --- local variables
c
c      sigElp   (dp,ar(3  ),l)            trial stress
c      dsdeEl   (dp,ar(3,3),l)            elastic moduli
c      pleq_t   (dp,sc     ,l)            equivalent plastic strain at beginnig of time increment
c      pleq     (dp,sc     ,l)            equivalent plastic strain at end of time increment
c      dpleq    (dp,sc     ,l)            incremental equivalent plastic strain
c      gamma    (dp,sc     ,l)            variable for solving incremental equivalent plastic strain
c      dgamma   (dp,sc     ,l)            correction of gamma
c      sigy_t   (dp,sc     ,l)            yield stress at beginnig of time increments
c      sigy     (dp,sc     ,l)            yield stress at end of time increment
c      young    (dp,sc     ,l)            Young's modulus
c      posn     (dp,sc     ,l)            Poiss's ratio
c      sigy0    (dp,sc     ,l)            initial yield stress
c      dsigdep  (dp,sc     ,l)            plastic slop
c      twoG     (dp,sc     ,l)            two time of shear moduli
c      funcf    (dp,sc     ,l)            nonlinear function to be solved for gamma
c      dFdep    (dp,sc     ,l)            derivative of nonlinear function over gamma
c
c --- temperary variables for solution purpose
c      i, j
c      c1, c2, c3, fratio
c      wk1(3), wk2(3), wk3(3), wk4(3) vector working arrays
c
      EXTERNAL         vmove, vzero, vapb1, vamb1,get_ElmData
      DOUBLE PRECISION sigElp(mcomp), dsdeEl(mcomp,mcomp), 
     &                 wk1(3), wk2(3), wk3(3), wk4(3)

      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7, var8

      INTEGER          i, j, k
      DOUBLE PRECISION pleq_t,  sigy_t , sigy,
     &                 cpleq, dpleq,   pleq,    twoG,    et,
     &                 young, posn,    sigy0,   dsigdep, 
     &                 gamma, dgamma,  dfdga,   dplga,   fratio,
     &                 funcFb,funcFb2, funcf,   dFdep,
     &                 c1, c2, c3, c4, c5
      DOUBLE PRECISION pv(3)
      data pv/TWOTHIRD, TWO, TWO/
c*************************************************************************
c
      keycut   = 0
c *** equivalent plastic strain at beginning of time step
      pleq_t   = ustatev(ncomp+1)
      pleq     = pleq_t
c *** get Young's modulus and Poisson's ratio, initial yield stress and slope of stress-strain
      young    = prop(1)
      posn     = prop(2)
      sigy0    = prop(3)
      et       = prop(4)
c *** calculate plastic slope
      dsigdep  = young * et/(young - et)
      twoG     = young / (ONE+posn)
c *** define tsstif(1) since it is used for calculation of hourglass stiffness
      tsstif(1) = HALF * twoG
c
c *** calculate elastic stiffness matrix
c
      call vzero(dsdeEl(1,1), ncomp * ncomp)
      c1 = twoG * HALF
      dsdeEl (1,1) = young
      dsdeEl (2,2) = c1
      dsdeEl (3,3) = c1
      DO i = 1, ncomp
         wk3(i) = dsdeEl(i,i)
      END DO
c *** calculate predicted strain 
      call vmove(Strain(1), wk1(1), ncomp)
      call vapb1(wk1(1), dStrain(1), ncomp)
      call vamb1(wk1(1), ustatev(1), ncomp)

c
c *** get initial stress
      call vzero(sigi(1),ncomp)
      i = ncomp
      call get_ElmData ('ISIG', elemId,kDomIntPt, i, sigi)

c
c *** calculate the trial stress and 
c     copy elastic moduli dsdeEl to material Jacobian matrix
      call vmove(dsdeEl(1,1), dsdePl(1,1), ncomp * ncomp)
      do i=1,ncomp
         sigElp(i) = wk3(i) * wk1(i) + sigi(i)
      end do
c
      funcFb2 = ZERO
      DO i = 1, ncomp
        funcFb2 = funcFb2 + pv(i) * sigElp(i) * sigElp(i)
      END DO
      funcFb = sqrt(funcFb2)

c *** compute current yield stress
      sigy    = sigy0 + dsigdep * pleq
c
      fratio = funcFb/sigy - SQTWOTHIRD
c *** check for yielding
      IF (fratio .LE. -SMALL) GO TO 500
      sigy_t  = sigy

      DO i = 1, ncomp
         wk3(i) = wk3(i) * pv(i)
      END DO

      gamma    = ZERO
      dplga    = ZERO
      dfdga    = ZERO
      dpleq    = ZERO
      pleq     = pleq_t 

c *** Local New-Raphson procedure for solving the gamma and 
c     thus the incremental equivalent plastic strain
      DO k=1,NEWTON
         funcFb2 = ZERO
         dfdga   = ZERO
         DO j = 1 , ncomp
            c1 = ONE + gamma * wk3(j)
            c1 = ONE / c1
            c2 = sigElp(j) * c1
            wk4(j) = c2
            funcFb2 = funcFb2 + pv(j) * c2 * c2
            c2 = c2 * c2 * c1 * wk3(j) * pv(j)
            dfdga   = dfdga - c2
         END DO
         funcFb   = sqrt(funcFb2)
c ***    derivative of funcFb w.r.t. gamma
         dfdga   = dfdga / funcFb

c ***    calculate the incremental equivalent plastic strain
         dpleq    = gamma * SQTWOTHIRD * funcFb
c ***    update the total equivalent plastic strain
         pleq     = pleq_t + dpleq
c ***    current yield stress
         sigy     = sigy0 + dsigdep * pleq
c ***    calculate the residual
         funcf    = funcFb - SQTWOTHIRD * sigy
c ***    derivative of incremental equivalent plastic strain w.r.t. gamma
         dplga    = SQTWOTHIRD * (gamma * dfdga + funcFb)
c ***    derivative of residual function w.r.t. gamma
         dFdep    = dfdga - SQTWOTHIRD * dsigdep * dplga
c ***    correction of gamma
         dgamma   = -funcf / dFdep 
         gamma    = gamma   + dgamma
c ***    check for negative gamma
         gamma    = max (gamma, sqTiny)
         fratio   = funcf/ sigy
c
c ***    Check for convergence of local New-Raphson iteration
         IF (((abs(fratio) .LT. ONEDM05 ) .AND.
     &        (abs(dgamma) .LT. ONEDM02*gamma)) .OR.
     &       ((abs(fratio) .LT. ONEDM05 ) .AND.
     &        (dgamma      .LE. sqTiny  ) .AND.
     &        ( gamma      .LE. sqTiny  )))  GO TO 100

      END DO
c
c *** Uncovergence, set keycut to 1 for bisection/cutback the load increment
      keycut   = 1
      GO TO 990
 100  CONTINUE
c
c *** update stresses
      call vmove(wk4(1), stress(1), ncomp)

c *** calculate incremental plastic strain 
      DO j = 1, ncomp
         wk2(j) = gamma * pv(j) * wk4(j)
      END DO
c *** update plastic strains
      call vapb1(epsPl(1),wk2(1),ncomp)

c *** Update state variables
      ustatev(ncomp+1) = pleq
      do i=1,ncomp
         ustatev(i) = epsPl(i)
      end do

c *** update plastic work
      sedPl     = sedPl + HALF * (sigy_t+sigy) * dpleq

      c1     = TWOTHIRD * dsigdep
      c3     = c1 * funcFb2 / (ONE - c1 * gamma)
      DO j = 1 , ncomp
         c1 = ONE / (ONE + gamma * wk3(j))
         wk3(j) = wk3(j) * c1 / pv(j)
      END DO
      DO j = 1 , ncomp
         wk4(j) = wk4(j) * pv(j)
      END DO
      DO j = 1 , ncomp
         c3 = c3 + wk4(j) * wk4(j) * wk3(j)
      END DO
      DO j = 1 , ncomp
         wk4(j) = wk4(j) * wk3(j)
      END DO

      c3     = ONE / c3
      DO i=1,ncomp
         dsdePl(i,i) = wk3(i)
      END DO

c *** Calculate the plastic Jacobian

      DO i=1,ncomp
         DO j=1,ncomp
            dsdePl(i,j) =    dsdePl(i,j) - c3 * wk4(i) * wk4(j)
         END DO
      END DO

      goto 600

  500 continue

c *** Update stress in case of elastic/unloading
      call vmove(sigElp(1),stress(1),ncomp)

  600 continue
c *** elastic strain energy
      sedEl = ZERO
      DO i = 1 , ncomp
         sedEl = sedEl + stress(i)*(Strain(i)+dStrain(i)-epsPl(i))
      END DO
      sedEl = sedEl * HALF
      ustatev(nStatev) = funcFb / SQTWOTHIRD

 990  CONTINUE
c
      return
      end
*deck,usermatps    USERDISTRIB  parallel                                gal
      subroutine usermatps(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ, 
     &                   var1, var2, var3, var4, var5, 
     &                   var6, var7, var8)
c*************************************************************************
c     *** primary function ***
c
c           user defined material constitutive model
c
c      Attention:
c           User must define material constitutive law properly
c           according to the stress state such as 3D, plane strain
c           and axisymmetry, plane stress and beam.
c
c           a 3D material constitutive model can use for
c           plane strain and axisymmetry cases.
c
c           When using shell elements, a plane stress algorithm
c           must be use.
c
c                                             gal July, 1999
c
c       The following demonstrates a USERMAT subroutine for
c       a plasticity model of plane stress state (such as PLANE182,
c       PLANE183 or SHELL181). The plasticity model is the same 
c       as ANSYS TB, BISO.
c
c       See "ANSYS user material subroutine USERMAT" for detailed
c       description of how to write a USERMAT routine.
c
c*************************************************************************
c
c     input arguments
c     ===============
c      matId     (int,sc,i)               material #
c      elemId    (int,sc,i)               element #
c      kDomIntPt (int,sc,i)               "k"th domain integration point
c      kLayer    (int,sc,i)               "k"th layer
c      kSectPt   (int,sc,i)               "k"th Section point
c      ldstep    (int,sc,i)               load step number
c      isubst    (int,sc,i)               substep number
c      nDirect   (int,sc,in)              # of direct components
c      nShear    (int,sc,in)              # of shear components
c      ncomp     (int,sc,in)              nDirect + nShear
c      nStatev   (int,sc,l)               Number of state variables
c      nProp     (int,sc,l)               Number of material ocnstants
c
c      Temp      (dp,sc,in)               temperature at beginning of
c                                         time increment
c      dTemp     (dp,sc,in)               temperature increment 
c      Time      (dp,sc,in)               time at beginning of increment (t)
c      dTime     (dp,sc,in)               current time increment (dt)
c
c      Strain   (dp,ar(ncomp),i)          Strain at beginning of time increment
c      dStrain  (dp,ar(ncomp),i)          Strain increment
c      prop     (dp,ar(nprop),i)          Material constants defined by TB,USER
c      coords   (dp,ar(3),i)              current coordinates
c      defGrad_t(dp,ar(3,3),i)            Deformation gradient at time t
c      defGrad  (dp,ar(3,3),i)            Deformation gradient at time t+dt
c
c     input output arguments              
c     ======================             
c      stress   (dp,ar(nTesn),io)         stress
c      ustatev   (dp,ar(nStatev),io)      user state variables
c            ustatev(1)                     - equivalent plastic strain
c            ustatev(2) - ustatev(1+ncomp)   - plastic strain vector
c            ustatev(nStatev)               - von-Mises stress
c      sedEl    (dp,sc,io)                elastic work
c      sedPl    (dp,sc,io)                plastic work
c      epseq    (dp,sc,io)                equivalent plastic strain
c      tsstif   (dp,ar(2),io)             transverse shear stiffness
c                                         tsstif(1) - Gxz
c                                         tsstif(2) - Gyz
c                                         tsstif(1) is also used to calculate hourglass
c                                         stiffness, this value must be defined when low
c                                         order element, such as 181, 182, 185 with uniform 
c                                         integration is used.
c      var?     (dp,sc,io)                not used, they are reserved arguments 
c                                         for further development
c
c     output arguments
c     ================
c      keycut   (int,sc,io)               loading bisect/cut control
c                                         0 - no bisect/cut
c                                         1 - bisect/cut 
c                                         (factor will be determined by ANSYS solution control)
c      dsdePl   (dp,ar(ncomp,ncomp),io)   material jacobian matrix
c      epsZZ    (dp,sc,o)                 strain epsZZ for plane stress, 
c                                         define it when accounting for thickness change 
c                                         in shell and plane stress states
c
c*************************************************************************
c
c      ncomp   6   for 3D  (nshear=3)
c      ncomp   4   for plane strain or axisymmetric (nShear = 1)
c      ncomp   3   for plane stress (nShear = 1)
c      ncomp   3   for 3d beam      (nShear = 2)
c      ncomp   1   for 1D (nShear = 0)
c
c      stresss and strains, plastic strain vectors
c          11, 22, 33, 12, 23, 13    for 3D
c          11, 22, 33, 12            for plane strain or axisymmetry
c          11, 22, 12                for plane stress
c          11, 13, 12                for 3d beam
c          11                        for 1D
c
c      material jacobian matrix
c        3D
c           dsdePl    |  1111   1122   1133   1112   1123   1113 |
c           dsdePl    |  2211   2222   2233   2212   2223   2213 |
c           dsdePl    |  3311   3322   3333   3312   3323   3313 |
c           dsdePl    |  1211   1222   1233   1212   1223   1213 |
c           dsdePl    |  2311   2322   2333   2312   2323   2313 |
c           dsdePl    |  1311   1322   1333   1312   1323   1313 |
c        plane strain or axisymmetric (11, 22, 33, 12)
c           dsdePl    |  1111   1122   1133   1112 |
c           dsdePl    |  2211   2222   2233   2212 |
c           dsdePl    |  3311   3322   3333   3312 |
c           dsdePl    |  1211   1222   1233   1212 |
c        plane stress (11, 22, 12)
c           dsdePl    |  1111   1122   1112 |
c           dsdePl    |  2211   2222   2212 |
c           dsdePl    |  1211   1222   1212 |
c        3d beam (11, 13, 12)
c           dsdePl    |  1111   1113   1112 |
c           dsdePl    |  1311   1313   1312 |
c           dsdePl    |  1211   1213   1212 |
c        1d
c           dsdePl    |  1111 |
c
c*************************************************************************
#include "impcom.inc"
c
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp), sigi(ncomp),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),       
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
c
c***************** User defined part *************************************
c
c --- parameters
c
      INTEGER          NEWTON, mcomp
      DOUBLE PRECISION HALF, THIRD, ONE, TWO, SMALL, 
     &                 SQTWOTHIRD, SQTWO1,
     &                 ZERO, TWOTHIRD, ONEDM02, ONEDM05, sqTiny
      PARAMETER       (ZERO       = 0.d0,
     &                 HALF       = 0.5d0,
     &                 THIRD      = 1.d0/3.d0,
     &                 ONE        = 1.d0,
     &                 TWO        = 2.d0,
     &                 SMALL      = 1.d-08,
     &                 sqTiny     = 1.d-20,
     &                 ONEDM02    = 1.d-02,
     &                 ONEDM05    = 1.d-05,
     &                 TWOTHIRD   = 2.0d0/3.0d0,
     &                 SQTWOTHIRD = 0.816496580927726030d0,
     &                 SQTWO1     = 0.707106769084930420d0,
     &                 NEWTON     = 20,
     &                 mcomp      = 6
     &                 )
c
c --- local variables
c
c      sigElp   (dp,ar(6  ),l)            trial stress
c      dsdeEl   (dp,ar(6,6),l)            elastic moduli
c      pleq_t   (dp,sc     ,l)            equivalent plastic strain at beginnig of time increment
c      pleq     (dp,sc     ,l)            equivalent plastic strain at end of time increment
c      dpleq    (dp,sc     ,l)            incremental equivalent plastic strain
c      gamma    (dp,sc     ,l)            variable for solving incremental equivalent plastic strain
c      dgamma   (dp,sc     ,l)            correction of gamma
c      sigy_t   (dp,sc     ,l)            yield stress at beginnig of time increments
c      sigy     (dp,sc     ,l)            yield stress at end of time increment
c      young    (dp,sc     ,l)            Young's modulus
c      posn     (dp,sc     ,l)            Poiss's ratio
c      sigy0    (dp,sc     ,l)            initial yield stress
c      dsigdep  (dp,sc     ,l)            plastic slop
c      twoG     (dp,sc     ,l)            two time of shear moduli
c      funcf    (dp,sc     ,l)            nonlinear function to be solved for dpleq
c      dFdep    (dp,sc     ,l)            derivative of nonlinear function over dpleq
c
c --- temperary variables for solution purpose
c      i, j
c      con1eOv2qEl, oneOv3G, qElOv3G, con1, con2, fratio
c
      EXTERNAL         vmove, vzero, vapb1, get_ElmData
      DOUBLE PRECISION sigElp(mcomp), dsdeEl(mcomp,mcomp), 
     &                 wk1(3), wk2(3), wk3(3), wk4(3)

      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7, var8

      INTEGER          i, j, k
      DOUBLE PRECISION pleq_t,  sigy_t , sigy,
     &                 dpleq,   pleq,    twoG,    et,
     &                 young, posn,    sigy0,   dsigdep, tEo1pm,
     &                 gamma, dgamma,  dfdga,   dplga, 
     &                 funcFb,funcFb2, funcf,   dFdep,   fratio,
     &                 con1,  con2,    con3,  con4,
     &                 con2p1, ocon2p1,
     &                 ocon2p2, con4p1, ocon4p1, ocon4p2,
     &                 c1, c2, c3,c4, c5,dperr(3)
c*************************************************************************
c
      keycut   = 0
c *** equivalent plastic strain at beginning of time step
      pleq_t   = ustatev(1)
      pleq     = pleq_t
c *** get Young's modulus and Poisson's ratio, initial yield stress and slope of stress-strain
      young    = prop(1)
      posn     = prop(2)
      sigy0    = prop(3)
      et       = prop(4)
c *** calculate plastic slope
      dsigdep  = young * et/(young - et)
      twoG     = young / (ONE+posn)
      tEo1pm   = THIRD * young /(ONE - posn)
c *** define tsstif(1) since it is used for calculation of hourglass stiffness
      tsstif(1) = HALF * twoG
c
c *** calculate elastic stiffness matrix (3d)
c
      c1 = ONE - posn * posn
      c2 = young / c1
      c3 = posn * c2
      dsdeEl (1,1) = c2
      dsdeEl (1,2) = c3
      dsdeEl (1,3) = ZERO
      dsdeEl (2,2) = c2
      dsdeEl (2,3) = ZERO
      dsdeEl (3,3) = HALF * twoG
      do i=1,ncomp-1
        do j=i+1,ncomp
          dsdeEl(j,i)=dsdeEl(i,j)
        end do
      end do
      call vmove(ustatev(2), epsPl(1), ncomp)
c *** calculate elastic strain
      do i=1,ncomp
         wk1(i) = Strain(i) - epsPl(i) + dStrain(i)
      end do

c
c *** get initial stress
      call vzero(sigi(1),ncomp)
      i = ncomp
      call get_ElmData ('ISIG', elemId,kDomIntPt, i, sigi)

c
c *** calculate the trial stress and 
c     copy elastic moduli dsdeEl to material Jacobian matrix
      do i=1,ncomp
         sigElp(i) = ZERO
         do j=1,ncomp
            dsdePl(j,i) = dsdeEl(j,i)
            sigElp(i) = sigElp(i) + dsdeEl(j,i) * wk1(j)
         end do
         sigElp(i) = sigElp(i) + sigi(i)
      end do
c
      wk1(1)   = SQTWO1 * ( sigElp(1) +  sigElp(2) )
      wk1(2)   = SQTWO1 * (-sigElp(1) +  sigElp(2) )
      wk1(3)   = sigElp(3)
c
      funcFb2 =  THIRD * wk1(1)  * wk1(1) +
     &   wk1(2)  * wk1(2) + TWO *  wk1(3)  * wk1(3)

      funcFb  = sqrt(funcFb2)

c *** compute current yield stress
      sigy    = sigy0 + dsigdep * pleq
c
      fratio = funcFb/sigy - SQTWOTHIRD
c *** check for yielding
      IF (fratio .LE. -SMALL) GO TO 500
      sigy_t  = sigy

      gamma    = ZERO
      dplga    = ZERO
      dfdga    = ZERO
      con1     = THIRD * wk1(1) * wk1(1)
      con2     = wk1(2) * wk1(2) + TWO * wk1(3) * wk1(3)
      con3     = - con1 * tEo1pm
      con4     = - con2 * twoG
      funcf    = funcFb - SQTWOTHIRD * sigy
      dfdga    =  (con3 + con4) / funcFb
      dFdep    = dfdga - TWOTHIRD * dsigdep * funcFb
      dgamma   = -funcf / dFdep
      gamma    = gamma   + dgamma
      gamma    = max (gamma, sqTiny)

      DO k=1,NEWTON

         con2p1   = ONE + tEo1pm * gamma
         con4p1   = ONE + twoG * gamma
         ocon2p1  = ONE / con2p1
         ocon2p2  = ocon2p1 * ocon2p1
         ocon4p1  = ONE / con4p1
         ocon4p2  = ocon4p1 * ocon4p1
         funcFb2  =  con1 * ocon2p2 + con2 * ocon4p2
         funcFb   = sqrt(funcFb2)
c
c ***    calculate the incremental equivalent plastic strain
         dpleq    = gamma * SQTWOTHIRD * funcFb
c ***    update the total equivalent plastic strain
         pleq     = pleq_t + dpleq
c ***    current yield stress
         sigy     = sigy0 + dsigdep * pleq

c ***    calculate the residual
         funcf    = funcFb - SQTWOTHIRD * sigy
c ***    derivative of funcFb w.r.t. gamma
         dfdga    =  ( con3 * ocon2p2 * ocon2p1
     &               + con4 * ocon4p2 * ocon4p1) / funcFb
c ***    derivative of incremental equivalent plastic strain w.r.t. gamma
         dplga    = SQTWOTHIRD * (gamma * dfdga + funcFb)
c ***    derivative of residual function w.r.t. gamma
         dFdep  = dfdga - SQTWOTHIRD * dsigdep * dplga
         fratio   = funcf/ sigy

         dgamma   = -funcf / dFdep
         gamma    = gamma   + dgamma
         gamma    = max (gamma, sqTiny)
c
         IF (((abs(fratio) .LT. ONEDM05 ) .AND.
     &        (abs(dgamma) .LT. ONEDM02*gamma)) .OR.
     &       ((abs(fratio) .LT. ONEDM05 ) .AND.
     &        (dgamma      .LE. sqTiny  ) .AND.
     &        ( gamma      .LE. sqTiny  )))  GO TO 100

      END DO
c
c *** Uncovergence, set keycut to 1 for bisect/cut
      keycut   = 1
      GO TO 990
 100  CONTINUE
c
c *** update stresses
      c1 = ONE / ( ONE + gamma * tEo1pm)
      c2 = ONE / ( ONE + gamma * twoG  )
      stress(1) = SQTWO1 * (c1 * wk1(1) - c2 * wk1(2))
      stress(2) = SQTWO1 * (c1 * wk1(1) + c2 * wk1(2))
      stress(3) = c2 * wk1(3)
c
      wk1(1)    = SQTWO1 * ( stress(1) +  stress(2))
      wk1(2)    = SQTWO1 * (-stress(1) +  stress(2))
      wk1(3)    = stress(3)
      con1  = SQTWO1 * gamma
c *** calculate incremental plastic strain, wk2(i)
      wk2(1)    = con1 * ( THIRD * wk1(1) -  wk1(2))
      wk2(2)    = con1 * ( THIRD * wk1(1) +  wk1(2))
      wk2(3)    = TWO * gamma * wk1 (3)

c ***  update plastic strains
      call vapb1(epsPl(1),wk2(1),ncomp)
      epseq     = pleq
      ustatev(1) = pleq
c *** calculate plastic work
      sedPl     = sedPl + HALF * (sigy_t+sigy) * dpleq

c *** consistent tangent stiffness matrix
      con3      = TWOTHIRD*dsigdep*funcFb2/(ONE-TWOTHIRD*dsigdep*gamma)
      call vmove(wk1(1),wk4(1),3)
      c1 = wk4(1)
      c2 = wk4(2)
      c3 = wk4(3)
      c4 = tEo1pm / (ONE + tEo1pm * gamma)
      c5 = twoG / (ONE + twoG * gamma)
      wk4(1) = SQTWO1 * ( c4 * c1 - c5 * c2)
      wk4(2) = SQTWO1 * ( c4 * c1 + c5 * c2)
      wk4(3) = c5 * c3

      call vmove(wk4(1),wk3(1),3)
      c1 = wk3(1)
      c2 = wk3(2)
      wk3(1) = SQTWO1 * ( c1+c2)
      wk3(2) = SQTWO1 * (-c1+c2)

      con3     = con3 + THIRD*wk1(1)*wk3(1)+
     &           wk1(2)*wk3(2)+TWO*wk1(3)*wk3(3)
      c1  = 3.0d0 * tEo1pm / (ONE + tEo1pm * gamma)
      c2  = twoG / (ONE + twoG * gamma)
      call vzero(dsdePl(1,1),9)
      dsdePl(1,1) = HALF * (c1 + c2)
      dsdePl(2,1) = HALF * (c1 - c2)
      dsdePl(1,2) = dsdePl(2,1)
      dsdePl(2,2) = dsdePl(1,1)
      dsdePl(3,3) = HALF * c2
      con3     = ONE / con3
      DO i=1,3
         DO j=1,3
            dsdePl(i,j) = dsdePl(i,j) - con3 * wk4(i) * wk4(j)
         END DO
      END DO

      goto 600

  500 continue

c *** Update stress in case of elastic/unloading
      call vmove(sigElp(1),stress(1),ncomp)

  600 continue
c *** elastic strain and put into epsZZ
      epsZZ    = -posn/young * (stress(1) + stress(2))
c *** add plastic strain to total strain epsZZ
      epsZZ    = epsZZ - (epsPl(1) + epsPl(2))

c *** elastic strain energy
      sedEl = ZERO
      DO i = 1 , ncomp
         sedEl = sedEl + stress(i)*(Strain(i)+dStrain(i)-epsPl(i))
      END DO
      sedEl    = sedEl * HALF
      ustatev(nStatev) = sigy
      do i=1,ncomp
         ustatev(i+1) = epsPl(i)
      end do

 990  CONTINUE
c
      return
      end
