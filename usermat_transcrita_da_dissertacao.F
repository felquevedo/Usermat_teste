      subroutine usermat3d(
     & matId, elemId,kDomIntPt, kLayer, kSectPt,
     & ldstep,isubst,keycut,
     & nDirect,nShear,ncomp,nStatev,nProp,
     & Time,dTime,Temp,dTemp,
     & stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     & Strain,dStrain, epsPl, prop, coords,
     & var0, defGrad_t, defGrad,
     & tsstif, epsZZ, cutFactor,
     & var1, var2, var3, var4, var5,
     & var6, var7)
c
c input arguments
c ===============
c matId (int,sc,i) material #
c elemId (int,sc,i) element #
c kDomIntPt (int,sc,i) "k"th domain integration point
c kLayer (int,sc,i) "k"th layer
c kSectPt (int,sc,i) "k"th Section point
c ldstep (int,sc,i) load step number
c isubst (int,sc,i) substep number
c nDirect (int,sc,in) # of direct components
c nShear (int,sc,in) # of shear components
c ncomp (int,sc,in) nDirect + nShear
c nstatev (int,sc,l) Number of state variables
c nProp (int,sc,l) Number of material ocnstants
c
c Temp (dp,sc,in) temperature at beginning of
c time increment
c dTemp (dp,sc,in) temperature increment
c Time (dp,sc,in) time at beginning of increment (t)
c dTime (dp,sc,in) current time increment (dt)
c
c Strain (dp,ar(ncomp),i) Strain at beginning of time increment
c dStrain (dp,ar(ncomp),i) Strain increment
c prop (dp,ar(nprop),i) Material constants defined by TB,USER
c coords (dp,ar(3),i) current coordinates
c defGrad_t(dp,ar(3,3),i) Deformation gradient at time t
c defGrad (dp,ar(3,3),i) Deformation gradient at time t+dt
c
c input output arguments
c ======================
c stress (dp,ar(nTesn),io) stress
c ustatev (dp,ar(nstatev),io) user state variable
c ustatev(1) - equivalent plastic strain
c ustatev(2) - statev(1+ncomp) - plastic strain vector
c ustatev(nStatev) - von-Mises stress
c sedEl (dp,sc,io) elastic work
c sedPl (dp,sc,io) plastic work
c epseq (dp,sc,io) equivalent plastic strain
c tsstif (dp,ar(2),io) transverse shear stiffness
c tsstif(1) - Gxz
c tsstif(2) - Gyz
c tsstif(1) is also used to calculate
c hourglass
c stiffness, this value must be defined when
c low
c order element, such as 181, 182, 185 with
c uniform
c integration is used.
c var? (dp,sc,io) not used, they are reserved arguments
c for further development
c
c
c
c output arguments
c ================
c keycut (int,sc,io) loading bisect/cut control
c 0 - no bisect/cut
c 1 - bisect/cut
c (factor will be determined by ANSYS solution
c control)
c dsdePl (dp,ar(ncomp,ncomp),io) material jacobian matrix
c epsZZ (dp,sc,o) strain epsZZ for plane stress,
c define it when accounting for thickness
c change
c in shell and plane stress states
c cutFactor(dp,sc,o) time step size cut-back factor
c define it if a smaller step size is wished
c recommended value is 0~1
c
c*************************************************************************
c
c ncomp 6 for 3D (nshear=3)
c ncomp 4 for plane strain or axisymmetric (nShear = 1)
c ncomp 3 for plane stress (nShear = 1)
c ncomp 3 for 3d beam (nShear = 2)
c ncomp 1 for 1D (nShear = 0)
c
c stresss and strains, plastic strain vectors
c 11, 22, 33, 12, 23, 13 for 3D
c 11, 22, 33, 12 for plane strain or axisymmetry
c
c 3D
c dsdePl | 1111 1122 1133 1112 1123 1113 |
c dsdePl | 2211 2222 2233 2212 2223 2213 |
c dsdePl | 3311 3322 3333 3312 3323 3313 |
c dsdePl | 1211 1222 1233 1212 1223 1213 |
c dsdePl | 2311 2322 2333 2312 2323 2313 |
c dsdePl | 1311 1322 1333 1312 1323 1313 |
c plane strain or axisymmetric (11, 22, 33, 12)
c dsdePl | 1111 1122 1133 1112 |
c dsdePl | 2211 2222 2233 2212 |
c dsdePl | 3311 3322 3333 3312 |
c dsdePl | 1211 1222 1233 1212 |
c
c*************************************************************************
#include "impcom.inc"
c
      INTEGER
     & matId, elemId,
     & kDomIntPt, kLayer, kSectPt,
     & ldstep,isubst,keycut,
     & nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION
     & Time, dTime, Temp, dTemp,
     & sedEl, sedPl, epseq, epsZZ, cutFactor
      DOUBLE PRECISION
     & stress (ncomp ), ustatev (nStatev),
     & dsdePl (ncomp,ncomp),
     & Strain (ncomp ), dStrain (ncomp ),
     & epsPl (ncomp ), prop (nProp ),
     & coords (3),
     & defGrad (3,3), defGrad_t(3,3),
     & tsstif (2)
c
c***************** User defined part *************************************
c
c --- parameters
c
      INTEGER mcomp
      DOUBLE PRECISION HALF, THIRD, ONE, TWO, SMALL, ONEHALF,
     & ZERO, TWOTHIRD, ONEDM02, ONEDM05, sqTiny
      PARAMETER (ZERO = 0.d0,
     & HALF = 0.5d0,
     & THIRD = 1.d0/3.d0,
     & ONE = 1.d0,
     & TWO = 2.d0,
     & SMALL = 1.d-08,
     & sqTiny = 1.d-20,
     & ONEDM02 = 1.d-02,
     & ONEDM05 = 1.d-05,
     & ONEHALF = 1.5d0,
     & TWOTHIRD = 2.0d0/3.0d0,
     & mcomp = 6
     & )
c
c --- local variables
c
c sigElp (dp,ar(6 ),l) trial stress
c dsdeEl (dp,ar(6,6),l) elastic moduli
c sigDev (dp,ar(6 ),l) deviatoric stress tensor
c dfds (dp,ar(6 ),l) derivative of the yield function
c JM (dp,ar(6,6),l) 2D matrix for a 4 order tensor
c pEl (dp,sc ,l) hydrostatic pressure stress
c qEl (dp,sc ,l) von-Mises stress
c pleq_t (dp,sc ,l) equivalent plastic strain at beginnig of time increment
c pleq (dp,sc ,l) equivalent plastic strain at end of time increment
c dpleq (dp,sc ,l) incremental equivalent plastic strain
c sigy_t (dp,sc ,l) yield stress at beginnig of time increments
c sigy (dp,sc ,l) yield stress at end of time increment
c young (dp,sc ,l) Young's modulus
c posn (dp,sc ,l) Poiss's ratio
c sigy0 (dp,sc ,l) initial yield stress
c dsigdep (dp,sc ,l) plastic slop
c twoG (dp,sc ,l) two time of shear moduli
c threeG (dp,sc ,l) three time of shear moduli
c
c --- temperary variables for solution purpose
c i, j
c threeOv2qEl, oneOv3G, qElOv3G, con1, con2, fratio
c
      EXTERNAL vzero, vmove, get_ElmData
      DOUBLE PRECISION sigElp(mcomp), dsdeEl(mcomp,mcomp), G(mcomp),
     & sigDev(mcomp), JM (mcomp,mcomp), dfds(mcomp),
     & sigi (mcomp), strainEl(mcomp)
      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     & var6, var7
      DATA G/1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,0.0D0/
      INTEGER i, j, k
      DOUBLE PRECISION pEl, qEl, pleq_t, sigy_t , sigy,
     & dpleq, pleq,
     & young, posn, sigy0, dsigdep,
     & elast1,elast2,
     & twoG, threeG, oneOv3G, qElOv3G, threeOv2qEl,
     & fratio, con1, con2, dperr(3)
c
c --- variaveis locais do usuario necessárias ao cálculo da fluência e retração do
c concreto
c
      ! Variaveis de entrada do usuario TBDATA: relacionadas aos dados do concreto
      !
      ! fck [kN/cm²] resistencia caracteristica do concreto
      ! s [adm] coef. que depende do tipo de concreto
      ! ti [dias] tempo da analise apartir da qual o concreto iniciou
      ! rh [%] umidade relativa do ambiente
      ! hf [cm] espessura ficticia
      ! ts [dias] idade do concreto no inicio da secagem
      ! t0 [dias] idade do concreto no inicio do carregamento
      ! tinf [dias] tempo final da análise
      ! betasc [adm] coeficiente que depende do tipo de cimento
      ! alpha [adm] efeito do tipo de cimento durante a cura
      ! ndec [pontos por decada] No. de pontos para discretização
      ! ligafluencia 0 - desliga a fluência, 1 - liga fluência
      ! ligaretracao 0 - desliga a retração, 1 - liga a retração
      ! ktime [dias] identificacao do dia na saida de dados
      ! kelemid [un] identificacao do elemento na saida de dados
      ! kkDomIntPt [un] identificacao do ponto de integracao na saida de dados
      !
c
c
c
      ! Variaveis calculadas no interior da USERMAT
      !
      ! tmat [dias] idade do material durante a análise (time-ti)
      ! fcm28 [kN/cm²] resistência média do concreto aos 28 dias
      ! sigmalim_CEB90 [kN/cm²] limite de tensão para validade da viscoelasticidade do concreto
      ! Ec28 [kN/cm²] módulo de elasticidade do concreto aos 28 dias
      ! mkelvin [un] no. máximo de elementos de Kelvin
      ! ndt [un] No. de pontos para calcular os parâmetros da cadeia de Kelvin
      ! nkelvin [un] No. de camadas da cadeia de Kelvin (DIAS,2013: N da eq3.3, p44)
      ! tau [cm²/(kN*dia)] tempos de retardação dos elementos da cadeia de Kelvin
      ! Emu [kN/cm²] módulos de elasticidade dos elementos da cadeia de Kelvin
      ! tmeio [dias] idade do material somada meio incremento de tempo (tmat+dtime/2)
      ! tmatadj [dias] idade do material durante a análise corrigida com a temperatura
      ! fi0 [adm] coeficiente de fluência que depende do tempo t0
      ! V [adm] volume solidificado
      ! dsdeElinv [cm²/kN] matriz de rigidez inversa
      ! epsloncf1 [adm] primeira parcela deform. fluência de uma unidade de Kelvin
      ! epsloncf2 [adm] segunda parcela deform. por fluência de uma unidade de Kelvin
      ! epsloncf3 [adm] deformação por fluência da unidade i
      ! deltay [adm] parametro do modelo do Bazant & Prasannan
      ! lbd [adm] parametro do modelo do Bazant & Prasannan
      ! tensor_stress(3,3) [kN/cm²] tensor das tensoes
      ! stressp(3) [kN/cm²] vetor com as tensoes principais
      ! ev(3,3) [kN/cm²] matriz com as direcoes principais
      ! nrot [un] numero de rotacoes do metodo que encontra as tensoes principais
      ! saida 0 - nao faz saida de dados, 1 - faz saida de dados
      !
      ! Outras variaveis
      !
      ! iott, wrinqr variaveis que guardam o numero da output window do Ansys
      !
      ! Variaveis de estado
      !
      ! ustatev(1) primeirarodada ! 9999 – indica primeira rodada
      ! ustatev(2) fcm28
      ! ustatev(3) Ec28
      ! ustatev(4) nkelvin
      ! ustatev(5) t0
      ! ustatev(6) tmat
      ! ustatev(7) stressp(1)
      ! ustatev(8) stressp(2)
      ! ustatev(9) stressp(3)
      ! ustatev(10) FV = 0.4*fcm28/stressp(3)
      ! ustatev(11..) Emu
      ! ustatev(21..) tau
      ! ustatev(31..) epsloncsant
      ! ustatev(41..) epsloncs
      ! ustatev(51..) depsloncs
      ! ustatev(101..) epsloncfantmu
      ! (111..) epsloncfantmu
      ! (121..) epsloncfantmu
      ! (130..) epsloncfantmu
      ! (140..) epsloncfantmu
      ! (150..) epsloncfantmu
      ! (160..) epsloncfantmu
      ! ustatev(201..) epsloncfmu
      ! (211..) epsloncfmu
      ! (221..) epsloncfmu
      ! (230..) epsloncfmu
      ! (240..) epsloncfmu
      ! (250..) epsloncfmu
      ! (260..) epsloncfmu
      ! ustatev(301..) depsloncf
      ! ustatev(401..) dstress
      !
      DOUBLE PRECISION fck,s,ti,rh,hf,ts,t0,tinf,betasc,alpha
      INTEGER ndec,ligafluencia,ligaretracao
      INTEGER kelemId,kkDomIntPt
      DOUBLE PRECISION ktime
      DOUBLE PRECISION tmat,fcm28,Ec28
      INTEGER mkelvin,ndt,nkelvin
      PARAMETER (mkelvin=7)
      DOUBLE PRECISION tau(mkelvin),EMU(mkelvin),tmeio,tmatadj,fi0,V,
     & dsdeElinv(mcomp,mcomp),epsloncf1(mcomp),
     & epsloncf2(mcomp),epsloncf3(mcomp),
     & deltay(mkelvin),lbd(mkelvin)
      EXTERNAL wrinqr
      INTEGER iott,wrinqr
      DOUBLE PRECISION wr_output
      DOUBLE PRECISION tensor_stress(3,3)
      DOUBLE PRECISION stressp(3)
      DOUBLE PRECISION ev(3,3)
      INTEGER nrot
      INTEGER saida
      DOUBLE PRECISION CEB90_fcm
      DOUBLE PRECISION CEB90_Ec
      DOUBLE PRECISION CEB90_epsloncs
      DOUBLE PRECISION CEB90_tadj
      DOUBLE PRECISION CEB90_tT
      DOUBLE PRECISION CEB90_fi0
c
c*************************************************************************
c
      !
      !*************************************************************************!
      ! PASSO1) Entrada de dados do script !
      !*************************************************************************!
      ligafluencia = prop(1)
      ligaretracao = prop(2)
      ti = prop(3)
      ts = prop(4)
      t0 = prop(5)
      tinf = prop(6)
      fck = prop(7)
      posn = prop(8)
      s = prop(9)
      rh = prop(10)
      hf = prop(11)
      betasc = prop(12)
      alpha = prop(13)
      ndec = prop(14)
      ktime = prop(15)
      kelemid = prop(16)
      kkDomIntPt = prop(17)
      iott = wrinqr(WR_OUTPUT)
      !
      !*************************************************************************!
      ! PASSO2) Preparando arquivo de saida com os valores do vetor ustatev !
      !*************************************************************************!
      saida = 0
      if (kelemId.EQ.elemId.
     & AND.kkDomIntPt.EQ.kDomIntPt.
     & AND.ktime.EQ.time) THEN
      saida = 1
      endif
      !
      if (saida.EQ.1) then
          open(1,file='ustatev.txt',access='sequential',
     & form="formatted",status="unknown",
     & position="append")
          write(1,*)'----------------------------------------------'
          write(1,*)'SAIDA DE DADOS DO VETOR USTATEV'
          write(1,*)'----------------------------------------------'
          write(1,*)'ID: time,elemId,kDomIntPt:'
          write(1,*) time,elemId,kDomIntPt
          write(1,*)'ldstep,isubst'
          write(1,*) ldstep,isubst
      endif
      !
      !*************************************************************************!
      ! PASSO3) Obtendo parametros do modelo de Kelvin !
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
                  call KELPAR(t0,tinf,ndec,nkelvin,ndt)
                  !
                  ! Calculando os tempos de retardação e módulos do ajuste
                  call KELVIN_betactt0(nkelvin,t0,tinf,ndec,ndt,
     & tau,Emu,temp,rh,hf)
                  !
                  ! Guardando os parâmetros
                  ustatev(2) = fcm28
                  ustatev(3) = Ec28
                  ustatev(4) = nkelvin
                  ustatev(5) = t0
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
          t0 = ustatev(5)
          do i=1,nkelvin
              Emu(i) = ustatev(10+i)
              tau(i) = ustatev(20+i)
          enddo
      endif
      !
      !*************************************************************************!
      ! PASSO4) Calculo da idade do material !
      !*************************************************************************!
      !
      if (time.GT.ti) then
          tmat = time-ti
      else
          tmat = 0.0d0
      endif
      ustatev(6) = tmat
      !
      !*************************************************************************!
      ! PASSO5) Calculo do incremento de deformacao por retracao !
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
      ! PASSO6) Calculo do incremento de deformacao por fluencia !
      !*************************************************************************!
      !
      ! Calcula apenas quando a idade da peça atinge a idade em que entra em carga
      if (tmat.GT.t0) then
          !
              ! Guardando deformação por fluência do ultimo passo que convergiu
                  do i=1,nkelvin
                      do j=1,ncomp
                          ustatev(100+(i-1)*10+j) = ustatev(200+(i-1)*10+j)
                      enddo
                  enddo
              !
          ! Calculando o fator do envelhecimento (com a idade da peça em carga (tmat-t0))
          tmatadj = CEB90_tadj(tmat-t0,alpha,temp,tmat)
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
      ! PASSO7) Calculo do incremento de deformacao elastica !
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
      tmatadj = CEB90_tadj(tmat-t0,alpha,temp,tmat)
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
      ! Descontando o incremento de deformação por fluência e retração do incremento de deformação
      c total (Equação 5.60)
      do i=1,ncomp
          dstrain(i) = dstrain(i)
     & -ustatev(300+i)*ligafluencia
     & -ustatev(50+i)*ligaretracao
      end do
      !
      !*************************************************************************!
      ! PASSO8) Calculo do incremento de tensão !
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
      ! PASSO9) Atualizando a tensão !
      !*************************************************************************!
      !
      do i=1,ncomp
          stress(i) = stress(i)+sigElp(i)
      end do
      !
      !*************************************************************************!
      ! PASSO10) Atualizando a matriz Jacobiana !
      !*************************************************************************!
      !
      do i=1,ncomp
          do j=1,ncomp
              dsdePl(j,i) = dsdeEl(j,i)
          end do
      end do
      !
      !*************************************************************************!
      ! PASSO11) Guardando a tensão incremento de tensao e deformação retração !
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
      !*************************************************************************!
      ! PASSO12) Calculando o fator de validação !
      !*************************************************************************!
      !
      ! Encontrando as tensões e direções principais
      call jacobi(stress,3,3,stressp,ev,nrot)
      !
      ! Ordenando as tensões e direções principais
      call eigsrt(stressp,ev,3,3)
      !
      ustatev(7) = stressp(1)
      ustatev(8) = stressp(2)
      ustatev(9) = stressp(3)
      ustatev(10) = abs(0.4d0*ustatev(2)/ustatev(9))
      !
      !*************************************************************************!
      ! PASSO13) Escrevendo o vetor ustatev na saida de dados !
      !*************************************************************************!
      if (saida.EQ.1) then
          do i=1,nStatev
              write(1,*)i, ustatev(i)
          enddo
          close(1,STATUS='KEEP')
      endif
      !
      return
      End
      subroutine matrizD(E,Poisson,D)
      !**************************************************************************!
      !** Função: matrizD **!
      !** **!
      !** Objetivo: calcula a matriz consitutiva do material isotrópico **!
      !** **!
      !** Situação: (28-09-2016) OK **!
      !** **!
      !**************************************************************************!
      DOUBLE PRECISION E ! módulo de elasticidade
      DOUBLE PRECISION Poisson ! coeficiente de Poisson
      DOUBLE PRECISION D(6,6) ! matriz constitutiva elástica isotrópica
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
      !** Função: CEB90_tT **!
      !** **!
      !** Objetivo: considera o efeito da temperatura na maturidade do **!
      !** concreto **!
      !** **!
      !** Referencia: CEB-MC90 **!
      !** Situação: (28-09-2016) OK **!
      !** **!
      !**************************************************************************!
      DOUBLE PRECISION temp ! [oC] temperatura
      DOUBLE PRECISION deltat ! [dias] intervalo em dias em que a temp. prevalece
      DOUBLE PRECISION CEB90_tT
      !
      CEB90_tT = deltat*dexp(13.65d0-4000.0d0/(273.0d0+temp))
      !
      end
      function CEB90_tadj(t,alpha,temp,deltat)
      !**************************************************************************!
      !** Função: CEB90_tadj **!
      !** **!
      !** Objetivo: considera o efeito do tipo de cimento e temperatura de **!
      !** cura na idade de carregamento. **!
      !** **!
      !** Referencia: CEB-MC90 **!
      !** Situação: (28-09-2016) OK **!
      !** **!
      !**************************************************************************!
      DOUBLE PRECISION t ! [dias] idade
      DOUBLE PRECISION alpha ! [adm] coeficiente que depende do tipo de cimento
      DOUBLE PRECISION temp ! [oC] temperatura
      DOUBLE PRECISION deltat ! [dias] intervalo em dias em que a temp. prevalece
      DOUBLE PRECISION tT ! [dias] idade ajustada
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
      !** Função: CEB90_fcm **!
      !** **!
      !** Objetivo: calcular a resistência média do concreto considerando o **!
      !** o tempo t e a temperatura constante temp. **!
      !** **!
      !** Referencia: CEB-MC90 **!
      !** Situação: (28-09-2016) OK **!
      !** **!
      !**************************************************************************!
      DOUBLE PRECISION fck ! [kN/cm²] resistência característica a compressão
      DOUBLE PRECISION t ! [dias] idade
      DOUBLE PRECISION s ! [adm] coeficiente que depende do tipo de cimento
      DOUBLE PRECISION temp ! [oC] temperatura
      DOUBLE PRECISION Bcc ! [adm] coeficiente
      DOUBLE PRECISION fcm ! [kN/cm²] resistência média a compressão
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
      !** Função: CEB90_Eci **!
      !** **!
      !** Objetivo: calcular o módulo de elasticidade do concreto considerando **!
      !** o tempo t e a temperatura constante temp. **!
      !** **!
      !** Referencia: CEB-MC90 **!
      !** Situação: (28-09-2016) OK **!
      !** **!
      !**************************************************************************!
      DOUBLE PRECISION fck ! [kN/cm²] resistência característica a compressão
      DOUBLE PRECISION s ! [adm] coeficiente que depende do tipo de cimento
      DOUBLE PRECISION t ! [dias] idade
      DOUBLE PRECISION temp ! [oC] temperatura
      DOUBLE PRECISION Ec0 ! [kN/cm²] módulo de elasticidade de referência
      DOUBLE PRECISION fcm0 ! [kN/cm²] resistência de referência
      DOUBLE PRECISION Eci ! [kN/cm²] módulo de elasticidade inicial
      DOUBLE PRECISION Bcc ! [adm] coeficiente
      DOUBLE PRECISION Be ! [adm] coeficiente
      DOUBLE PRECISION Ec ! [kN/cm²] módulo de elasticidade do concreto
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
      !** Função: CEB90_fi0 **!
      !** **!
      !** Objetivo: calcular o coeficiente de fluência que depende do tempo t0 **!
      !** **!
      !** Referencia: CEB-MC90 **!
      !** Situação: (28-09-2016) OK **!
      !** **!
      !**************************************************************************!
      DOUBLE PRECISION t0 ! [dias] idade no instante de aplicação da carga
      DOUBLE PRECISION rh ! [%] umidade relativa
      DOUBLE PRECISION hf ! [cm] espessura ficticia
      DOUBLE PRECISION fcm ! [kN/cm²] resistência média a compressão
      DOUBLE PRECISION temp ! [oC] temperatura
      DOUBLE PRECISION fcm0 ! [kN/cm²^] resistência de referencia
      DOUBLE PRECISION firh ! [adm] coeficiente
      DOUBLE PRECISION fit ! [adm] coeficiente
      DOUBLE PRECISION betafcm ! [adm] coeficiente
      DOUBLE PRECISION betat0 ! [adm] coeficiente
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
      !** Função: CEB90_betactt0 **!
      !** **!
      !** Objetivo: calcular o coeficiente de fluência que depende do **!
      !** intervalo de tempo de aplicação da carga t-t0 **!
      !** **!
      !** Referencia: CEB-MC90 **!
      !** Situação: (28-09-2016) OK **!
      !** **!
      !**************************************************************************!
      DOUBLE PRECISION t ! [dias] idade
      DOUBLE PRECISION t0 ! [dias] idade na aplicação da carga
      DOUBLE PRECISION rh ! [%] umidade relativa
      DOUBLE PRECISION hf ! [cm] espessura fictícia
      DOUBLE PRECISION temp ! [oC] temperatura
      DOUBLE PRECISION betah ! [adm] coeficiente
      DOUBLE PRECISION betat ! [adm] coeficiente
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
      !** Função: CEB90_epsloncs **!
      !** **!
      !** Objetivo: calcular o deformação de retração **!
      !** **!
      !** Referencia: CEB-MC90 **!
      !** Situação: (28-09-2016) OK **!
      !** **!
      !**************************************************************************!
      DOUBLE PRECISION t ! [dias] idade
      DOUBLE PRECISION ts ! [dias] idade no final da cura
      DOUBLE PRECISION fcm ! [kN/cm²] resistência média a compressão
      DOUBLE PRECISION betasc ! [adm] coeficien que depende do tipo de cimento
      DOUBLE PRECISION rh ! [%] umidade relativa
      DOUBLE PRECISION hf ! [cm] espessura ficticia
      DOUBLE PRECISION temp ! [oC] temperatura
      DOUBLE PRECISION fcm0 ! [kN/cm²] resistência de referencia
      DOUBLE PRECISION betasrh ! [adm] coeficiente
      DOUBLE PRECISION betarh ! [adm] coeficiente
      DOUBLE PRECISION betast ! [adm] coeficiente
      DOUBLE PRECISION betastts ! [adm] coeficiente
      DOUBLE PRECISION epsloncsfcm ! [adm] deformação específica
      DOUBLE PRECISION epsloncs0 ! [adm] deformação específica
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
      subroutine KELPAR(t0,tinf,ndec,nkelvin,ndt)
      !**************************************************************************!
      !** Subrotina: KELPAR                                                    **!
      !** **!
      !** Objetivo: subrotina para determinar o numero de camadas de kelvin **!
      !** e o numero de intervalos de tempo **!
      !** **!
      !** Referencia: DIAS,M.M. Análise numérica de vigas mistas aço-concreto **!
      !** pelo método dos elementos finitos: efeitos de longa **!
      !** duração. 2013. Dissertação (Mestrado em Engenharia) **!
      !** Programa de Pós Graduação em Engenharia Civil, UFRGS, **!
      !** Porto Alegre. **!
      !** **!
      !** BAZANT,Z.P.; PRASANNAN, S. Solidification Theory for **!
      !** concrete creep. I: Formulation, Journal of Enginnering **!
      !** mechanics, v.115, n.8, p.1691-1703, 1989. **!
      !** **!
      !** BAZANT,Z.P.; PRASANNAN, S. Solidification Theory for **!
      !** concrete creep. II: verification and application. **!
      !** Journal of Enginnering Mechanics, v.115, n.8, **!
      !** p.1704-1725, 1989. **!
      !** **!
      !** Situação: (28-09-2016) OK **!
      !** **!
      !**************************************************************************!
      DOUBLE PRECISION t0 ! [dias] idade inicial para o ajuste
      DOUBLE PRECISION tinf ! [dias] idade final para o ajuste
      INTEGER ndec ! [un] numero de pontos por decada no ajuste
      INTEGER nkelvin ! [un] numero de elementos de Kelvin
      INTEGER ndt ! [un] total de pontos do ajuste
      DOUBLE PRECISION d4 ! auxiliar
      DOUBLE PRECISION d2 ! auxiliar
      !
      ! numero de camadas do modelo de Kelvin Generalizado
      nkelvin=0
      do
          nkelvin=nkelvin+1
          d4=(10.0d0**(nkelvin-1))*(0.01d0*t0)
          if(d4>=(tinf/2.0d0)) then
              exit
          end if
      end do
      if(nkelvin>=7) then
          nkelvin=6
      end if
      !
      ! numero de tempos discretos
      ndt=0
      d2=(t0)*0.1d0
      do
          ndt=ndt+1
          d2=(10.0d0**(1.0d0/ndec))*d2
          if(d2>=tinf-t0) then
              exit
          end if
      end do
      !
      end
      subroutine KELVIN_betactt0(nkelvin,t0,tinf,ndec,ndt,
     & tau,Emu,temp,rh,hf)
      !**************************************************************************!
      !** Subrotina: Kelvin_Generalizado **!
      !** **!
      !** Objetivo: calcula os parâmetros do modelo de Kelvin-Generalizado **!
      !** **!
      !** Referencia: DIAS,M.M. Análise numérica de vigas mistas aço-concreto **!
      !** pelo método dos elementos finitos: efeitos de longa **!
      !** duração. 2013. Dissertação (Mestrado em Engenharia) **!
      !** Programa de Pós Graduação em Engenharia Civil, UFRGS, **!
      !** Porto Alegre. **!
      !** Porto Alegre. **!
      !** **!
      !** BAZANT,Z.P.; PRASANNAN, S. Solidification Theory for **!
      !** concrete creep. I: Formulation, Journal of Enginnering **!
      !** mechanics, v.115, n.8, p.1691-1703, 1989. **!
      !** **!
      !** BAZANT,Z.P.; PRASANNAN, S. Solidification Theory for **!
      !** concrete creep. II: verification and application. **!
      !** Journal of Enginnering Mechanics, v.115, n.8, **!
      !** p.1704-1725, 1989. **!
      !** **!
      !** CEB MODEL CODE 1990 **!
      !** **!
      !** Situação: (29-09-2016) TESTE **!
      !** **!
      !**************************************************************************!
      INTEGER nkelvin ! [un] numero de elementos de Kelvin
      DOUBLE PRECISION t0 ! [dias] idade inicial para o ajuste
      DOUBLE PRECISION tinf ! [dias] idade final para o ajuste
      INTEGER ndec ! [un] numero de pontos por decada no ajuste
      INTEGER ndt ! [un] total de pontos do ajuste
      DOUBLE PRECISION tau(nkelvin) ! [cm²/(kN*dia)] tempos de retardação da cadeia
      DOUBLE PRECISION Emu(nkelvin) ! [kN/cm²] modulos de elasticidade da cadeia
      DOUBLE PRECISION temp ! [oC] temperatura
      DOUBLE PRECISION rh ! [%] umidade relativa
      DOUBLE PRECISION hf ! [cm] espessura ficticia
      !
      DOUBLE PRECISION dti(ndt) ! [dias] intervalos de tempo
      DOUBLE PRECISION fit ! [adm] guarda o coeficiente que será ajustado
      DOUBLE PRECISION C(nkelvin,nkelvin) ! matriz de coeficiente
      DOUBLE PRECISION b(nkelvin) ! vetor de termos independentes
      DOUBLE PRECISION x(nkelvin) ! vetor de incognitas
      INTEGER ntau ! auxiliar
      INTEGER ndti ! auxiliar
      INTEGER ndt0i ! auxiliar
      INTEGER i,j,k ! contador
      DOUBLE PRECISION a1,a2,aa ! auxiliar
      INTEGER er ! variavel que indica erro na solução do sistema
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
      dti(1)=0.1d0*t0
      DO ndti=2,ndt
          dti(ndti)=(10.0d0**(1.0d0/ndec))*dti(ndti-1)
      ENDDO
      dti(ndt)=tinf-t0
      !
      ! Calcula os tempos de retaração
      tau(1) = 0.01d0*t0
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
      END
      subroutine gauss_parcial(a,b,n,x,tol,er)
      !**************************************************************************!
      !** Subrotina: gauss_parcial **!
      !** **!
      !** Objetivo: resolve sistema de equações pelo metodo de Gauss. **!
      !** **!
      !** Referencia: CHAPRA, S. C.; CANALE, R. P. Métodos numéricos para **!
      !** engenharia. 5.ed, São Paulo: McGraw-Hill, 2008. **!
      !** **!
      !** Situação: (26-10-2016) OK **!
      !** **!
      !**************************************************************************!
      DOUBLE PRECISION a(n,n) ! matriz dos coeficientes
      DOUBLE PRECISION b(n) ! vetor de termos independentes
      INTEGER n ! dimensão do sistema
      DOUBLE PRECISION x(n) ! vetor solução
      DOUBLE PRECISION tol ! tolerância
      INTEGER er ! -1, indica que houve erro
      DOUBLE PRECISION s(n) ! vetor que guarda o maior valor da linha
      INTEGER i,j ! contador
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
      !** Subrotina: eliminate **!
      !** **!
      !** Objetivo: faz a eliminação do sistema de Gauss. **!
      !** **!
      !** Referencia: CHAPRA, S. C.; CANALE, R. P. Métodos numéricos para **!
      !** engenharia. 5.ed, São Paulo: McGraw-Hill, 2008. **!
      !** **!
      !** Situação: (26-10-2016) OK **!
      !** **!
      !**************************************************************************!
      DOUBLE PRECISION a(n,n) ! matriz dos coeficientes
      DOUBLE PRECISION s(n) ! vetor que guarda o maior valor da linha
      INTEGER n ! dimensão do sistema
      DOUBLE PRECISION b(n) ! vetor de termos independentes
      DOUBLE PRECISION x(n) ! vetor solução
      DOUBLE PRECISION tol ! tolerância
      INTEGER er ! -1, indica que houve erro
      INTEGER i,j,k ! contador
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
      !** Subrotina: pivot **!
      !** **!
      !** Objetivo: faz a troca de linhas para evitar pivot negativo ou zero **!
      !** **!
      !** Referencia: CHAPRA, S. C.; CANALE, R. P. Métodos numéricos para **!
      !** engenharia. 5.ed, São Paulo: McGraw-Hill, 2008. **!
      !** **!
      !** Situação: (26-10-2016) OK **!
      !** **!
      !**************************************************************************!
      DOUBLE PRECISION a(n,n) ! matriz dos coeficientes
      DOUBLE PRECISION b(n) ! vetor de termos independentes
      DOUBLE PRECISION x(n) ! vetor solução
      DOUBLE PRECISION s(n) ! vetor que guarda o maior valor da linha
      INTEGER n ! dimensão do sistema
      INTEGER k,p ! auxiliar
      INTEGER ii,jj ! contador
      DOUBLE PRECISION dummy ! auxiliar
      DOUBLE PRECISION maior ! auxiliar
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
      !** Subrotina: substitute **!
      !** **!
      !** Objetivo: faz a retrosubstituição do método de Gauss **!
      !** **!
      !** Referencia: CHAPRA, S. C.; CANALE, R. P. Métodos numéricos para **!
      !** engenharia. 5.ed, São Paulo: McGraw-Hill, 2008. **!
      !** **!
      !** Situação: (26-10-2016) OK **!
      !** **!
      !**************************************************************************!
      DOUBLE PRECISION a(n,n) ! matriz dos coeficientes
      INTEGER n ! dimensão do sistema
      DOUBLE PRECISION b(n) ! vetor de termos independentes
      DOUBLE PRECISION x(n) ! vetor solução
      INTEGER i,j ! contador
      DOUBLE PRECISION soma ! auxiliar
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
      !** Subrotina: matinv **!
      !** **!
      !** Objetivo: inverte uma matriz pela técnica de pivotamento **!
      !** **!
      !** Situação: (26-10-2016) OK **!
      !** **!
      !**************************************************************************!
      INTEGER n ! dimensão do sistema
      DOUBLE PRECISION a(n,n) ! matriz dos coeficientes
      DOUBLE PRECISION ainv(n,n) ! matriz inversa
      DOUBLE PRECISION b(n,2*n) ! matriz aumentada
      DOUBLE PRECISION pivot ! pivô
      DOUBLE PRECISION xnum ! auxiliar
      INTEGER i,j,k ! contador
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
      ! Passo 2: mudando o restante da coluno do pivo para 0,
      c adicionando a cada linha um multiplo adequado do pivot
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
      subroutine jacobi(vect,n,np,d,v,nrot)
      !**************************************************************************!
      !** Subrotina: jacobi **!
      !** **!
      !** Objetivo: encontra os autovalores e autovetores de uma matriz Anxn **!
      !** simétrica Anxn através do método de Jacobi **!
      !** **!
      !** Referencia: PRESS, W.H.; TEUKOLSKY, S.A.; VETTERLING, W.T.; FLANNERY **!
      !** B.P. Numerical Recipes in Fortran 77. 2ed. V1. 1997. **!
      !** **!
      !** Situação: (08-01-2017) OK **!
      !** **!
      !**************************************************************************!
      !
      INTEGER n ! tamanho da matriz
      INTEGER np ! tamanho do array que armazena a matriz
      INTEGER nrot ! numero de rotações
      INTEGER nmax ! tamanho máximo da matriz
      PARAMETER (nmax = 500) ! define tamanho máximo 500
      DOUBLE PRECISION vect(6) ! vetor com as componentes de tensao
      DOUBLE PRECISION a(np,np) ! matriz
      DOUBLE PRECISION d(np) ! vetor contendo os autovalores da matriz a
      DOUBLE PRECISION v(np,np) ! matriz contendo os autovetores da matriz a
      INTEGER i,ip,iq,j ! contadores
      DOUBLE PRECISION c,g,h,s,sm,t,tau,theta,tresh ! auxiliares
      DOUBLE PRECISION b(nmax),z(nmax) ! vetores auxiliares
      !
      ! montando a matriz
      a(1,1) = vect(1)
      a(2,2) = vect(2)
      a(3,3) = vect(3)
      a(1,2) = vect(4)
      a(1,3) = vect(5)
      a(2,3) = vect(6)
      a(2,1) = a(1,2)
      a(3,1) = a(1,3)
      a(3,2) = a(2,3)
      !
      ! Inicializando para a matriz identidade
      do ip=1,n
          do iq=1,n
              v(ip,iq)=0d0
          enddo
          v(ip,ip)=1
      enddo
      !
      ! Inicializando b e d para a diagonal de a
      do ip=1,n
          b(ip)=a(ip,ip)
          d(ip)=b(ip)
          z(ip)=0d0
      enddo
      !
      nrot=0
      !
      do i=1,50
          !
          ! Soma dos elementos fora da diagonal
          sm=0
          do ip=1,n-1
              do iq=ip+1,n
                  sm=sm+abs(a(ip,iq))
              enddo
          enddo
          !
          if (sm.eq.0d0) then
              return
          endif
          !
          if (i.LT.4) then
              tresh = 0.2*sm/n**2
          else
              tresh=0d0
          endif
          !
          do ip=1,n-1
              do iq=ip+1,n
                  g=100d0*abs(a(ip,iq))
                  !
                  ! Após quatro varreduras, pula a rotação se o elemento fora da diagonal é pequeno
                  if((i.GT.4).AND.(abs(d(ip))+g.EQ.abs(d(ip))
     & .AND.(abs(d(iq))+g.EQ.abs(d(iq))))) then
                      a(ip,iq)=0d0
                  else if(abs(a(ip,iq)).GT.tresh) then
                      h=d(iq)-d(ip)
                      if (abs(h)+g.eq.abs(h)) then
                          t=a(ip,iq)/h
                      else
                          theta=0.5*h/a(ip,iq)
                          t=1d0/(abs(theta)+sqrt(1d0+theta**2))
                          if (theta.LT.0d0) t=-t
                      endif
                      c=1d0/sqrt(1+t**2)
                      s=t*c
                      tau=s/(1d0+c)
                      h=t*a(ip,iq)
                      z(ip)=z(ip)-h
                      z(iq)=z(iq)+h
                      d(ip)=d(ip)-h
                      d(iq)=d(iq)+h
                      a(ip,iq)=0d0
                      do j=1,ip-1
                          g=a(j,ip)
                          h=a(j,iq)
                          a(j,ip)=g-s*(h+g*tau)
                          a(j,iq)=h+s*(g-h*tau)
                      enddo
                      do j=ip+1,iq-1
                          g=a(ip,j)
                          h=a(j,iq)
                          a(ip,j)=g-s*(h+g*tau)
                          a(j,iq)=h+s*(g-h*tau)
                      enddo
                      do j=iq+1,n
                          g=a(ip,j)
                          h=a(iq,j)
                          a(ip,j)=g-s*(h+g*tau)
                          a(iq,j)=h+s*(g-h*tau)
                      enddo
                      do j=1,n
                          g=v(j,ip)
                          h=v(j,iq)
                          v(j,ip)=g-s*(h+g*tau)
                          v(j,iq)=h+s*(g-h*tau)
                      enddo
                      nrot=nrot+1
                  endif
              enddo
          enddo
          do ip=1,n
              b(ip)=b(ip)+z(ip)
              d(ip)=b(ip)
              z(ip)=0d0
          enddo
      enddo
      pause
      return
      !
      end subroutine
      subroutine eigsrt(d,v,n,np)
      !**************************************************************************!
      !** Subrotina: eigsrt **!
      !** **!
      !** Objetivo: ordena os autovalores e autovetores **!
      !** **!
      !** Referencia: PRESS, W.H.; TEUKOLSKY, S.A.; VETTERLING, W.T.; FLANNERY **!
      !** B.P. Numerical Recipes in Fortran 77. 2ed. V1. 1997. **!
      !** **!
      !** Situação: (08-01-2017) OK **!
      !** **!
      !**************************************************************************!
      !
      INTEGER n,np ! dimensão da matriz
      DOUBLE PRECISION d(np) ! vetor contendo os autovalores
      DOUBLE PRECISION v(np,np) ! matriz contendo os autovetores
      INTEGER i,j,k ! contador
      DOUBLE PRECISION p ! auxiliar
      !
      do i=1,n-1
          k=i
          p=d(i)
          do j=i+1,n
              if(d(j).GE.p) then
                  k=j
                  p=d(j)
              endif
          enddo
          if(k.NE.i)then
              d(k)=d(i)
              d(i)=p
              do j=1,n
                  p=v(j,i)
                  v(j,i)=v(j,k)
                  v(j,k)=p
              enddo
          endif
      enddo
      return
      end subroutine