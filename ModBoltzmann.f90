Module ModBoltzmann

!  Include 'mpif.h'
! ... data type
  INTEGER, PARAMETER :: RFREAL = SELECTED_REAL_KIND(PRECISION(1.0D0)), &
       DP     = SELECTED_REAL_KIND(PRECISION(1.0D0)), &
       WP     = SELECTED_REAL_KIND(PRECISION(1.0D0))


! ... data structure to store rate data
  TYPE t_rate
     REAL(RFREAL) :: Rrate,Rloss
     character(LEN=80) :: Rtype, Rmolecule
     integer :: indexX
  END TYPE t_rate

! ... data structure to store spline data
  TYPE t_spline
     INTEGER :: nxb, nvar,ib
     REAL(RFREAL) :: dxb
     REAL(RFREAL), POINTER :: xb(:)
     REAL(RFREAL), POINTER :: coeffs(:,:,:)
     REAL(RFREAL), POINTER :: yi(:)
     LOGICAL :: extrapolate=.false.
  END TYPE t_spline

! ... data structure to store cross section data
  TYPE t_crossSection
     character(LEN=80) :: type, molecule
     real(8) :: loss
     INTEGER :: nval, indexX
     real(8),pointer :: val(:,:)
     TYPE(t_spline), Pointer :: spline
  end type t_crossSection

  TYPE t_inelastic
     TYPE(t_crossSection), Pointer :: sigma(:)
  end type t_inelastic


  type(t_crossSection) :: tsigma
  type(t_rate),  pointer, dimension(:)  ::rates
  logical, pointer, dimension(:)  :: is3B(:)
  TYPE(t_spline) :: sVals, sOuts, Fvals
  TYPE(t_spline), Pointer :: m_table(:,:,:,:)
  Real(8), Pointer :: m_XOtable(:), m_Ttable(:), m_Etable(:),m_Wtable(:),m_Ntable(:), m_XOcoeffs(:,:,:)
  Real(8) :: m_DXO, m_DT,  m_DE, m_maxT, m_maxO2, m_DW, m_DN, m_maxE, m_minE, m_maxN, m_minN, m_maxW, nui, nui2

  type(t_inelastic),  pointer, dimension(:)  :: sigmaIN
  type(t_crossSection),  pointer, dimension(:)  :: sigmaEL, sigmaEF
  integer :: nsigmaIN(3), nsigmaEL, nsigmaEF, m_newtIter, n, myBrank=0
  integer, PARAMETER :: SHARING_OPTION=2,&      !parameter for ionization
                                                !1=energy sharing; 2=no energy sharing
                        ENERGY_OPTION=1,& 
                        ELECTRON_N_OPTION = 1,& !parameter for electron-electron collision
                                                !1=with e-e collision;
                                                !0=without e-e collision;
                        m_neqin=10,&
                        m_ntarget=1             !number of involved species
  integer,pointer :: IPIV(:)

  real(8), pointer :: ep(:), epm(:), dep(:), W(:),D0(:),D1(:), sepm(:), sigmaE(:), sigmaMM(:), sigmaM(:), sigmaM1(:), F(:), F0(:), &
       igr(:),Nuv(:,:),P0(:),Fl(:),Fr(:),b(:),A(:,:), Vals(:,:), tempVar(:,:), Outs(:,:), A1(:), A2(:), A3(:)
  real(8), pointer :: A1_(:,:), A2_(:,:), A3_(:,:), Fr_(:,:),Fl_(:,:), Te_(:)

  character(len=80), allocatable :: species(:)
  real(8), parameter :: & 
       m_e=1.602176565d-19, &          !electron charge [C]
       m_m = 9.10938215d-31, &         !electron mass [kg]
       m_GAMMA = sqrt(2d0*m_e/m_m), &  !constant=f(e,m)
       m_eps0 = 8.8541878176d-12, &    !permittivity [F/m]		
       m_kB=1.3806488d-23, &           !Boltzmann constant [J/K]
       m_kBeV = 8.6173324d-5, &        !Boltzmann constant [eV/K]				
       m_NA = 6.02214129d26, &         !Avogadro number [particles/kmole]
       m_ev2K = 1.1604522167d4, &      !conversion factor from eV to K
       m_P=1d5, &                      !Pa !Conditions read in by PlascomCM
       m_Td2Vmsq = 1d-21, &            !conversion factor from Td to Vm^2
       m_TOL = 1d-8,&                  !tollerance for Newton method
       m_ng0 = m_P/m_kB/3d2,&          !number density from ideal gas law
       m_PI = ACOS(-1d0),&             !greek pi
       m_giLim = 1d-1                  !limit on logarithmic slope
  real(8) :: m_T=3d2, &                !Temperature [K]
       m_X(m_ntarget), &               !mass fraction
       m_ng, &                         !gas number density
       m_ne, &                         !electron number density
       m_DC, &
       m_RD, &
       m_RT, &
       m_WIONEXP, &
       m_freq, m_delay, m_omega, m_Th, m_ETd_plascom, m_tf,m_tr, &
       m_TIMEint, m_eye(m_neqin,m_neqin), &
       m_scale(2) = [1d0,1d0], m_maxFactor = 1d0,m_facH2, m_kDT=0d0, m_KN2Hp(5), m_KH2Op(2), m_KH3Op(5), m_KH2O3m, &
       m_KOm(5), m_KN2p(1), &
       m_a, m_lambda                   !multiplicative constant & Coulomb logarithm
  REAL(8), PARAMETER :: D_QNAN = TRANSFER((/ Z'00000000', Z'7FF80000' /),1.0_8)
  logical, parameter :: m_debug = .false., m_write=.false., m_pause = .false.
  logical :: m_print = .false.
  integer, parameter :: m_iH2=1, m_iO2=2, m_iN2=3, m_iH2O=4, m_speciesCoupled = m_ntarget
  integer :: m_ikd=0,m_ikaH2=0,m_ikii=0,m_ikiH2=0,m_ikaox=0,m_ikiiox=0,m_ikiox=0, &
       m_ikrox=0,m_ikin=0,m_ikrn=0,m_ikiw=0,m_ikrw=0,m_ikaw=0,m_ikrn2h=0,m_ikrh3o=0,m_ikr=0,&
       m_ikdiss1 =0,m_ikdiss2 =0,m_ikdiss3 =0,m_ikdiss4 =0, m_ikdiss5=0, m_ikdiss6=0, m_ikdiss7=0,  &
       m_itotalLoss =0,m_invals=0, m_inouts=0,m_isvals
  integer :: m_nout, m_islogTab(4)=[1,0,0,0]  !nwto
  real(8), pointer :: m_yin(:)
  character(280) :: m_strVals

!parallel stuff
  integer :: ierr, myrank, numproc,mycomm

CONTAINS

!!!subroutine needed to read the cross section file and save it into variables
  Subroutine readCrossSections
    Implicit None
    character(len=80) :: fn,line
    integer :: i,j,k,kk,kkk, L,i3body
    logical :: found

! ... open the file within the Cross Section
    fn = '../input/argon_tot.tbf'
    open(123,file=trim(fn))


! ... allocate the memory
    allocate(sigmaIN(3))              !Excitation, Ionization and Attachment cross
                                      !section
    do k=1,3;
      allocate(sigmaIN(k)%sigma(0));
    end do
    allocate(sigmaEF(0), sigmaEL(0))  !Effective and elastic cross section
    j=0;
    do
      read(123,*, ERR=10,END=20) line
      tsigma%type = trim(line)               !collision type
      read(123,*, ERR=10,END=20) line
      tsigma%molecule = trim(line)           !involved specie
      j=j+1;
      read(123,*, ERR=10,END=20) tsigma%loss !threshold energy
      read(123,*, ERR=10,END=20) tsigma%nval !number of values
      allocate(tsigma%val(tsigma%nval,2))
      do k =1,tsigma%nval
        read(123,*) tsigma%val(k,1:2)        !energy and cross section vectors
      enddo
      allocate(tsigma%spline)
      read(123,*, ERR=10,END=20) tsigma%spline%nxb   !number of interpolated
                                                     !values
      tsigma%spline%nvar = 1
      allocate(tsigma%spline%xb(tsigma%spline%nxb));  
      do k =1, tsigma%spline%nxb
        read(123,*) tsigma%spline%xb(k)              !interpolated energy vector
      enddo
      tsigma%spline%dxb = tsigma%spline%xb(2)-tsigma%spline%xb(1) !delta energy
      read(123,*) i,k
      allocate(tsigma%spline%coeffs(tsigma%spline%nxb-1,4,tsigma%spline%nvar))
      Do k = 1, tsigma%spline%nxb-1
        Read (123,*) tsigma%spline%coeffs(k,:,1)  !spline polynomial
                                                  !coefficients
      End Do

      if(trim(tsigma%type) == 'EXCITATION') then
        nsigmaIN(1) = pointerT_increment(sigmaIN(1)%sigma,tsigma)
      elseif(trim(tsigma%type) == 'IONIZATION') then
        nsigmaIN(2) = pointerT_increment(sigmaIN(2)%sigma,tsigma)
      elseif(trim(tsigma%type) == 'ATTACHMENT') then
        nsigmaIN(3) = pointerT_increment(sigmaIN(3)%sigma,tsigma)
      elseif(trim(tsigma%type) == 'EFFECTIVE') then
        nsigmaEF = pointerT_increment(sigmaEF,tsigma)
      elseif(trim(tsigma%type) == 'ELASTIC') then
        nsigmaEL = pointerT_increment(sigmaEL,tsigma)
      end if

    end do
10  CONTINUE
20  CONTINUE

    close(123)

    allocate(is3B(nsigmaIN(3)))
    is3B = .false.
    i3body = 1
    allocate(species(m_ntarget))
    if(m_ntarget == 4) then
      species = ['H2 ','O2 ','N2 ','H2O']
    elseif(m_ntarget == 1) then
      !species =['N2']
      species = ['Ar']
    else
      stop ' No. of species urecognized'
    endif
    
    L=0!total counter on rates
    kkloop: do kk=1,size(sigmaIN)
      do kkk=1,size(sigmaIN(kk)%sigma)
        L=L+1
        found = .false.;
        do k = 1,size(species)
        if (trim(sigmaIN(kk)%sigma(kkk)%molecule)== trim(species(k))) then
            sigmaIN(kk)%sigma(kkk)%indexX = k
            found = .true.
            if(kk==3 .and. k==2 .and. i3body==1) then  !the first of the O2
              is3B(kkk) = .true.
              i3body = 0
!!print*,"3Body",sigmaIN(kk)%sigma(kkk)%molecule,sigmaIN(kk)%sigma(kkk)%type
            endif
            exit
          endif
        end do
        if(.not. found) then
!!print*,trim(sigmaIN(kk)%sigma(kkk)%molecule), 'not in ', (trim(species(k)),k=1,3)
          sigmaIN(kk)%sigma(kkk)%indexX = -1
        endif
      enddo
    end do kkloop

    allocate(rates(L))

    nullify(ep, epm, sigmaE, sigmaMM, F, F0, igr, Nuv, P0)

  end Subroutine readCrossSections

!!!subroutine needed to calculate the elastic and the momentum cross section
  Subroutine processCrossSections(Xin,Tin, electronN, Npoints,Ze,newtIter, TinitIn)
    Implicit None
    real(8), intent(IN) ::Xin(:),Tin
    real(8), intent(IN), OPTIONAL :: TinitIn
    integer :: i,j,k,Npoints, newtIter
    real(8) :: Ze,dep1,Te, electronN,n0

    n=Npoints;            !number of grid point
    m_newtIter = newtIter !number of iteration
    m_X = max(Xin,0d0)    !partial mass fraction
    m_X = m_X/sum(m_X);   !total mas fraction
    m_T = Tin             !temperature
    m_ng = m_P/m_kB/m_T   !gas number density
    m_ne = electronN      !electron number density
    n0=1d18;
!!!!print*,m_X,Xin        !print on video the total and partial mass fraction


    m_a = m_e**2d0*m_GAMMA/24d0/m_pi/m_eps0**2d0                      !multiplicative constant
    m_lambda = 12d0*m_pi*(m_eps0*m_kb)**(1.5d0)/m_e**3d0/sqrt(n0)   !Coulomb logarithm

    if(.not. associated(ep)) then


      allocate(ep(n),W(n), epm(n), sepm(n), dep(n), sigmaE(n), sigmaMM(n), sigmaM1(n), sigmaM(n), F(n), F0(n), igr(n), Nuv(n,n), &
           P0(n), D0(n), D1(n), Fl(n), Fr(n), tempVar(n,2), A1(n),A2(n),A3(n))
      allocate(A1_(n,n), A2_(n,n), A3_(n,n), Fr_(n,n), Fl_(n,n), Te_(n))

      A1=0d0;A2=0d0;A3=0d0;           !energy integrals-eq.34,35,36 Hagelaar article
      A1_=0d0;A2_=0d0;A3_=0d0;        !derivate of integrals for Jacobian matrix
      allocate(IPIV(n), A(n,n), b(n)) !parameter for solve the linear system
      Fvals%nvar = 1
      Fvals%nxb =  ubound(F,1)
      Fvals%ib = 1
      Fvals%extrapolate = .true.
      allocate(Fvals%xb(Fvals%nxb),Fvals%yI(Fvals%nvar));
      allocate(Fvals%coeffs(Fvals%nxb,4,Fvals%nvar))
    endif
    A=0d0
    b=0d0
    F = 0d0

! ... definition of the energy grid
    dep1=Ze/dble(n-1)          !deltaE [eV]
    ep=(/((k-1)*dep1,k=1,n)/)  !e_i+1/2 with first element=0
    epm = ep + dep1/2d0;         !e_i shifted of deltaE/2, first element=/0
    ep = ep + dep1;            !e_i+1/2 shifted of deltaE, first element=/0
    sepm = sqrt(epm)
    tempVar(:,1) = ep;

    do k =1,n
      if (k==1) then
        igr(k) = ep(k)*sepm(k)           !parameter usefull for the integration
      else                               !of the distribution function
        igr(k) = (ep(k)-ep(k-1))*sepm(k)
      end if
      if(k<n) then
        dep(k) = (ep(k+1)-ep(k))         !update of deltaE
      else
        dep(k) = (ep(k)-ep(k-1))
      end if
    enddo

    sigmaM1 = 1d0/(m_GAMMA*sqrt(ep)) !coefficient from eq.11 Hagelaar article
    sigmaE = 0d0
    sigmaMM = 0d0

    if(maxval(abs(F)) < 1d-20) then
      if(.not. present(TinitIn)) then
        Te = 1d0*m_ev2K;
      else
        Te = TinitIn
      endif
      F=exp(-m_e/m_kB/Te*epm);    !Initial distribution: Maxwellian distribution
      F=F/sum(igr*F);             !normalization of F
    endif


    do k = 1,ubound(species,1)
      do j = 1,ubound(sigmaEF,1)

        if (trim(species(k)) == trim(sigmaEF(j)%molecule) ) then
          do i = 1,ubound(ep,1)
            !total momentum transfer cross section, eq.7 Hagelaar article
            sigmaMM(i) = sigmaMM(i) + evaluateSplineExp(sigmaEF(j)%spline%dxb, &
                 sigmaEF(j)%spline%xb, sigmaEF(j)%spline%coeffs, sigmaEF(j)%spline%nvar, ep(i))*m_X(k);
            !elastic momentum transfer cross section, eq.42 Hagelaar article
            sigmaE(i) = sigmaE(i) + evaluateSplineExp(sigmaEL(j)%spline%dxb, &
                 sigmaEL(j)%spline%xb, sigmaEL(j)%spline%coeffs, sigmaEL(j)%spline%nvar, ep(i))*m_X(k);
           ! sigmaE(i)=0d0;
          enddo
        end if

      end do
    end do

    W=-m_GAMMA* ep**2d0*sigmaE; !first addend eq.40 Hagelaar article
    D0=-W*m_kB*Tin/m_e         !second addend eq.41 Hagelaar article
    tempVar(:,2) = sigmaE;

  end Subroutine processCrossSections

!!!subroutine needed to solve the Boltzmann equation
  Subroutine solveBoltzmann(E_N,isTdin,Meanenergy,mean_energy_K,sigmaM)
    use solver, only: DGESV
    Implicit None
    real(8), intent(IN) ::E_N
    logical, intent(IN), optional :: isTdin !isTdin is a logical variable
                                            !true=if E_N is in Td;false=if E_N is Vm^2
    real(8), intent(OUT) ::Meanenergy,mean_energy_K,sigmaM(n)

    integer :: i,ii,j,k,kk,kkk, iteration, endV, iL, iH, iH0
    real(8) :: snu(3), factor, loss, epsI, epsJ,gi,gj,eps1,eps2, nA, pik, E_Nsq, oldval, ionFactor(2)   !,nu
    real(8) :: D, z, epz, emz, omemz, omepz, uk, epsj1, epsj2, epsi1, epsi2, Te, epsH, epsL, sigmaH
    real(8) :: sigmaL, errorF, errorO, smallA, lambda,Wtilde, logLambda
    real(8),parameter :: smallNo=1d-99

!!!!!print*,'m_ne=',m_ne,',m_a=',m_a,',m_lambda=',m_lambda

    E_Nsq = E_N*E_N
    if(present(isTdin)) then
      if(isTdin) E_Nsq = (E_N*m_Td2Vmsq)**2d0
    endif
print*,E_Nsq
    if(m_debug) then
      print*,(trim(species(k)),k=1,m_ntarget);
      print*,'epm,      F,    F0'
      do k =1,n
        print*,epm(k),F(k), F0(k)
      end do

    endif

    if(isTdin .and. E_N < 1d-3) then
      F=exp(-m_e/m_kB/m_T*epm);
      F=F/sum(igr*F);
      return
    endif

    errorO = 1d99

! ... start of the Newton loop for the solution of a non linear system
    iterationLoop: do iteration=1,m_newtIter

      F0=F;
      A = 0d0
      b=0d0

      Nuv=0d0;
      P0 = 0d0;
      snu = [0d0,1d0,-1d0];
      nui2=0d0

      kkloop: do kk=1,size(sigmaIN)
        do kkk=1,size(sigmaIN(kk)%sigma)

          k = sigmaIN(kk)%sigma(kkk)%indexX
          if(k <= 0) cycle
          factor = m_GAMMA*m_X(k);
          !if(kk == 3 .and. is3B(kkk))  factor=factor*m_ng/m_ng0
          endV = ubound(sigmaIN(kk)%sigma(kkk)%val,1)

          if (kk < 3) then
            loss=sigmaIN(kk)%sigma(kkk)%loss;
          else
            loss = 0d0  !kk=3 --> attachment
          end if
          do i = 1,n
            epsI = epm(i);
            if (i == 1) then
              eps1=0d0;
              gi = -log(F(i+1)/F(i))/(epm(i+1)-epm(i));
            else
              eps1=ep(i-1);
              if (i < n .and. F(min(i+1,n)) > smallNo .and. F(i-1) > smallNo) then
                gi = -log(F(i+1)/F(i-1))/(epm(i+1)-epm(i-1)); !logarithmic slope
              elseif(F(i-1) > smallNo) then                   !eq.51 Hagelaar article
                 gi = -log(F(i)/F(i-1))/(epm(i)-epm(i-1));
              elseif(i < n .and. F(min(i+1,n)) > smallNo ) then
                gi = -log(F(i+1)/F(i))/(epm(i+1)-epm(i));
              else
                gi = 0d0
              end if
            end if
!if (isnan(gi)  .or. .not. abs(gi) < 1d99) cycle
            eps2=ep(i);
            if (kk <= 2 .and. eps2 < loss) cycle;
!
            iL=-1;iH=-1;
            do ii = 1,ubound(sigmaIN(kk)%sigma(kkk)%val,1)
              if(sigmaIN(kk)%sigma(kkk)%val(ii,1) < eps1) iL=ii
              if(sigmaIN(kk)%sigma(kkk)%val(ii,1) > eps2) then
                iH = ii;exit
              end if
            enddo

            if (iL < 0 .and. epsI > sigmaIN(kk)%sigma(kkk)%val(1,1)) iL=1
            if (iH < 0 .and. epsI <= sigmaIN(kk)%sigma(kkk)%val(endV,1)) iH=size(sigmaIN(kk)%sigma(kkk)%val(:,1))
            if (iL < 0 .or. iH < 0 .or. iL == iH) cycle


            pik=0d0;
            iH0 = iH
            if (iH > iL+1) then
              do while (iL+1 < iH0) 
                iH = iL+1;
                !cross section assumed to be linear in between the points in a
                !table of cross section versus energy
                epsH=sigmaIN(kk)%sigma(kkk)%val(iH,1);
                sigmaH=sigmaIN(kk)%sigma(kkk)%val(iH,2);
                epsL=sigmaIN(kk)%sigma(kkk)%val(iL,1);
                sigmaL=sigmaIN(kk)%sigma(kkk)%val(iL,2);

                oldval = eps2;
                eps2=epsH;  
                if (eps1 < loss .and. eps2 > loss) then
                  sigmaL = 0d0;
                  epsL = loss
                  eps1=loss;
                elseif (epsL < loss .and. epsH > loss) then
                  sigmaL = 0d0;
                  epsL = loss
                end if
                if(gi > m_giLim) then
                  !integral eq.47 Hagelaar article
                  pik=pik+(exp((-eps2 + epsI)*gi)*((-2d0 + gi*(epsL + eps2*(-2d0 - eps2*gi + epsL*gi)))* &
                       sigmaH + (2d0 + gi*(2d0*eps2 - epsH + eps2*(eps2 - epsH)*gi))*sigmaL) +   &
                       exp((-eps1 + epsI)*gi)*((2d0 + gi*(2d0*eps1 - epsL + eps1*(eps1 - epsL)*gi))* &
                       sigmaH + (-2d0 + gi*(epsH + eps1*(-2d0 - eps1*gi + epsH*gi)))*sigmaL))/  &
                       ((epsH - epsL)*gi**3d0)
                else
                  pik=pik+(2d0*eps1**3d0*(-sigmaH + sigmaL) + 3d0*eps1**2d0*(epsL*sigmaH - epsH*sigmaL) + &
                       eps2**2d0*(-3d0*epsL*sigmaH + 2d0*eps2*(sigmaH - sigmaL) + 3d0*epsH*sigmaL))/&
                       (6d0*(epsH - epsL))
                endif

                iH = iH+1;
                iL = iH-1;
                epsH=sigmaIN(kk)%sigma(kkk)%val(iH,1);
                sigmaH=sigmaIN(kk)%sigma(kkk)%val(iH,2);
                epsL=sigmaIN(kk)%sigma(kkk)%val(iL,1);
                sigmaL=sigmaIN(kk)%sigma(kkk)%val(iL,2);
                eps1=epsL;
                eps2=epsH;
              end do
              eps2=oldval;
            else
              epsH=sigmaIN(kk)%sigma(kkk)%val(iH,1);
              sigmaH=sigmaIN(kk)%sigma(kkk)%val(iH,2);
              if (iL > 0) then
                epsL=sigmaIN(kk)%sigma(kkk)%val(iL,1);
                sigmaL=sigmaIN(kk)%sigma(kkk)%val(iL,2);
              end if
            end if

            if (eps1 < loss .and. eps2 > loss) then
              eps1=loss;
              sigmaL = 0d0;
              epsL = loss;
            elseif (epsL < loss .and. epsH > loss) then
              sigmaL = 0d0;
              epsL = loss
            end if
!
            if(gi > m_giLim) then
              pik=pik+(exp((-eps2 + epsI)*gi)*((-2d0 + gi*(epsL + eps2*(-2d0 - eps2*gi + epsL*gi)))* &
                   sigmaH + (2d0 + gi*(2d0*eps2 - epsH + eps2*(eps2 - epsH)*gi))*sigmaL) +   &
                   exp((-eps1 + epsI)*gi)*((2d0 + gi*(2d0*eps1 - epsL + eps1*(eps1 - epsL)*gi))* &
                   sigmaH + (-2d0 + gi*(epsH + eps1*(-2d0 - eps1*gi + epsH*gi)))*sigmaL))/  &
                   ((epsH - epsL)*gi**3d0)
            else
              pik=pik+(2d0*eps1**3d0*(-sigmaH + sigmaL) + 3d0*eps1**2d0*(epsL*sigmaH - epsH*sigmaL) + &
                   eps2**2d0*(-3d0*epsL*sigmaH + 2d0*eps2*(sigmaH - sigmaL) + 3d0*epsH*sigmaL))/&
                   (6d0*(epsH - epsL))
            end if

            if(pik < 0d0)cycle
            A(i,i) = A(i,i) + factor*pik;
            Nuv(i,i)=Nuv(i,i) - factor*pik;
            nui2 = nui2 + factor*pik*F(i)*snu(KK)
            if (SHARING_OPTION == 2 .and. kk ==2) then !only for ionization, if the energy is not shared 
              P0(i) = P0(i)+factor*pik;
            end if
          end do
        end do
      end do kkloop

      if(m_debug) then
        write(465,'(1p9923D14.4)')(A(k,k),k=1,n)
        flush (465)
      end if

      

      if (SHARING_OPTION == 1) then
        ionFactor = [1d0,2d0];
      else
        ionFactor = [1d0,1d0];
      end if

      kkloopQ: do kk=1,size(sigmaIN)-1 !no attachement
        do kkk=1,size(sigmaIN(kk)%sigma)

          k = sigmaIN(kk)%sigma(kkk)%indexX
          if(k <= 0) cycle
          factor = m_GAMMA*m_X(k);

          uk=sigmaIN(kk)%sigma(kkk)%loss;
          loss=sigmaIN(kk)%sigma(kkk)%loss;
          endV = ubound(sigmaIN(kk)%sigma(kkk)%val,1)
          jloopQ: do j = 1,n

            if (j == 1) then
              epsj1=0d0;
              gj = -log(F(j+1)/F(j))/(epm(j+1)-epm(j));
            else
              epsj1=ep(j-1);
              if (j < n .and. F(min(j+1,n)) > smallNo .and. F(j-1) > smallNo) then
                gj = -log(F(j+1)/F(j-1))/(epm(j+1)-epm(j-1));
              elseif(F(j-1) > smallNo) then
                gj = -log(F(j)/F(j-1))/(epm(j)-epm(j-1));
              elseif(j < n .and. F(min(j+1,n)) > smallNo ) then
                gj = -log(F(j+1)/F(j))/(epm(j+1)-epm(j));
              else
                gj=0d0
              end if
            end if
!            if (isnan(gj)  .or. .not. abs(gj) < 1d99) cycle

            epsj2=ep(j);
            epsJ = epm(j);

            iloopQ: do i = 1,n

              if (i == 1)then
                epsi1=0d0;
              else
                epsi1=ep(i-1);
              end if
              !energy points defined by eq.49,50 Hagelaar article
              epsi1 = ionFactor(kk)*epsi1 + uk;
              epsi2 = ionFactor(kk)*ep(i) + uk;
              eps1 = max(epsi1, epsj1);
              eps2 = min(epsi2, epsj2);

              if (eps2 <= eps1) cycle
              if (eps2 < uk) cycle
!
              iL=-1;iH=-1;
              do ii = 1,ubound(sigmaIN(kk)%sigma(kkk)%val,1)
                if(sigmaIN(kk)%sigma(kkk)%val(ii,1) < eps1) then
                  iL=ii
                  cycle
                end if
                if(sigmaIN(kk)%sigma(kkk)%val(ii,1) > eps2) then
                  iH = ii;exit
                end if
              enddo
              if (iL < 0 .and. epsJ > sigmaIN(kk)%sigma(kkk)%val(1,1)) iL=1
              if (iH < 0 .and. epsJ < sigmaIN(kk)%sigma(kkk)%val(endV,1)) iH=size(sigmaIN(kk)%sigma(kkk)%val(:,1))
              if (iL < 0 .or. iH < 0 .or. iL == iH) cycle


              pik=0d0;
              iH0 = iH
              if (iH > iL+1) then
                do while (iL+1 < iH0) 
                  iH = iL+1;
                  epsH=sigmaIN(kk)%sigma(kkk)%val(iH,1);
                  sigmaH=sigmaIN(kk)%sigma(kkk)%val(iH,2);
                  epsL=sigmaIN(kk)%sigma(kkk)%val(iL,1);
                  sigmaL=sigmaIN(kk)%sigma(kkk)%val(iL,2);
                  oldval = eps2;
                  eps2=epsH;  
                  if (eps1 < loss .and. eps2 > loss) then
                    sigmaL = 0d0;
                    epsL = loss
                    eps1=loss;
                  elseif (epsL < loss .and. epsH > loss) then
                    sigmaL = 0d0;
                    epsL = loss
                  end if
                  if(gj  > m_giLim) then
                  !integral eq.48 Hagelaar article 
                    pik=pik+(exp((-eps2 + epsJ)*gj)*((-2d0 + gj*(epsL + eps2*(-2d0 - eps2*gj + epsL*gj)))* &
                         sigmaH + (2d0 + gj*(2d0*eps2 - epsH + eps2*(eps2 - epsH)*gj))*sigmaL) +   &
                         exp((-eps1 + epsJ)*gj)*((2d0 + gj*(2d0*eps1 - epsL + eps1*(eps1 - epsL)*gj))* &
                         sigmaH + (-2d0 + gj*(epsH + eps1*(-2d0 - eps1*gj + epsH*gj)))*sigmaL))/  &
                         ((epsH - epsL)*gj**3d0)
                  else
                    pik=pik+(2d0*eps1**3d0*(-sigmaH + sigmaL) + 3d0*eps1**2d0*(epsL*sigmaH - epsH*sigmaL) + &
                         eps2**2d0*(-3d0*epsL*sigmaH + 2d0*eps2*(sigmaH - sigmaL) + 3d0*epsH*sigmaL))/&
                         (6d0*(epsH - epsL))
                  endif
                  iH = iH+1;
                  iL = iH-1;
                  epsH=sigmaIN(kk)%sigma(kkk)%val(iH,1);
                  sigmaH=sigmaIN(kk)%sigma(kkk)%val(iH,2);
                  epsL=sigmaIN(kk)%sigma(kkk)%val(iL,1);
                  sigmaL=sigmaIN(kk)%sigma(kkk)%val(iL,2);
                  eps1=epsL;
                  eps2=epsH;
                end do
                eps2=oldval;
              else
                epsH=sigmaIN(kk)%sigma(kkk)%val(iH,1);
                sigmaH=sigmaIN(kk)%sigma(kkk)%val(iH,2);
                if (iL > 0) then
                  epsL=sigmaIN(kk)%sigma(kkk)%val(iL,1);
                  sigmaL=sigmaIN(kk)%sigma(kkk)%val(iL,2);
                end if
              end if
              if (eps1 < loss .and. eps2 > loss) then
                eps1=loss;
                sigmaL = 0d0;
                epsL = loss;
              elseif (epsL < loss .and. epsH > loss) then
                sigmaL = 0d0;
                epsL = loss
              end if
!
              if(gj  > m_giLim) then
                pik=pik+(exp((-eps2 + epsJ)*gj)*((-2d0 + gj*(epsL + eps2*(-2d0 - eps2*gj + epsL*gj)))* &
                     sigmaH + (2d0 + gj*(2d0*eps2 - epsH + eps2*(eps2 - epsH)*gj))*sigmaL) +   &
                     exp((-eps1 + epsJ)*gj)*((2d0 + gj*(2d0*eps1 - epsL + eps1*(eps1 - epsL)*gj))* &
                     sigmaH + (-2d0 + gj*(epsH + eps1*(-2d0 - eps1*gj + epsH*gj)))*sigmaL))/  &
                     ((epsH - epsL)*gj**3d0)
              else
                pik=pik+(2d0*eps1**3d0*(-sigmaH + sigmaL) + 3d0*eps1**2d0*(epsL*sigmaH - epsH*sigmaL) + &
                     eps2**2d0*(-3d0*epsL*sigmaH + 2d0*eps2*(sigmaH - sigmaL) + 3d0*epsH*sigmaL))/&
                     (6d0*(epsH - epsL))
              endif
              if(pik < 0d0) cycle
              A(i,j)=A(i,j) -factor*pik;
              Nuv(i,j)=Nuv(i,j) + factor*pik;
            end do iloopQ
         end do jloopQ
        end do
      end do kkloopQ


      if (SHARING_OPTION == 2) then
        A(1,:) = A(1,:) - P0;
        Nuv(1,:) = Nuv(1,:) + P0;
      end if


      nui=sum(matmul(Nuv,F));

      do i = 1,n
          A(i,i) = A(i,i) + igr(i) * nui
       do j =1,n
          A(i,j) = A(i,j) + igr(i) * F(i) *sum(Nuv(:,j))
      enddo
         b(i) = b(i) + igr(i) * F(i) * nui
   enddo

      if(ELECTRON_N_OPTION == 1) then
        !Te = min(max(getTe(),m_T/10d0),30d0*m_ev2K);
        call getTe(mean_energy_K,Meanenergy)
        Te = min(max(mean_energy_K,m_T/10d0),30d0*m_ev2K);
        call getTeDer                  !evaluation of temperature derivate for
                                       !Jacobian matrix
        Call EvaluateA1A2A3            !evaluation of eq.34,35,36,Hagelaar article
        lambda = m_lambda*Te**(1.5d0)
        logLambda = log(lambda)
        !print*,'coulomb logarithm',logLambda
        smallA=m_a*logLambda
        Te_ = 1.5d0*m_a/Te*Te_!this is smallA_
        A1 = -A1*3d0*smallA*(m_ne/m_ng)              !second addend of eq.40,Hagelaar article
        A2=A2*2d0* smallA* (m_ne/m_ng)               !third addend of eq.41,Hagelaar article
        A3=A3*2d0* smallA* (m_ne/m_ng) * ep**(1.5d0) !fourth addend of eq.41,Hagelaar article

        A1_ = -A1_*3d0*smallA*(m_ne/m_ng)            !derivative of A1, for Jacobian matrix
        A2_=A2_*2d0* smallA* (m_ne/m_ng)             !derivative of A2, for Jacobian matrix
        do k=1,n                                     !derivative of A3, for Jacobian matrix
          A3_(k,:)=A3_(k,:)*2d0* smallA* (m_ne/m_ng) * ep(k)**(1.5d0)
        enddo
! Te variation
        if(smallA > 0d0) then
          do k=1,n
            A1_(k,:) = A1_(k,:) + A1(k)/smallA*Te_(:)
            A2_(k,:) = A2_(k,:) + A2(k)/smallA*Te_(:)
            A3_(k,:) = A3_(k,:) + A3(k)/smallA*Te_(:)
          enddo
        end if
      endif

      sigmaM = sigmaMM + nui*sigmaM1; !eq.12, Hagelaar article
      D1 = m_GAMMA/3d0*ep/sigmaM     !first addend of eq.41,Hagelaar article
      do k=1,n
        if(ELECTRON_N_OPTION == 0) then
          D = D1(k)*E_Nsq + D0(k);
          Wtilde = W(k)
        else
          D = D1(k)*E_Nsq + D0(k) + A2(k) + A3(k); !eq.41,Hagelaar article
          Wtilde = W(k)+A1(k)                      !eq.40,Hagelaar article 
                                                   !plus A1, becasue W is negative
        endif
        z=Wtilde*dep(k)/D;     !eq. in the text pag 727,Hagelaar article
        epz=exp(z);            !positive exp in eq.45,Hagelaar article
        emz=exp(-z);           !negative exp in eq.45,Hagelaar article
        omemz = 1d0-emz;       !denominator in eq.45,Hagelaar article
        omepz = 1d0-epz;       !denominator in eq.45,Hagelaar article
        Fr(k) = Wtilde/omepz;  !second addend coefficient in eq.45,Hagelaar article
        Fl(k) = Wtilde/omemz;  !first addend coefficient in eq.45,Hagelaar article
        if(ELECTRON_N_OPTION == 1) then
          do kk=1,n
            !derivative of Fr and Fl, for Jacobian matrix
            Fr_(k,kk) =-(((-1d0 + epz)*A1_(k,kk) + (dep(k)*epz*Wtilde*(-(D*A1_(k,kk)) + &
                 Wtilde*(A2_(k,kk) + A3_(k,kk))))/D**2d0)/(-1d0 + epz)**2d0)
            Fl_(k,kk) =  ((1d0 - emz)*A1_(k,kk) + (dep(k)*Wtilde*(-(D*A1_(k,kk)) + &
                 Wtilde*(A2_(k,kk) + A3_(k,kk))))/(D**2d0*epz))/(-1d0 + emz)**2d0
          enddo
        endif
      end do

      k=1
      A(k,k) = A(k,k)+Fl(k)
      A(k,k+1) = A(k,k+1)+ Fr(k)
      do k=2,n-1
        A(k,k) = A(k,k) + Fl(k) - Fr(k-1)
        A(k,k-1) = A(k,k-1) - Fl(k-1)
        A(k,k+1) = A(k,k+1) + Fr(k)
      enddo


      if(ELECTRON_N_OPTION == 1) then
        k=1
        kk = k
        A(k,:) = A(k,:)+Fl_(k,:)*F(kk)
        b(k) = b(k) + sum(Fl_(k,:)*F(:))*F(kk)
        kk = k+1
        A(k,:) = A(k,:)+ Fr_(k,:)*F(kk)
        b(k) = b(k) + sum(Fr_(k,:)*F(:))*F(kk)
        do k=2,n-1
          kk = k
          A(k,:) = A(k,:) + (Fl_(k,:) - Fr_(k-1,:))*F(kk)
          b(k) = b(k) + sum((Fl_(k,:) - Fr_(k-1,:))*F(:))*F(kk)
          kk = k-1
          A(k,:) = A(k,:) - Fl_(k-1,:)*F(kk)
          b(k) = b(k) - sum(Fl_(k-1,:)*F(:))*F(kk)
          kk = k+1
          A(k,:) = A(k,:) + Fr_(k,:)*F(kk)
          b(k) = b(k) + sum(Fr_(k,:)*F(:))*F(kk)
        enddo
      endif




!norm of A
      nA = 0d0
      do j=1,n
        nA = max(nA,sum(abs(A(:,j))) )
      enddo


      A(n,:) = igr*nA;

      b(n)=1d0*nA;


      if(m_debug) then
        do i = 1,n
          write(355,'(1p9923D15.5)')A(i,:)
        enddo
        flush(355)
!pause "A matrix"
      end if

!print*,i
      CALL DGESV( N, 1, A, ubound(A,1) , IPIV, B, ubound(B,1), i )
      if(m_debug) print'("FLAG ",i5,1p9999e12.4)',i,B(1:5)


      F= F + (B-F)*min(dble(iteration)/5d0,1d0)
      P0 = F
      F(1:n) = abs(F(1:n))
      do i = 2,n;
        if(P0(i)<0d0) F(i)=F(i)/abs(dble((i-1d0)**2d0));
      enddo
      if(m_debug) then
        write(455,'(1p9923D16.7)')F
        flush (455)
      end if

! ... if cycle to add a Maxwellian decay if the distribution function is too small

      ii=-1
     do i = 1,n-1
       if( ( F(i+1)-F(i) > 0 .and. (F(i+1)+F(i))/2d0 < 1d-3*maxval(F)) .or. max(P0(i), P0(i+1)) < 0d0 )then
          ii=i+1
          exit
        end if
      end do

!print*,'ii',ii
      if (ii>0) then
        F(ii+1:) = 0d0;
        F=F/sum(igr*F);
        !Te = max(getTe(),m_T);
        Te = max(mean_energy_K,m_T);
        do k = ii+1,n
print*,F(ii), m_e, m_kB, Te, epm(k), epm(ii)
!pause
          F(k) = F(ii)*dexp(-m_e/m_kB/Te*(epm(k)-epm(ii)));
print*,F(k)
     !     print*,epm(k)-epm(ii)
        enddo
      !pause
      endif
      !print*,'e,kb,Te',m_e,m_kB,Te
      !print*,'epm',epm
      !pause
      F(1:n) = abs(F(1:n))
      !print*,F(1:n)
 !      pause


      F=F/sum(igr*F);
      if(m_debug) then
        write(455,'(1p9923D14.4)')F
        flush (455)
!pause"->F"
      end if

! ... check on the error for the Newton method
      if(ii <=1) then
        errorF = maxval(abs(F-F0))
      else
        errorF = maxval(abs(F(1:ii)-F0(1:ii)))
      end if
      if(m_write)print'(i4,1p123e12.4)', iteration,E_N, errorF, Te
      if (errorF < m_TOL)exit

      if(errorF > errorO .and. m_pause .and. m_debug) then
        do i = 1,n;
          print'(1p123e12.4)',epm(i),F(i),P0(i),F0(i);
        enddo
      !  pause "errorF > errorO "
      endif

      errorO = errorF

    end do iterationLoop

    if(iteration >= m_newtIter) print*,'Iteration Count Exceeded'

    if(m_write) print'(i4,1p123e12.4)', iteration,E_N, errorF

! ... if cycle to add a Maxwellian decay if the distribution function is too small

    ii=-1
    do i = 1,n-1
!do i = 1,n
!if( F(i+1)-F(i) > 0 .and. (min(P0(i), P0(i+1)) < 0 .or.  (F(i+1)+F(i))/2d0 < 1e-3) )then
     if( (epm(i) > 4d0 .and. F(i+1)-F(i) > 0 .and. (F(i+1)+F(i))/2d0 < 1d-3) .or. min(P0(i), P0(i+1)) < 0d0 )then
        ii=i+1
        exit
      end if
    end do

    if (ii>0) then
      F(ii+1:) = 0d0;
      F=F/sum(igr*F);
     ! Te = max(getTe(),m_T);
     Te = max(mean_energy_K,m_T);
!Te = 1d0*m_ev2K;
      do k = ii+1,n
        F(k) = F(ii)*exp(-m_e/m_kB/Te*(epm(k)-epm(ii)));
      enddo
      F=F/sum(igr*F);
    end if
    if(m_debug) then
      call getrate(rates)
      write(455,'(1p9923D16.7)')P0
      write(455,'(1p9923D16.7)')F
      flush (455)
      if(m_pause) then
        print*,'Pausing for',E_N,'Rates below'
        do i = 1, size(rates)
          if (trim(rates(i)%Rtype) /= 'EXCITATION') then
            print*,trim(rates(i)%Rmolecule), ' ',trim(rates(i)%Rtype),rates(i)%Rrate;
          endif
        enddo
       ! pause 'End solveBoltzmann'
      endif
    end if

  end Subroutine solveBoltzmann

!!!subroutine needed to evaluate elastic and momentum croos section through spline
  function evaluateSpline(dxb, xb, coeffs, nvar, x)

    Implicit None

! ... global variables
    Integer :: nvar
    Real(RFREAL) :: dxb, x, evaluateSpline, xx
    Real(RFREAL), Pointer :: xb(:), coeffs(:,:,:)

! ... local variables
    Integer :: index, k
    Real(RFREAL) :: val

! ... find lower bound of breakpoint interval
    index = floor((x-xb(1))/dxb)+1

! ... local coordinate
    xx = x - xb(index)

! evaluate spline
    val = coeffs(index,1,nvar)
    do k = 2, 4
      val = val * xx + coeffs(index,k,nvar)
    end do

    evaluateSpline = val

  end function evaluateSpline

!!!subroutine needed to add all the information about the cross section
  function evaluateSplineExp(dxb, xb, coeffs, nvar, x)

    Implicit None

! ... global variables
    Integer :: nvar
    Real(RFREAL) :: dxb, x, evaluateSplineExp, xx
    Real(RFREAL), Pointer :: xb(:), coeffs(:,:,:)

! ... local variables
    Integer :: index, k
    Real(RFREAL) :: val, lx

! ... find lower bound of breakpoint interval
    lx=log(max(x,1d-40))
    index = floor((lx-xb(1))/dxb)+1

! ... local coordinate
    xx = lx - xb(index)

! evaluate spline
    val = coeffs(index,1,nvar)
    do k = 2, 4
      val = val * xx + coeffs(index,k,nvar)
    end do

    evaluateSplineExp = exp(val)

  end function evaluateSplineExp

!!!subroutine needed to add all the information about the cross section
  function pointerT_increment(ptr,val) result(k)
    type(t_crossSection),pointer   :: ptr(:)
    type(t_crossSection) :: val
    type(t_crossSection),allocatable   :: tmp(:)
    integer::k

    k = ubound(ptr,1)
    allocate(tmp(k))
    tmp(1:k) = ptr(1:k)
    allocate(ptr(k+1))
    ptr(1:k) = tmp(1:k)
    k = ubound(ptr,1)
    ptr(k) = val
    deallocate(tmp)

  end function pointerT_increment

!!!subroutine needed to calculated the Temperature usign the eq.37, Hagelaar article
  subroutine getTe(mean_energy_K, Meanenergy)
    Implicit None
    integer :: i
    real(8), INTENT(OUT) :: mean_energy_K, Meanenergy
    real(8) :: factor, epsI, eps1, eps2, gi,agi, rooteps1gi, rooteps2gi, pik, erfrooteps1gi, erfrooteps2gi

    Meanenergy=0d0;
    do i = 1,n
      factor = F(i);
      epsI = epm(i);
      if(F(i) <= 0d0 .or. F(min(i+1,n)) <= 0d0 .or. F(max(i-1,1)) <= 0d0) cycle
      if (i == 1) then
        eps1=0d0;
        gi = log(F(i+1)/F(i))/(epm(i+1)-epm(i));
      else
        eps1=ep(i-1);
        if (i < n .and. F(min(i+1,n)) > 0d0) then
          gi = -log(F(i+1)/F(i-1))/(epm(i+1)-epm(i-1));
        else
          gi = -log(F(i)/F(i-1))/(epm(i)-epm(i-1));
        end if
      end if
      gi = min(abs(gi),7d0)
      agi = abs(gi);
      eps2=ep(i);
      rooteps1gi=sqrt(abs(eps1*gi));
      rooteps2gi=sqrt(abs(eps2*gi));
      erfrooteps1gi = erf(rooteps1gi)
      erfrooteps2gi = erf(rooteps2gi)
      pik =  (2d0*Sqrt(gi)*(exp(eps2*gi)*Sqrt(eps1)*(3d0 + 2d0*eps1*gi) - &
           exp(eps1*gi)*Sqrt(eps2)*(3d0 + 2d0*eps2*gi)) + &
           3d0*exp((eps1 + eps2)*gi)*Sqrt(M_Pi)*(-erfrooteps1gi + erfrooteps2gi))/&
           (4d0*exp((eps1 + eps2 - epsI)*gi)*gi**2.5d0)
      Meanenergy=Meanenergy+factor*pik;
    end do

    mean_energy_K=Meanenergy*2d0/3d0*m_ev2K;

    if(isnan(mean_energy_K) ) stop

  end subroutine getTe

  subroutine getTeDer
    Implicit None
    integer :: i,m,i_(3)
    real(8) :: Meanenergy, factor, epsI, eps1, eps2, gi,agi, rooteps1gi, rooteps2gi, &
         pik, erfrooteps1gi, erfrooteps2gi, gi_(3),pik_

    Meanenergy=0d0;
    Te_=0d0
    do i = 1,n
      factor = F(i);
      epsI = epm(i);
      if(F(i) <= 0d0 .or. F(min(i+1,n)) <= 0d0 .or. F(max(i-1,1)) <= 0d0) cycle
      gi_ = 0d0
      i_=[i-1,i,i+1]
      if (i == 1) then
        eps1=0;
        gi = log(F(i+1)/F(i))/(epm(i+1)-epm(i));
        gi_(2) = 1d0/(epm(i+1)-epm(i))/F(i)
        gi_(3) = -1d0/(epm(i+1)-epm(i))/F(i+1)
      else
        eps1=ep(i-1);
        if (i < n .and. F(min(i+1,n)) > 0d0) then
          gi = -log(F(i+1)/F(i-1))/(epm(i+1)-epm(i-1));
          gi_(1) = 1d0/(epm(i+1)-epm(i-1))/F(i-1)
          gi_(3) = -1d0/(epm(i+1)-epm(i-1))/F(i+1)
        else
          gi = -log(F(i)/F(i-1))/(epm(i)-epm(i-1));
          gi_(2) = -1d0/(epm(i)-epm(i-1))/F(i)
          gi_(1) = 1d0/(epm(i)-epm(i-1))/F(i-1)
        end if
      end if
      gi = min(abs(gi),7d0)
      agi = abs(gi);
      eps2=ep(i);
      rooteps1gi=sqrt(abs(eps1*gi));  !gi.0
      rooteps2gi=sqrt(abs(eps2*gi));
      erfrooteps1gi = erf(rooteps1gi)
      erfrooteps2gi = erf(rooteps2gi)
      pik =  (2d0*Sqrt(gi)*(exp(eps2*gi)*Sqrt(eps1)*(3d0 + 2d0*eps1*gi) - &
           exp(eps1*gi)*Sqrt(eps2)*(3d0 + 2d0*eps2*gi)) + &
           3d0*exp((eps1 + eps2)*gi)*Sqrt(M_Pi)*(-erfrooteps1gi + erfrooteps2gi))/&
           (4d0*exp((eps1 + eps2 - epsI)*gi)*gi**2.5d0)

      pik_ = (-2d0*exp(eps2*gi)*(-4d0*eps1**1.5d0*epsI*gi**2.5d0 + 15d0*Sqrt(eps1*gi) + &
           10d0*(eps1*gi)**1.5d0 + 4d0*(eps1*gi)**2.5d0 - 6d0*epsI*Sqrt(eps1*gi**3d0)) + &
           2d0*exp(eps1*gi)*(-4d0*eps2**1.5d0*epsI*gi**2.5d0 + 15d0*Sqrt(eps2*gi) + &
           10d0*(eps2*gi)**1.5d0 + 4d0*(eps2*gi)**2.5d0 - 6d0*epsI*Sqrt(eps2*gi**3d0)) +& 
           3d0*exp((eps1 + eps2)*gi)*(-5d0 + 2d0*epsI*gi)*Sqrt(m_Pi)*&
           (-erfrooteps1gi + erfrooteps2gi ))/&
           (8d0*exp((eps1 + eps2 - epsI)*gi)*gi**3.5d0)
      pik_ = pik_* factor
      Te_(i) = Te_(i) + pik
      do m = 1,3
        if(gi_(m) /= 0d0) then
          Te_(i_(m)) = Te_(i_(m)) + pik_*gi_(m)
        endif
      enddo
    end do
!Electron energy in K
    Te_=Te_*2d0/3d0*m_ev2K;
  end subroutine getTeDer


!!!subroutine needed to calculated eq.34,35,36, Hagelaar article
  subroutine EvaluateA1A2A3
    Implicit None
    integer :: i,m,i_(3)    !, nval,k
    real(8) :: summedValue(3), factor, epsI, eps1, eps2, gi,agi, rooteps1gi, rooteps2gi, pik(3),&
         erfrooteps1gi, erfrooteps2gi, gi_(3),pik_(3) 
    real(8),parameter :: smallNo=1d-91

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! verification of A1, A2, A3 integrals !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!    real(8), pointer :: bolsig(:,:)
!!!!    character(len=80) :: fileF0

!!!!CAMBIO    
!!!!    fileF0 = '../input/F0_500.dat'
!!!!    open(123,file=trim(fileF0))
!!!!      read(123,*,ERR=10,END=20) nval
!!!!      allocate(bolsig(nval,2))
!!!!      do k = 1,nval
!!!!       read(123,*,ERR=10,END=20) bolsig(k,1:2)
!!!!    print*,F(k)
!!!!      end do
!!!!   pause
!!!!10 CONTINUE
!!!!20 CONTINUE
!!!!    close(123)

!!!!open(3999,file='bolsigvsluca',action='write',status='unknown')
!!!!open(4999,file='A1A2A3.dat',action='write',status='unknown')
    summedValue=0d0;
    A1_=0d0;A2_=0d0;A3_=0d0;
    do i = 1,n       !nval-1
!!!!      factor = bolsig(i,2);
      factor = F(i);
!!!!print*,'factor e F(i)',factor,bolsig(i,2),F(i)
!!!!write(3999,*)bolsig(i,1),bolsig(i,2),epm(i),F(i) 
      if(factor < smallNo) CYCLE
          epsI = epm(i);
!!!!      epsI = bolsig(i,1);
!!!!print*,'epsI',epsI      
      if(F(i) <= 0d0 .or. F(min(i+1,n)) <= 0d0 .or. F(max(i-1,1)) <= 0d0) cycle
!!!!      if(bolsig(i,2) <= 0d0 .or. bolsig(min(i+1,n),2) <= 0d0 .or. bolsig(max(i-1,1),2) <= 0d0) cycle
      gi_ = 0d0
      i_=[i-1,i,i+1]
      if (i == 1) then
        eps1=0d0;
!!!!        gi = -log(bolsig(i+1,2)/bolsig(i,2))/(bolsig(i+1,1)-bolsig(i,1));
        gi = -log(F(i+1)/F(i))/(epm(i+1)-epm(i));
!!!!        gi_(2) = 1d0/(bolsig(i+1,1)-bolsig(i,1))/bolsig(i,2)
        gi_(2) = 1d0/(epm(i+1)-epm(i))/F(i)
!!!!        gi_(3) = -1d0/(bolsig(i+1,1)-bolsig(i,1))/bolsig(i+1,2)
        gi_(3) = -1d0/(epm(i+1)-epm(i))/F(i+1)
!!!!print*,'gi- i=1',gi
      else
!!!!        eps1=(bolsig(i,1)+bolsig(i-1,1))/2;
        eps1=ep(i-1);
        if (i < n .and. F(min(i+1,n)) > 0d0) then
!!!!        if (i < n .and. bolsig(min(i+1,n),2) > 0d0) then
          gi = -log(F(i+1)/F(i-1))/(epm(i+1)-epm(i-1));
!!!!          gi = -log(bolsig(i+1,2)/bolsig(i-1,2))/(bolsig(i+1,1)-bolsig(i-1,1));
          gi_(1) = 1d0/(epm(i+1)-epm(i-1))/F(i-1)
!!!!          gi_(1) = 1d0/(bolsig(i+1,1)-bolsig(i-1,1))/bolsig(i-1,2)
          gi_(3) = -1d0/(epm(i+1)-epm(i-1))/F(i+1)
!!!!          gi_(3) = -1d0/(bolsig(i+1,1)-bolsig(i-1,1))/bolsig(i+1,2)
!!!!print*,'gi i<n e f(min(i+1,n)<0',gi
!!!!print*,'logluca logbolsig',log(F(i+1)/F(i-1)),log(bolsig(i+1,2)/bolsig(i-1,2))
!!!!print*,'expluca expbolsig',1/(epm(i+1)-epm(i-1)),1/(bolsig(i+1,1)-bolsig(i-1,1))
!!!!write(3999,*)'logaritmi ed esponenziali',log(F(i+1)/F(i-1)),log(bolsig(i+1,2)/bolsig(i-1,2)),&
!!!!1/(epm(i+1)-epm(i-1)),1/(bolsig(i+1,1)-bolsig(i-1,1))
        else
          gi = -log(F(i)/F(i-1))/(epm(i)-epm(i-1));
!!!!          gi = -log(bolsig(i,2)/bolsig(i-1,2))/(bolsig(i,1)-bolsig(i-1,1));
          gi_(2) = -1d0/(epm(i)-epm(i-1))/F(i)
!!!!          gi_(2) = -1d0/(bolsig(i,1)-bolsig(i-1,1))/bolsig(i,2)
          gi_(1) = 1d0/(epm(i)-epm(i-1))/F(i-1)
!!!!          gi_(1) = 1d0/(bolsig(i,1)-bolsig(i-1,1))/bolsig(i-1,2)
!!!!print*,'gi else',gi
        end if
      end if
      gi = min(abs(gi),7d0)
!!!!print*,'gi',gi
      agi = abs(gi);
!!!!      eps2=(bolsig(i+1,1)+bolsig(i,1))/2;
      eps2=ep(i);
!!!!print*,'eps1 e eps2',eps1,eps2!,ep(i-1),ep(i)
      rooteps1gi=sqrt(abs(eps1*gi));  !gi>0
      rooteps2gi=sqrt(abs(eps2*gi));
      erfrooteps1gi = erf(rooteps1gi)
      erfrooteps2gi = erf(rooteps2gi)
      pik(1) = (2d0*exp((-eps1 + epsI)*gi)*rooteps1gi - &
           2d0*exp((-eps2 + epsI)*gi)*rooteps2gi - &
           exp(epsI*gi)*Sqrt(m_Pi)*erfrooteps1gi + &
           exp(epsI*gi)*Sqrt(m_Pi)*erfrooteps2gi)/(2d0*gi**1.5d0)
!!!!print*,'pik1',pik(1)
      pik(2) = (2d0*Sqrt(gi)*(exp(eps2*gi)*Sqrt(eps1)*(3d0 + 2d0*eps1*gi) - &
           exp(eps1*gi)*Sqrt(eps2)*(3d0 + 2d0*eps2*gi)) + &
           3d0*exp((eps1 + eps2)*gi)*Sqrt(M_Pi)*(-erfrooteps1gi + erfrooteps2gi))/&
           (4d0*exp((eps1 + eps2 - epsI)*gi)*gi**2.5d0)
!!!!print*,'pik2',pik(2)
      pik(3) = (exp((-eps1 + epsI)*gi) - exp((-eps2 + epsI)*gi))/gi
!!!!print*,'pik3',pik(3)
      pik_(1) = -(2d0*(exp(eps2*gi)* &
           (3d0*Sqrt(eps1*gi) + 2d0*(eps1*gi)**1.5d0 - 2d0*epsI*Sqrt(eps1*gi**3d0)) + &
           exp(eps1*gi)*(-3d0*Sqrt(eps2*gi) - 2d0*(eps2*gi)**1.5d0 + &
           2d0*epsI*Sqrt(eps2*gi**3d0))) +  &
           exp((eps1 + eps2)*gi)*(-3d0 + 2d0*epsI*gi)*Sqrt(M_Pi)* &
           (erfrooteps1gi -erfrooteps2gi ))/&
           (4d0*exp((eps1 + eps2 - epsI)*gi)*gi**2.5d0)

      pik_(2) = (-2d0*exp(eps2*gi)*(-4d0*eps1**1.5d0*epsI*gi**2.5d0 + 15d0*Sqrt(eps1*gi) + &
           10d0*(eps1*gi)**1.5d0 + 4d0*(eps1*gi)**2.5d0 - 6d0*epsI*Sqrt(eps1*gi**3d0)) + &
           2d0*exp(eps1*gi)*(-4d0*eps2**1.5d0*epsI*gi**2.5d0 + 15d0*Sqrt(eps2*gi) + &
           10d0*(eps2*gi)**1.5d0 + 4d0*(eps2*gi)**2.5d0 - 6d0*epsI*Sqrt(eps2*gi**3d0)) +& 
           3d0*exp((eps1 + eps2)*gi)*(-5d0 + 2d0*epsI*gi)*Sqrt(m_Pi)*&
           (-erfrooteps1gi + erfrooteps2gi ))/&
           (8d0*exp((eps1 + eps2 - epsI)*gi)*gi**3.5d0)

      pik_(3) = (exp(eps1*gi)*(1d0 + eps2*gi - epsI*gi) + exp(eps2*gi)*(-1d0 - eps1*gi + epsI*gi))/ &
           (exp((eps1 + eps2 - epsI)*gi)*gi**2d0)
      pik_ = pik_*factor

      summedValue=summedValue+factor*pik;

      A1(i) = summedValue(1)
!!!!print*,'A1',A1(i)
      A2(i) = summedValue(2)
!!!!print*,'A2',A2(i)
      A3(i) = summedValue(3)
!!!!print*,'A3',A3(i)
!!!!write(4999,*) bolsig(i,1),A1(i),A2(i),A3(i)
      A1_(i:,i) = A1_(i:,i)+pik(1)
      A2_(i:,i) = A2_(i:,i)+pik(2)
      A3_(i:,i) = A3_(i:,i)+pik(3)
      do m = 1,3
        if(gi_(m) /= 0d0) then
          A1_(i:,i_(m)) = A1_(i:,i_(m)) + pik_(1)*gi_(m)
          A2_(i:,i_(m)) = A2_(i:,i_(m)) + pik_(2)*gi_(m)
          A3_(i:,i_(m)) = A3_(i:,i_(m)) + pik_(3)*gi_(m)
        endif
      enddo
    end do

!!!!    A3 = A3(nval-1) - A3
    A3 = A3(n) - A3
!!!!write(4999,*)'a3 aggiornato',A3
    A3(n) = 0d0
    do i = 1,n-1
      A3_(i,:) = A3_(n,:) - A3_(i,:)
    enddo
    A3_(n,:) = 0d0

!!!!close(3999)
!!!!close(4999)
!!!!pause
  end subroutine EvaluateA1A2A3

  function singleRate(val,loss, iorderIn) result(tempRate)
    Implicit None
    real(8), pointer :: val(:,:)
    integer, intent(IN), optional :: iorderIn
    integer :: i,iL, iH, ii, iorder, endV, iH0
    real(8) :: loss, oldVal
    real(8) :: tempRate, factor,pik, epsI, gi, eps1, eps2, epsL, epsH, sigmaL, sigmaH
    real(8),parameter :: smallNo=1d-80

    iorder=1
    if(present(iorderIn)) then
      iorder=iorderIn
    end if
    endV = ubound(val,1)

    tempRate=0d0;
    do i = 1,n
      factor = m_GAMMA*F(i);
      epsI = epm(i);

      if(factor < smallNo) CYCLE

      if (i == 1) then
        eps1=0d0;
        if(F(i+1) < smallNo) CYCLE
        gi = -log(F(i+1)/F(i))/(epm(i+1)-epm(i));
      else
        eps1=ep(i-1);
        if (i < n .and. F(min(i+1,n)) > smallNo .and. F(i-1) > smallNo) then
          gi = -log(F(i+1)/F(i-1))/(epm(i+1)-epm(i-1));
        elseif(F(i-1) > smallNo) then
          gi = -log(F(i)/F(i-1))/(epm(i)-epm(i-1));
        elseif(i < n .and. F(min(i+1,n)) > smallNo ) then
          gi = -log(F(i+1)/F(i))/(epm(i+1)-epm(i));
        else
          CYCLE
        end if
      end if
      eps2=ep(i);
      if (eps2 < loss) cycle;
!
      iL=-1;iH=-1;
      do ii = 1,ubound(val,1)
        if(val(ii,1) < eps1) iL=ii
        if(val(ii,1) > eps2) then
          iH = ii;exit
        end if
      enddo

      if (iL < 0 .and. epsI > val(1,1)) iL=1
      if (iH<0 .and.   epsI <= val(endV,1)) iH=size(val(:,1))
      if (iL < 0 .or. iH < 0 .or. iL == iH) cycle


      pik=0d0;
      iH0 = iH
      if (iH > iL+1) then
        do while (iL+1 < iH0) 
          iH = iL+1;
          epsH=val(iH,1);
          sigmaH=val(iH,2);
          epsL=val(iL,1);
          sigmaL=val(iL,2);
          oldval = eps2;
          eps2=epsH;
          if (eps1 < loss .and. eps2 > loss) then
            sigmaL = 0d0;
            epsL = loss
            eps1=loss;
          elseif (epsL < loss .and. epsH > loss) then
            sigmaL = 0d0;
            epsL = loss
          end if
          if( iorder == 1) then
            if(gi > m_giLim) then
              pik=pik+(exp((-eps2 + epsI)*gi)*((-2d0 + gi*(epsL + eps2*(-2d0 - eps2*gi + epsL*gi)))* &
                   sigmaH + (2d0 + gi*(2d0*eps2 - epsH + eps2*(eps2 - epsH)*gi))*sigmaL) +   &
                   exp((-eps1 + epsI)*gi)*((2d0 + gi*(2d0*eps1 - epsL + eps1*(eps1 - epsL)*gi))* &
                   sigmaH + (-2d0 + gi*(epsH + eps1*(-2d0 - eps1*gi + epsH*gi)))*sigmaL))/  &
                   ((epsH - epsL)*gi**3d0)
            else
              pik=pik+(2d0*eps1**3d0*(-sigmaH + sigmaL) + 3d0*eps1**2d0*(epsL*sigmaH - epsH*sigmaL) + &
                   eps2**2d0*(-3d0*epsL*sigmaH + 2d0*eps2*(sigmaH - sigmaL) + 3d0*epsH*sigmaL))/&
                   (6d0*(epsH - epsL))
            endif

          else
            if(gi > m_giLim) then
              pik=pik+(-(exp((-eps1 + epsI)*gi)*((-6d0 + epsL*gi*(2d0 + eps1*gi*(2d0 + eps1*gi)) - &
                   eps1*gi*(6d0 + eps1*gi*(3d0 + eps1*gi)))*sigmaH + &
                   (6d0 - epsH*gi*(2d0 + eps1*gi*(2d0 + eps1*gi)) + &
                   eps1*gi*(6d0 + eps1*gi*(3d0 + eps1*gi)))*sigmaL)) + &
                   exp((-eps2 + epsI)*gi)*((-6d0 + epsL*gi*(2d0 + eps2*gi*(2d0 + eps2*gi)) - &
                   eps2*gi*(6d0 + eps2*gi*(3d0 + eps2*gi)))*sigmaH + &
                   (6d0 - epsH*gi*(2d0 + eps2*gi*(2d0 + eps2*gi)) + &
                   eps2*gi*(6d0 + eps2*gi*(3d0 + eps2*gi)))*sigmaL))/((epsH - epsL)*gi**4d0)
            else
              pik=pik+(3d0*eps1**4d0*(-sigmaH + sigmaL) + 4d0*eps1**3d0*(epsL*sigmaH - epsH*sigmaL) + &
                   eps2**3d0*(-4d0*epsL*sigmaH + 3d0*eps2*(sigmaH - sigmaL) + 4d0*epsH*sigmaL))/  &
                   (12d0*(epsH - epsL))
            endif
          endif
          iH = iH+1;
          iL = iH-1;
          epsH=val(iH,1);
          sigmaH=val(iH,2);
          epsL=val(iL,1);
          sigmaL=val(iL,2);
          eps1=epsL;
          eps2=epsH;
        enddo
        eps2=oldval;
      else
        epsH=val(iH,1);
        sigmaH=val(iH,2);
        if (iL > 0) then
          epsL=val(iL,1);
          sigmaL=val(iL,2);
        end if
      end if
      if (eps1 < loss .and. eps2 > loss) then
        eps1=loss;
        sigmaL = 0d0;
        epsL = loss;
      elseif (epsL < loss .and. epsH > loss) then
        sigmaL = 0d0;
        epsL = loss
      end if
!
      if( iorder == 1) then

        if(gi > m_giLim) then
          pik=pik+(exp((-eps2 + epsI)*gi)*((-2d0 + gi*(epsL + eps2*(-2d0 - eps2*gi + epsL*gi)))* &
               sigmaH + (2d0 + gi*(2d0*eps2 - epsH + eps2*(eps2 - epsH)*gi))*sigmaL) +   &
               exp((-eps1 + epsI)*gi)*((2d0 + gi*(2d0*eps1 - epsL + eps1*(eps1 - epsL)*gi))* &
               sigmaH + (-2d0 + gi*(epsH + eps1*(-2d0 - eps1*gi + epsH*gi)))*sigmaL))/  &
               ((epsH - epsL)*gi**3d0)
        else
          pik=pik+(2d0*eps1**3d0*(-sigmaH + sigmaL) + 3d0*eps1**2d0*(epsL*sigmaH - epsH*sigmaL) + &
               eps2**2d0*(-3d0*epsL*sigmaH + 2d0*eps2*(sigmaH - sigmaL) + 3d0*epsH*sigmaL))/&
               (6d0*(epsH - epsL))
        endif
      else
        if(gi > m_giLim) then
          pik=pik+(-(exp((-eps1 + epsI)*gi)*((-6d0 + epsL*gi*(2d0 + eps1*gi*(2d0 + eps1*gi)) - &
               eps1*gi*(6d0 + eps1*gi*(3d0 + eps1*gi)))*sigmaH + &
               (6d0 - epsH*gi*(2d0 + eps1*gi*(2d0 + eps1*gi)) + &
               eps1*gi*(6d0 + eps1*gi*(3d0 + eps1*gi)))*sigmaL)) + &
               exp((-eps2 + epsI)*gi)*((-6d0 + epsL*gi*(2d0 + eps2*gi*(2d0 + eps2*gi)) - &
               eps2*gi*(6d0 + eps2*gi*(3d0 + eps2*gi)))*sigmaH + &
               (6d0 - epsH*gi*(2d0 + eps2*gi*(2d0 + eps2*gi)) + &
               eps2*gi*(6d0 + eps2*gi*(3d0 + eps2*gi)))*sigmaL))/((epsH - epsL)*gi**4d0)
        else
          pik=pik+(3d0*eps1**4d0*(-sigmaH + sigmaL) + 4d0*eps1**3d0*(epsL*sigmaH - epsH*sigmaL) + &
               eps2**3d0*(-4d0*epsL*sigmaH + 3d0*eps2*(sigmaH - sigmaL) + 4d0*epsH*sigmaL))/  &
               (12d0*(epsH - epsL))
        endif
      endif

      tempRate=tempRate+factor*pik
    end do

  end function singleRate

  subroutine getrate(rate)
    Implicit None
    type(t_rate),  pointer, dimension(:)  ::rate
    integer :: kk,kkk,L
    real(8) :: tempRate


    L=0
    kkloop: do kk=1,size(sigmaIN)
      do kkk=1,size(sigmaIN(kk)%sigma)
        L = L+1
        rate(L)%Rloss=sigmaIN(kk)%sigma(kkk)%loss
        rate(L)%Rmolecule=sigmaIN(kk)%sigma(kkk)%molecule
        rate(L)%Rtype=sigmaIN(kk)%sigma(kkk)%type

        tempRate = singleRate(sigmaIN(kk)%sigma(kkk)%val, rate(L)%Rloss);
        rate(L)%Rrate=tempRate
        rate(L)%indexX = sigmaIN(kk)%sigma(kkk)%indexX

      end do
    end do kkloop
  end subroutine getrate

  subroutine printVals(typeIN, moleculeIN,rateOUT)
    implicit none
    character(LEN=*), INTENT(IN) :: typeIN, moleculeIN
    real(8), INTENT(OUT) :: rateOUT
    integer :: k
    call getrate(rates)

!get the rates used in the ODE solution
    ratesLoop: do k=1, size(rates)
      if (trim(rates(k)%Rmolecule) == trim(moleculeIN)) then
        if (trim(rates(k)%Rtype) == trim(typeIN)) then
          print*,'Rate: ',rates(k)%Rrate
          rateOUT = rates(k)%Rrate
        endif
      end if
    end do ratesLoop

  end subroutine printVals

  subroutine seval(u, splineIn)
    !integer n
!
!  this subroutine evaluates the cubic spline function
!
!    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
!
!    where  x(i) .lt. u .lt. x(i+1), using horner's rule
!
!  if  u .lt. x(1) then  i = 1  is used.
!  if  u .ge. x(n) then  i = n  is used.
!
!  input..
!
!    n = the number of data points
!    u = the abscissa at which the spline is to be evaluated
!    x,y = the arrays of data abscissas and ordinates
!    b,c,d = arrays of spline coefficients computed by spline
!
!  if  u  is not in the same interval as the previous call, then a
!  binary search is performed to determine the proper interval.
!
    integer i, j, k,nxb
    TYPE(t_spline), intent(INOUT) :: splineIn
    REAL(RFREAL) dx,u
    REAL(RFREAL), POINTER :: xb(:), yi(:), coeffs(:,:,:)

    nxb=splineIn%nxb
    XB => splineIn%xb
    i = splineIn%ib
    yI => splineIn%yI
    coeffs => splineIn%coeffs


    if ( i .ge. nxb ) i = 1
    if ( u .lt. xb(i) ) go to 10
    if ( u .le. xb(i+1) ) go to 30

    if ( u .lt. xb(1) ) then
      i=1;
      go to 30
    endif
    if ( u .ge. xb(nxb) ) then
      i=nxb;
      go to 30
    endif
!
!  binary search
!
10  i = 1
    j = nxb+1
20  k = (i+j)/2
    if ( u .lt. xb(k) ) j = k
    if ( u .ge. xb(k) ) i = k
    if ( j .gt. i+1 ) go to 20
!
!  evaluate spline
!
30  dx = u - xb(i)
    splineIn%ib = i
    if(.not. splineIn%extrapolate) then
      if (i == 1) dx = max(dx,0d0)
      if (i == nxb) dx = min(dx,0d0)
    endif
    do k=1,splineIn%nvar
      yI(k) = coeffs(i,4,k)
      if(i < nxb) then
!yI(k) = coeffs(i,4,k) + dx*(coeffs(i,3,k) + dx*(coeffs(i,2,k) + dx*coeffs(i,1,k)))
        yI(k) = yI(k) + dx*(coeffs(i+1,4,k) - coeffs(i,4,k))/(xb(i+1) - xb(i))
      endif
    enddo
    return
  end subroutine seval

  subroutine resetGrid(Ymax)
    Use Solver, only: spline
    real(8) :: Ymax
    integer :: i, k

    k=0
    if(maxval(abs(F)) > 1d-30) then
      Fvals%xb = epm
      Fvals%coeffs(:,4,1) = log(max(F,1d-300))
      k=1
      call spline(Fvals%nxb, Fvals%xb, Fvals%coeffs(:,4,k), Fvals%coeffs(:,3,k) , Fvals%coeffs(:,2,k), Fvals%coeffs(:,1,k))
    endif
    !call processCrossSections(m_X,m_T, m_ne, ubound(F,1) ,Ymax,m_newtIter,1000,sigmaM)
    if(k > 0) then
      do i = 1, ubound(F,1)
        call seval(epm(i), Fvals)
        F(i) = exp(FVals%yI(1))
      enddo
    endif


  end subroutine resetGrid


  Subroutine GeometricStretchFactor(x1,xn,n,dx,guess,safety_fac,alpha)
    Implicit None

! ... input data
    Real(KIND=8) :: x1,xn,dx
    Integer :: n
    Real(KIND=8) :: guess, safety_fac

! ... output data
    Real(KIND=8) :: alpha

! ... local variables
   ! Integer :: i
    Real(KIND=8) :: alpha_old, alpha_new, f_old, f_new, err, slope
    Real(KIND=8), Parameter :: err_tol = 1D-12

! ... use secant method to solve the equation dx = (alpha - 1)/(alpha^{n-1}-1)
    alpha_old = guess
    alpha_new = 1.01_8 * guess

    f_old = (xn-x1) - (alpha_old**(n-1) - 1.0_8) / (alpha_old-1.0_8) * dx
    f_new = (xn-x1) - (alpha_new**(n-1) - 1.0_8) / (alpha_new-1.0_8) * dx
    err = dabs(f_new)

! ... error check that we got it right by guess
    If (err < err_tol) Then
      alpha = alpha_new
      Return
    End If

! ... iterate
    Do While (err > err_tol)

! ... secant line slope
      If (dabs(alpha_new - alpha_old) < err_tol) Then
        alpha = alpha_new
        Return
      End If
      slope = (f_new-f_old)/(alpha_new-alpha_old)

! ... save
      alpha_old = alpha_new
      f_old = f_new

! ... estimate new
      If (dabs(slope) >= err_tol) Then
        alpha_new = alpha_new - safety_fac * f_new / slope
      Else
        alpha = alpha_new
        Return
      End If

! ... evaluate new
      f_new = (xn-x1) - (alpha_new**(n-1) - 1.0_8) / (alpha_new-1.0_8) * dx
      err = dabs(f_new)

    End Do

    alpha = alpha_new
    Return

  End Subroutine GeometricStretchFactor


End Module ModBoltzmann


program testBoltz

  Use ModBoltzmann
  Implicit None

  real(8) :: esp(1)
  real(8) :: EbyN(1)
  real(8) :: Temperature, Composition(m_speciesCoupled)
  real(8) :: Meanenergy,mean_energy_K,rateEXC,rateION
  real(8) :: electronDensity, degreeIonization(1)
  integer :: i,j,k,m,l
  character(120) :: rate_file_ion, rate_file_exc
  character(120) :: mu_e_f_file, mu_sm_file, energy_file
!-----

  Call readCrossSections

!plascomCM inputs
  Composition = [1d0]
!  EbyN(1)=0.1d0
!  degreeIonization(1)=0.33113112148259077d0
! degreeIonization(2)=0.57543993733715659d0
! degreeIonization(3)=1d0 
  esp(1)=-6d0
!  do l=1,10
!     EbyN(l)=0d0 + 0.1d0*l
!  end do

!  do k=1,50
k=2
     EbyN(1) = 0.d0+10.d0*k;
     !EbyN(k+10) = 0.d0+10.d0*k;
!     esp(k+1) = esp(1) + 0.12d0*k;
!  enddo
!  do i=1,51
i=1
     degreeIonization(i) = 10d0**(esp(i))
!  enddo
!print*,degreeIonization(1)
mu_e_f_file = 'check_fit.dat'
mu_sm_file = 'check_fit_sM.dat' 
energy_file = 'energy_file_2.dat'
rate_file_ion = 'ratecoeff_ion.dat'
rate_file_exc = 'ratecoeff_exc.dat'

  open(2000,file = trim(energy_file),action='write',status='unknown')
  open(3000,file = trim(mu_e_f_file),action='write',status='unknown')
  open(4000,file = trim(mu_sm_file),action='write',status='unknown')
  
  open(5001,file = trim(rate_file_ion),action='write',status='unknown')
  open(5002,file = trim(rate_file_exc),action='write',status='unknown')

  write(2000,*) size(EbyN)
  write(2000,*) size(degreeIonization)
  write(2000,*) EbyN
  write(2000,*) degreeIonization
  write(3000,*) size(EbyN)
  write(3000,*) size(degreeIonization)
  write(3000,*) EbyN
  write(3000,*) degreeIonization
  write(4000,*) size(EbyN)
  write(4000,*) size(degreeIonization)
  write(4000,*) EbyN
  write(4000,*) degreeIonization

  write(5001,*) size(EbyN)
  write(5001,*) size(degreeIonization)
  write(5001,*) EbyN
  write(5001,*) degreeIonization
  write(5002,*) size(EbyN)
  write(5002,*) size(degreeIonization)
  write(5002,*) EbyN
  write(5002,*) degreeIonization
 
  Temperature = 3d2

  do m=1,size(degreeIonization)
     electronDensity = degreeIonization(m) * m_P/m_kB/Temperature

     write(2000,*) m
     write(3000,*) m, 'Ionization degree:',degreeIonization(m)
     write(4000,*) m, 'Ionization degree:',degreeIonization(m)

     write(5001,*) m
     write(5002,*) m

     do j = 1,size(EbyN) 

        Call processCrossSections(Composition, Temperature, electronDensity, 500, 50d0, 2500, TinitIn=100000d0)
        print*,'Solving E/N=',EbyN(j),'ion_degree',degreeIonization(m)
        CALL solveBoltzmann(EbyN(j),.true.,Meanenergy,mean_energy_K,sigmaM)
        !CALL solveBoltzmann(EbyN,.true.,Meanenergy)

        write(3000,*) j, 'Electric field:',EbyN(j)
        write(4000,*) j, 'Electric field:',EbyN(j)
        do i=1,n
           write(3000,*) ep(i),epm(i),F(i),Meanenergy,mean_energy_K
        enddo
        write(4000,*) sigmaM
  
        Call printVals('EXCITATION', 'N2',rateEXC)
        Call printVals('IONIZATION', 'N2',rateION)

        write(2000,*) ,Meanenergy,rateION
        write(5001,*) rateION
        write(5002,*) rateEXC

     end do
  end do

  close(2000)
  close(3000)
  close(4000)
  close(5001)
  close(5002)

  stop
end program testBoltz


