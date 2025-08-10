!####################################################################
!
!   Computation of finite field (FF) central derivatives from a
!   general electric-field dependent property.
!   Romberg-Rutishauser treatmer for the higher-order derivatives.
!
!   Code adapted from Pau Besalú-Sala (2020) to compute both high-order
!   derivatives and isotropic non-linear optical properties from
!   analytic electronic energies, dipole moments and polarizabilities.
!   It is also generalized in order to consider non-power of two steps
! 
!       - Coded by Guillem Pey, 2025   
!
!####################################################################
Module GrebmoR
implicit none
integer :: inlop,onlop,positiveFields,negativeFields,totalFields,mF,pF,nF,tF,derivative_order

Contains
    Subroutine PrintP(TensorP,numberComponents,totalFields,level)
    implicit none
    double precision, allocatable, intent(inout), dimension(:,:,:) :: TensorP
    integer, intent(in) :: level,numberComponents,totalFields
    integer :: i,j

        !-- Print the tensor: ({field strength},{component of the tensor},{X/Y/Z})
    write(*,*)
    if (level.eq.1) write(*,*) "    Printing the components of the tensor for the X-derivative"
    if (level.eq.2) write(*,*) "    Printing the components of the tensor for the Y-derivative"
    if (level.eq.3) write(*,*) "    Printing the components of the tensor for the Z-derivative"
    do i=1,totalFields
        write(*,*) (TensorP(i,j,level),j=1,numberComponents)
    end do

    End subroutine PrintP

    Subroutine ReadValues(inlop,line,size,isEnergy,fchk1,fchk2,substr1,substr2,direction,strength,lineP)
    implicit none
    logical,intent(in) :: isEnergy
    integer,intent(in) :: size,inlop
    character,intent(in) :: substr1,substr2
    character*1 :: substr11,substr22
    character*6,intent(in) :: fchk1,fchk2
    character*200,intent(in) :: line
    character*1,intent(out) :: direction
    character*200,intent(out) :: lineP
    double precision,intent(out) :: strength
    integer :: Findex,indexF,ios

        !-- Get the direction of the applied field, intent(out)
    Findex=index(line,substr1)
    Findex=Findex+index(line,substr2)
    read(line(Findex:Findex),*) direction

        !-- Get the field strength, intent(out)
    if (isEnergy.eqv..FALSE.) then
        indexF=index(line,fchk1)
    else
        indexF=index(line,fchk2)
    end if
    read(line(Findex+1:indexF-1),*,iostat=ios) strength

        !-- Get the substring of only properties, intent(out)
    if (inlop.ne.0) lineP=line(indexF+6:size)
    if (inlop.eq.0) lineP=line(indexF+7:size)

    End subroutine Readvalues

    Subroutine DeclareProperties(isEnergy,substr1,substr2,inlop,mF,totalFields,TensorP)
    implicit none
    logical, intent(inout) :: isEnergy
    character*200 :: line,postLine
    character*1 :: field_direction
    character, intent(in) :: substr1,substr2
    integer, intent(in) :: inlop,mF,totalFields
    integer :: ios,nF,pF,level
    double precision :: field_strength
    double precision, dimension(5,2) :: P_general
    double precision, allocatable, intent(out), dimension(:,:,:) :: TensorP

    !-- Acabar d'adaptar per pillar d'input l'energia i el moment dipolar

    if (inlop.eq.0) allocate(TensorP(totalFields,2,3))    
    if (inlop.eq.1) allocate(TensorP(totalFields,4,3))   
    if (inlop.eq.2) allocate(TensorP(totalFields,7,3))   
    if (inlop.eq.3) allocate(TensorP(totalFields,11,3))    

        !-- Get field direction for the later property parsing
    if (substr1.eq."x".or.substr2.eq."X") level=1
    if (substr1.eq."y".or.substr2.eq."Y") level=2
    if (substr1.eq."z".or.substr2.eq."Z") level=3
    
    nF=0; pF=0
    !if (inlop.eq.0) backspace(4)
    do while(nF+pF+1.ne.totalFields)

            !-- Get the necessary properties to parse them into the ComputeDerivatives subroutine
        read(4,'(A)',iostat=ios,end=1523) line
            1523 continue
            !write(*,*) line
           !call ReadValues(inlop,line,size,isEnergy,fchk1,fchk2,substr1,substr2,direction,strength,lineP)
            call ReadValues(inlop,line,len_trim(line),isEnergy,".fchk-",".fchk:",substr1,substr2,field_direction,field_strength,postLine)
            if (field_strength.lt.0.0d0) nF=nF+1
            if (field_strength.gt.0.0d0) pF=pF+1
            read(postLine,*,iostat=ios,end=1927) P_general(1,1),P_general(2,1),P_general(3,1),P_general(4,1),P_general(5,1)
        if (inlop.ne.0) read(4,'(A)',iostat=ios) line
            1927 continue
            !write(*,*) line
            if (inlop.ne.0) call ReadValues(inlop,line,len_trim(line),isEnergy,".fchk-",".fchk:",substr1,substr2,field_direction,field_strength,postLine)

                !-- Condition only applicable for alpha (only P_general(1,2)) and beta (P_general(1:5,2))
            if (inlop.ge.2) read(postLine,*,iostat=ios) P_general(1,2),P_general(2,2),P_general(3,2),P_general(4,2),P_general(5,2)
        read(4,*,iostat=ios,end=1612) line
            !write(*,*) line

            !-- Has to be specifically adapted to consider the change of sign of the beta tensor.
        1612 continue
        if (nF.gt.pF) then
            TensorP(mF-nF,1,level)=field_strength
            if (inlop.eq.0) TensorP(mF-nF,2,level)=P_general(1,1)
            if (inlop.eq.1) TensorP(mF-nF,2:4,level)=P_general(1:3,1)
            if (inlop.ge.2) TensorP(mF-nF,2:6,level)=P_general(1:5,1)
            if (inlop.eq.2) TensorP(mF-nF,7,level)=P_general(1,2)
            if (inlop.eq.3) TensorP(mF-nF,7:11,level)=-1.0d0*P_general(1:5,2)
        else if (nF.eq.pF) then
            TensorP(mF+nF,1,level)=field_strength
            if (inlop.eq.0) TensorP(mF+nF,2,level)=P_general(1,1)
            if (inlop.eq.1) TensorP(mF+nF,2:4,level)=P_general(1:3,1)
            if (inlop.ge.2) TensorP(mF+nF,2:6,level)=P_general(1:5,1)
            if (inlop.eq.2) TensorP(mF+nF,7,level)=P_general(1,2)
            if (inlop.eq.3) TensorP(mF+nF,7:11,level)=-1.0d0*P_general(1:5,2)
        end if

    end do

    End subroutine DeclareProperties

    Subroutine ComputeDerivatives(axis,component,inlop,onlop,field_direction,o_derivative,mF,totalFields,F,P,secRombergP,mainRombergP)
    implicit none
    character*1,dimension(3) :: dipoleComponents
    character*2,dimension(6) :: alphaComponents
    character*3,dimension(10) :: betaComponents
    character*1, intent(in) :: field_direction
    integer, intent(in) :: axis,component,totalFields,inlop,onlop,mF,o_derivative 
    double precision, intent(in), dimension(totalFields) :: F,P
    double precision, intent(out), dimension(3) :: mainRombergP,secRombergP
    double precision :: a,main_errRomberg,sec_errRomberg
    double precision :: energyDerivative = 1.0d0
    double precision, allocatable, dimension(:,:) :: RombergP    
    double precision, allocatable, dimension(:,:) :: errorP    
    double precision, allocatable, dimension(:) :: FF
    integer, dimension(2,2) :: minloc_errRomberg
    integer :: i,j

    !-- Reference: M. Medved et al. / Journal of Molecular Structure: THEOCHEM 847 (2007) 39–46

    if(allocated(FF)) deallocate(FF)
    allocate(FF(mF-1))
    if (inlop.eq.0) energyDerivative=-1.0d0 

        !-- String variables to parse into Romberg for printing
    dipoleComponents=(/"X","Y","Z"/)
    alphaComponents=(/"XX","XY","YY","XZ","YZ","ZZ"/)
    betaComponents=(/"XXX","XXY","XYY","YYY","XXZ","XYZ","YZZ","XZZ","YZZ","ZZZ"/)

    a=F(totalFields)/F(totalFields-1)

        !-- Change the formula for the general ones for a general step
    if (o_derivative.eq.1) then !compute finite-field's first derivative
        do i=1,mF-1
            FF(i)=energyDerivative*(P(i)-P(mF+i))/(2.0d0*F(mF+I))
        end do
    else if (o_derivative.eq.2) then !compute finite-field's second derivative
        do i=1,mF-1
            FF(i)=energyDerivative*(P(i)-2.0d0*P(mF)+P(mF+i))/F(mF+i)**2.0d0
        end do
    else if (o_derivative.eq.3) then !compute finite-field's third derivative
        do i=1,mF-1
            FF(i)=energyDerivative*(P(i)-2.0d0*P(i)+2.0d0*P(mF+i)-P(mF+i+1))/(2.0d0*F(mF+i)**3)
        end do
    else if (o_derivative.eq.4) then !compute finite-field's fourth derivative
        do i=1,mF-1
            FF(i)=energyDerivative*(P(i+2)-4.0d0*P(i)+6.0d0*P(mF)-4.0d0*P(mF+i)+P(mF+i+2))/F(mF+i)**4
        end do
    end if

    do i=1,totalFields
        write(*,*) F(i),P(i)
    end do

        ! E    ! μ    ! α    ! β    ! γ
    if (inlop.eq.0) write(*,'(" RomberG - Computing the derivative of the Energy with respect to ",A1)') field_direction
    if (inlop.eq.1) write(*,'(" RomberG - Computing the derivative of μ",A2," with respect to ",A1)') dipoleComponents(component-1),field_direction
    if (inlop.eq.2) write(*,'(" RomberG - Computing the derivative of α ",A2," with respect to ",A1)') alphaComponents(component-1),field_direction
    if (inlop.eq.3) write(*,'(" RomberG - Computing the derivative of β",A3," with respect to ",A1)') betaComponents(component-1),field_direction
    call RombergProcedure(a,derivative_order,FF,minloc_errRomberg,errorP,main_errRomberg,sec_errRomberg,RombergP)
    
        !-- Get the output values in vector form for the later handling
    mainRombergP(1)=RombergP(minloc_errRomberg(1,1),minloc_errRomberg(2,1)); mainRombergP(2)=main_errRomberg; mainRombergP(3)=1.0d2*mainRombergP(2)/mainRombergP(1)
    secRombergP(1)=RombergP(minloc_errRomberg(1,2),minloc_errRomberg(2,2)); secRombergP(2)=sec_errRomberg; secRombergP(3)=1.0d2*secRombergP(2)/secRombergP(1)

    End subroutine ComputeDerivatives

    Subroutine RombergProcedure(a,derivative_order,P,err_minloc,errRombergT,abs_errRomberg1,abs_errRomberg2,RombergT)
    implicit none
    integer, intent(in) :: derivative_order
    double precision, intent(in) :: a
    double precision, allocatable, intent(in), dimension(:) :: P
    double precision, allocatable, intent(out), dimension(:,:) :: RombergT,errRombergT
    double precision, allocatable, dimension(:,:) :: errRomberg
    integer, intent(out), dimension(2,2) :: err_minloc
    double precision :: abs_errRomberg1,abs_errRomberg2
    integer :: iterRR
    integer :: i,j

    !-- Reference: M. Medved et al. / Journal of Molecular Structure: THEOCHEM 847 (2007) 39–46

    if (allocated(RombergT)) deallocate(RombergT)
    allocate(RombergT(mF-1,mF-1))
    if (allocated(errRombergT)) deallocate(errRombergT)
    allocate(errRombergT(mF-2,mF-2))
    if (allocated(errRomberg)) deallocate(errRomberg)
    allocate(errRomberg(mF-2,mF-2))

        !-- Define Romberg matrix
    RombergT(1:mF-1,1)=P
    RombergT(1:mF-1,2:mF-1)=9.9d99
    errRomberg=9.9d99

        !-- Compute & print the Romberg triangle
    do i=2,mF-1     !-- Compute from the second column, the first is the "zero-th" Romberg iteration
        iterRR=i-1  !-- Compute the Romberg iteration
        do j=1,mF-i !-- Compute the Romberg for each row
            RombergT(j,i)=(a**(2*iterRR)*RombergT(j,i-1)-RombergT(j+1,i-1))/(a**(2*iterRR)-1.0d0)
        end do
    end do
    
    write(*,*) ("                            ROMBERG TRIANGLE",i=1,3)
    write(*,*) "Iteration:",(j-1,"              ",j=1,mF-1)
    do i=1,mF-1
        do j=1,mF-1
            if (RombergT(i,j).eq.9.9d99) then
                write(*,'(A)',advance='NO') "   "
            else
                write(*,'(xxxxxxxxF20.6)',advance="no") RombergT(i,j)
            end if
        end do
        write(*,*)
    end do

        !-- Compute the Romberg absolute error matrix (errRomberg,intent(out)) and save:
        !       # the lowest (errRomberg1) error
        !       # the second-to-lowest (errRomberg2) error 
        !       # its respective value within 'errRomberg'
    write(*,*)
    do i=1,mF-2
        do j=1,mF-2
            if (RombergT(j,i).ge.9.9d99.or.Rombergt(j+1,i).ge.9.9d99) cycle
            errRomberg(j,i)=abs(RombergT(j,i)-RombergT(j+1,i))
        end do
    end do
    errRombergT=errRomberg

    write(*,*) ("                        ROMBERG ERROR TRIANGLE",i=1,3)
    write(*,*) "Iteration:",(j-1,"              ",j=1,mF-1)
    do i=1,mF-2
        do j=1,mF-2
            if (errRombergT(i,j).eq.9.9d99) then
                write(*,'(A)',advance='NO') "   "
            else
                write(*,'(xxxxxxxxF20.6)',advance="no") errRomberg(i,j)
            end if
        end do
        write(*,*)
    end do

        !-- Get the position of the value with the lowest and second-to-lowest errors
    abs_errRomberg1=minval(errRomberg)
    err_minloc(1:2,1)=minloc(errRomberg)
        errRomberg(err_minloc(1,1),err_minloc(2,1))=9.9d99
    abs_errRomberg2=minval(errRomberg) 
    err_minloc(1:2,2)=minloc(errRomberg)
    deallocate(errRomberg)

    write(*,*)
    write(*,'(" RomberG - Minimum value:",xF20.10,xx"and second best:",xF20.10)') RombergT(err_minloc(1,1),err_minloc(2,1)),RombergT(err_minloc(1,2),err_minloc(2,2))
    write(*,'(" RomberG - Minimum absolute errors for the properties:",xF20.10,xx"and:",xF20.10)') abs_errRomberg1,abs_errRomberg2
    write(*,'(" RomberG - Romberg errors:",xF20.10,"% &&",xF20.10,"%")') 1.0d2*abs_errRomberg1/RombergT(err_minloc(1,1),err_minloc(2,1)),1.0d2*abs_errRomberg2/RombergT(err_minloc(1,2),err_minloc(2,2))
    write(*,*)

    End subroutine RombergProcedure

    Subroutine Reps(isEnergy,printProperties,indexString,inlop,onlop,components,derivative_order,mF,totalFields,TensorP,secRombergP,mainRombergP)
    implicit none
    character*1,dimension(6) :: derivative_substr
    logical, intent(inout) :: isEnergy
    logical, intent(in) :: printProperties
    character*3, intent(in), dimension(components) :: indexString
    double precision, allocatable, intent(inout), dimension(:,:,:) :: mainRombergP,secRombergP,TensorP
    integer, intent(in) :: inlop,onlop,components,derivative_order,mF,totalFields
    double precision, dimension(3) :: secP,mainP
    integer :: i,j,k,n_dim
    
    derivative_substr=(/"x","y","z","X","Y","Z"/)

    do k=1,3
            !-- Get the properties into a three-layer tensor for the respective derivative computation
                       !call DeclareProperties(isEnergy,substr1,substr2,inlop,mF,totalFields,TensorP)
        call DeclareProperties(isEnergy,derivative_substr(k),derivative_substr(k+3),inlop,mF,totalFields,TensorP)
        if (printProperties.eqv..TRUE.) call PrintP(TensorP,components+1,totalFields,k)

        do i=2,components

                !-- Compute the finite-field approach for the given property and perform the Romberg procedure to reduce the error
            call ComputeDerivatives(k,i,inlop,onlop,derivative_substr(k+3),derivative_order,mF,totalFields,TensorP(1:totalFields,1,k),TensorP(1:totalFields,i,k),secP,mainP)
           !call ComputeDerivatives(axis,component,inlop,onlop,field_direction,o_derivative,stepSQRT,mF,totalFields,F,P,secRombergTensor,mainRombergTensor)

                !-- Get the full output values of the Romberg iterations
            mainRombergP(i-1,1,k)=mainP(1); mainRombergP(i-1,2,k)=mainP(2); mainRombergP(i-1,3,k)=mainP(3)
            secRombergP(i-1,1,k)=secP(1); secRombergP(i-1,2,k)=secP(2); secRombergP(i-1,3,k)=secP(3)

        end do
        write(*,*)
        write(*,*) "######################################################"
        write(*,*) "######################################################"
        write(*,*) "######################################################"
        write(*,*)
    end do
    
        !-- Print the mainRomberG output values
    do k=1,3
        write(*,*) ("                         ",indexString(i),i=1,components-1) 
        write(*,*)
        write(*,*) derivative_substr(k+3),"             Value=",(mainRombergP(j,1,k),j=1,components-1)
        write(*,*) derivative_substr(k+3),"    Absolute error=",(mainRombergP(j,2,k),j=1,components-1)
        write(*,*) derivative_substr(k+3),"  RomberG error(%)=",(mainRombergP(j,3,k),j=1,components-1)
    end do
    if (printProperties.eqv..TRUE.) then
        do k=1,3
            write(*,*) ("                       ",indexString(i),i=1,n_dim) 
            write(*,*)
            write(*,*) derivative_substr(k+3),"             Value=",(secRombergP(j,1,k),j=1,components-1)
            write(*,*) derivative_substr(k+3),"    Absolute error=",(secRombergP(j,2,k),j=1,components-1)
            write(*,*) derivative_substr(k+3),"  RomberG error(%)=",(secRombergP(j,3,k),j=1,components-1)
        end do
    end if
    write(*,*)

    End subroutine

End module GrebmoR
!###################################################################################################
!###################################################################################################
!###################################################################################################
!###################################################################################################
Program RomberG
use f90getopt   !-- Argument parser: https://github.com/haniibrahim/f90getopt
use GrebmoR     !-- Utilities module for RomberG
implicit none
character*1 :: arg,field_direction
character*1,dimension(6) :: derivative_substr
character :: dipoleComponents(3)*1, alphaComponents(6)*2,betaComponents(10)*3
character*200 :: line,postLine
character*100 :: basename,mol_name
logical :: doLongitudinal = .FALSE.
logical :: doIsotropic = .FALSE.
logical :: doSQRTstep = .FALSE.
logical :: printProperties = .FALSE.
logical :: isEnergy,isDipole,isAlpha,isBeta,isGamma
    !-- Tensors in which the properties are gathered
double precision, allocatable, dimension (:,:,:) :: P_energy,P_dipole,P_alpha,P_beta
    !-- Output tensors for each property
double precision, dimension(3) :: mainP,secP
double precision, allocatable, dimension (:,:,:) :: mainRombergEnergy,secRombergEnergy
double precision, allocatable, dimension (:,:,:) :: mainRombergDipole,secRombergDipole
double precision, allocatable, dimension (:,:,:) :: mainRombergAlpha,secRombergAlpha
double precision, allocatable, dimension (:,:,:) :: mainRombergBeta,secRombergBeta
    !-- Dummy matrix to store the properties 
double precision, dimension (5,2) :: P_general
double precision :: field_strength
integer :: indexEnergy,indexDipole,indexAlpha,indexBeta,indexGamma
integer :: Xindex,Yindex,Zindex,indexBase,indexField
integer :: stat1,ios,err
integer :: i,j,k,dummy

type(option_s) :: opts(6)
!-- opts(i)=option_s(long_name,arguments?,short_name)
opts(1) = option_s("input",.TRUE.,"i")
opts(2) = option_s("output",.TRUE.,"o")
opts(3) = option_s("longitudinal",.FALSE.,"l")
opts(4) = option_s("isotropic",.FALSE.,"I")
opts(5) = option_s("positive-fields",.TRUE.,"F")
opts(6) = option_s("P",.FALSE.,"p")

derivative_substr=(/"x","y","z","X","Y","Z"/)
dipoleComponents=(/"X  ","Y  ","Z  "/)
alphaComponents=(/"XX ","XY ","YY ","XZ ","YZ ","ZZ "/)
betaComponents=(/"XXX","XXY","XYY","YYY","XXZ","XYZ","YZZ","XZZ","YZZ","ZZZ"/)

if (command_argument_count().le.1) then
    write(*,*) "$ man RomberG.f95"
    write(*,*) "./ROMBERG.exe {options} {name_of_the_molecule}"
    stop
end if

    !-- Get command line arguments
do
    arg=getopt("i:o:F:lIph",opts)
    select case(arg)
        case(char(0))
            exit
        case("i","input")
                !-- Catch what are the OUTPUT values
            isEnergy=.FALSE.; isDipole=.FALSE.; isAlpha=.FALSE.; isBeta=.FALSE.; isGamma=.FALSE.
            call CatchProperties(optarg,"Energy","energy","E     ",indexEnergy,isEnergy)
            call CatchProperties(optarg,"Dipole","dipole","M     ",indexDipole,isDipole)
            call CatchProperties(optarg,"Alpha ","alpha ","A     ",indexAlpha,isAlpha)
            call CatchProperties(optarg,"Beta  ","beta  ","B     ",indexBeta,isBeta)

            if (indexEnergy.eq.0.and.indexDipole.eq.0.and.indexAlpha.eq.0.and.indexBeta.eq.0) then
                write(*,*) " Not valid INPUT option. Legal options:"
                write(*,*) " for energy:    Energy, energy, E"
                write(*,*) " for dipole:    Dipole, dipole, M"
                write(*,*) " for alpha:     Alpha,  alpha,  A"
                write(*,*) " for beta:      Beta,   beta,   B"
                stop !Intentar millorar aquest error handling
            end if
            
            if (isEnergy.eqv..TRUE.) inlop=0
            if (isDipole.eqv..TRUE.) inlop=1
            if (isAlpha.eqv..TRUE.)  inlop=2
            if (isBeta.eqv..TRUE.)   inlop=3

        case("o","output")
                !-- Catch what are the OUTPUT values
            isEnergy=.FALSE.; isDipole=.FALSE.; isAlpha=.FALSE.; isBeta=.FALSE.; isGamma=.FALSE.
            call CatchProperties(optarg,"Energy","energy","E     ",indexEnergy,isEnergy)
            call CatchProperties(optarg,"Dipole","dipole","M     ",indexDipole,isDipole)
            call CatchProperties(optarg,"Alpha ","alpha ","A     ",indexAlpha,isAlpha)
            call CatchProperties(optarg,"Beta  ","beta  ","B     ",indexBeta,isBeta)
            call CatchProperties(optarg,"Gamma ","gamma ","G     ",indexGamma,isGamma)

            if (indexDipole.eq.0.and.indexAlpha.eq.0.and.indexBeta.eq.0.and.indexGamma.eq.0) then
                write(*,*) " Not valid OUTPUT option. Legal options:"
                write(*,*) " for dipole:    Dipole, dipole, M"
                write(*,*) " for alpha:     Alpha,  alpha,  A"
                write(*,*) " for beta:      Beta,   beta,   B"
                write(*,*) " for gamma:     Gamma,  gamma,  G"
                stop !Intentar millorar aquest error handling
            else if(indexEnergy.gt.0) then
                write(*,*) " Energy is not a valid output. Stop"
                stop
            end if

            if (isEnergy.eqv..TRUE.) onlop=0
            if (isDipole.eqv..TRUE.) onlop=1
            if (isAlpha.eqv..TRUE.)  onlop=2
            if (isBeta.eqv..TRUE.)   onlop=3
            if (isGamma.eqv..TRUE.)  onlop=4

        case("F","total-fields")
            read(optarg,*) totalFields
            positiveFields=(totalFields-1)/2
            negativeFields=positiveFields
            mF=positiveFields+1

        case("l","longitudinal")
            doLongitudinal=.TRUE.

        case("I","isotropic")
            doIsotropic=.TRUE.

        case("p","P")
            printProperties=.TRUE.

        case("h")
            !-- To do
            write(*,*) "Help"
    end select
end do

    !-- Compute the derivative order to be performed
derivative_order=onlop-inlop
if (derivative_order.le.0) then
    write(*,*) " ERROR! Not valid input line! Stop."
    stop
end if

    !-- Get the name of the molecule from the command execution line
call get_command_argument(command_argument_count(),value=mol_name,status=stat1)

    !-- Get the field-dependent properties from the .fchk files
if (inlop.eq.0) then
    isEnergy=.TRUE.
    call system("rm isEnergy.nlop")
    call system("cd $(pwd)/"//mol_name//" ; grep -A 0 'Total Energy' *.fchk > isEnergy.nlop; sed -i 's/Total Energy                               R/T/g' isEnergy.nlop; cp isEnergy.nlop ../")
    open (unit=4,file="isEnergy.nlop",status="old")
else if (inlop.eq.1) then
    isEnergy=.FALSE.
    call system("rm isDipole.nlop")
    call system("cd $(pwd)/"//mol_name//" ; grep -A 1 'Dipole' *.fchk > isDipole.nlop; sed -i '/Dipole/d' isDipole.nlop; cp isDipole.nlop ../")
    open (unit=4,file="isDipole.nlop",status="old")
else if (inlop.eq.2) then
    isEnergy=.FALSE.
    call system("rm isAlpha.nlop")
    call system("cd $(pwd)/"//mol_name//" ; grep -A 2 'Polarizability' *.fchk > isAlpha.nlop; sed -i '/Polarizability/d' isAlpha.nlop; cp isAlpha.nlop ../")
    open (unit=4,file="isAlpha.nlop",status="old")
else if (inlop.eq.3) then
        !-- Remember that for beta the sign is changed with respect to the regular convetion!
    isEnergy=.FALSE.
    call system("rm isBeta.nlop")
    call system("cd $(pwd)/"//mol_name//" ; grep -A 2 'Hyperpolarizability' *.fchk > isBeta.nlop; sed -i '/Hyperpolarizability/d' isBeta.nlop; cp isBeta.nlop ../")
    open (unit=4,file="isBeta.nlop",status="old")
end if

    !-- Read the zero-field name and NLOPs
                read(4,*,iostat=ios) line,P_general(1,1),P_general(2,1),P_general(3,1),P_general(4,1),P_general(5,1)
if (inlop.ge.2) read(4,*,iostat=ios) line,P_general(1,2),P_general(2,2),P_general(3,2),P_general(4,2),P_general(5,2)

    !-- Declare the zero-field properties according to Gaussian's .fchk format
        !-- He d'intentar millorar aquesta part del codi
dummy=1
if (inlop.eq.3) then
    allocate(P_beta(totalFields,11,3))
    P_beta(mF,1,1:3)=0.0d0
    do j=1,2; do i=1,5
            dummy=dummy+1
            P_beta(mF,dummy,1:3)=P_general(i,j)
    end do; end do
else if (inlop.eq.2) then
    allocate(P_alpha(totalFields,7,3))
    P_alpha(mF,1,1:3)=0.0d0
    do i=1,5
        dummy=dummy+1
        P_alpha(mF,dummy,1:3)=P_general(i,1)
    end do
    P_alpha(mF,dummy+1,1:3)=P_general(1,2)
else if (inlop.eq.1) then
    allocate(P_dipole(totalfields,4,3))
    P_dipole(mF,1,1:3)=0.0d0
    do i=1,3
        dummy=dummy+1
        P_dipole(mF,dummy,1:3)=P_general(i,1)
    end do
else if (inlop.eq.0) then
    allocate(P_energy(totalfields,2,3))
    P_energy(mF,1,1:3)=0.0d0
    P_energy(mF,2,1:3)=P_general(1,1)
end if

write(*,*)
if (inlop.eq.0.and.onlop.eq.1) write(*,*) "                     COMPUTING THE DIPOLE MOMENT FROM THE ELECTORNIC ENERGY"
if (inlop.eq.0.and.onlop.eq.2) write(*,*) "                     COMPUTING THE DIAGONAL COMPONENTS OF THE POLARIZABILITY MATRIX FROM THE ELECTORNIC ENERGY"
if (inlop.eq.0.and.onlop.eq.3) write(*,*) "                     COMPUTING THE DIAGONAL COMPONENTS OF THE FIRST HYPERPOLARIZABILITY TENSOR FROM THE ELECTORNIC ENERGY"
if (inlop.eq.0.and.onlop.eq.4) write(*,*) "                     COMPUTING THE DIAGONAL COMPONENTS OF THE SECOND HYPERPOLARIZABILITY TENSOR FROM THE ELECTORNIC ENERGY"
if (inlop.eq.1.and.onlop.eq.2) write(*,*) "                     COMPUTING THE COMPONENTS OF THE POLARIZABILITY MATRIX FROM THE DIPOLE MOMENT"
if (inlop.eq.1.and.onlop.eq.3) write(*,*) "                     COMPUTING THE COMPONENTS OF THE FIRST HYPERPOLARIZABILITY TENSOR FROM THE DIPOLE MOMENT"
if (inlop.eq.1.and.onlop.eq.4) write(*,*) "                     COMPUTING THE COMPONENTS OF THE SECOND HYPERPOLARIZABILITY TENSOR FROM THE DIPOLE MOMENT"
if (inlop.eq.2.and.onlop.eq.3) write(*,*) "                     COMPUTING THE COMPONENTS OF THE FIRST HYPERPOLARIZABILITY TENSOR FROM THE POLARIZABILITY MATRIX"
if (inlop.eq.2.and.onlop.eq.4) write(*,*) "                     COMPUTING THE COMPONENTS OF THE SECOND HYPERPOLARIZABILITY TENSOR FROM THE POLARIZABILITY MATRIX"
if (inlop.eq.3.and.onlop.eq.4) write(*,*) "                     COMPUTING THE COMPONENTS OF THE SECOND HYPERPOLARIZABILITY TENSOR FROM THE FIRST HYPERPOLARIZABILITY TENSOR"
write(*,*) "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
write(*,*)

    !-- (components+1) because the Romberg data is stored from the second column onwards
if (inlop.eq.2) then !-- Compute the derivatives of alpha

    if (allocated(mainRombergAlpha)) deallocate(mainRombergAlpha)
    if (allocated(secRombergAlpha)) deallocate(secRombergAlpha)
    allocate(mainRombergAlpha(6,3,3)); allocate(secRombergAlpha(6,3,3))
   !call Reps(isEnergy,printProperties,indexString,inlop,onlop,components,derivative_order,mF,totalFields,TensorP,secRombergP,mainRombergP)
    call Reps(isEnergy,printProperties,alphaComponents,inlop,onlop,7,derivative_order,mF,totalFields,P_alpha,secRombergAlpha,mainRombergAlpha)

else if (inlop.eq.1) then !-- Compute the derivatives of the dipole moment

        !-- This section has to be revised
    if (allocated(mainRombergDipole)) deallocate(mainRombergDipole)
    if (allocated(secRombergDipole)) deallocate(secRombergDipole)
    allocate (mainRombergDipole(3,3,3)); allocate(secRombergDipole(3,3,3))
   !call Reps(isEnergy,printProperties,indexString,inlop,onlop,components,derivative_order,mF,totalFields,TensorP,secRombergP,mainRombergP)
    call Reps(isEnergy,printProperties,dipoleComponents,inlop,onlop,4,derivative_order,mF,totalFields,P_dipole,secRombergDipole,mainRombergDipole)

else if (inlop.eq.0) then !-- Compute the derivatives of the energy

    if (allocated(mainRombergEnergy)) deallocate(mainRombergEnergy)
    if (allocated(secRombergEnergy)) deallocate(secRombergEnergy)
    allocate (mainRombergEnergy(3,3,3)); allocate(secRombergEnergy(3,3,3))
   !call Reps(isEnergy,printProperties,indexString,inlop,onlop,components,derivative_order,mF,totalFields,TensorP,secRombergP,mainRombergP)
    call Reps(isEnergy,printProperties,dipoleComponents,inlop,onlop,2,derivative_order,mF,totalFields,P_energy,secRombergEnergy,mainRombergEnergy)

end if

    !-- Computing the isotropic values
write(*,*) "Polla"
End program

Subroutine CatchProperties(optarg,substr1,substr2,substr3,indexP,isP)
implicit none
character*30,intent(in) :: optarg
character*6,intent(in) :: substr1,substr2,substr3
integer, intent(out) :: indexP
logical, intent(out) :: isP

indexP=index(optarg,trim(substr1))
indexP=indexP+index(optarg,trim(substr2))
indexP=indexP+index(optarg,trim(substr3))
if(indexP.gt.0) isP = .TRUE.

End subroutine CatchProperties
