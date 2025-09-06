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
    do while(nF+pF+1.ne.totalFields)

            !-- Get the necessary properties to parse them into the ComputeDerivatives subroutine
                !-- This should be reworked for each kind of inlop to "simplify"
        if (inlop.le.2) then
            read(4,'(A)',iostat=ios,end=1523) line
                1523 continue
               !call ReadValues(inlop,line,size,isEnergy,fchk1,fchk2,substr1,substr2,direction,strength,lineP)
                call ReadValues(inlop,line,len_trim(line),isEnergy,".fchk-",".fchk:",substr1,substr2,field_direction,field_strength,postLine)
                if (field_strength.lt.0.0d0) nF=nF+1
                if (field_strength.gt.0.0d0) pF=pF+1
                read(postLine,*,iostat=ios,end=1927) P_general(1,1),P_general(2,1),P_general(3,1),P_general(4,1),P_general(5,1)
            if (inlop.ne.0) read(4,'(A)',iostat=ios) line
                1927 continue
                if (inlop.ne.0) call ReadValues(inlop,line,len_trim(line),isEnergy,".fchk-",".fchk:",substr1,substr2,field_direction,field_strength,postLine)

                    !-- Condition only applicable for alpha (only P_general(1,2)) and beta (P_general(1:5,2))
                if (inlop.ge.2) read(postLine,*,iostat=ios) P_general(1,2),P_general(2,2),P_general(3,2),P_general(4,2),P_general(5,2)
            read(4,*,iostat=ios,end=1612) line

        else if (inlop.eq.3) then
            read(4,'(A)',iostat=ios,end=1337) line
            1337 continue
               !call ReadValues(inlop,line,size,isEnergy,fchk1,fchk2,substr1,substr2,direction,strength,lineP)
                call ReadValues(inlop,line,len_trim(line),isEnergy,".fchk-",".fchk:",substr1,substr2,field_direction,field_strength,postLine)
                if (field_strength.lt.0.0d0) nF=nF+1
                if (field_strength.gt.0.0d0) pF=pF+1
                read(postLine,*,iostat=ios,end=1339) P_general(1,1),P_general(2,1),P_general(3,1),P_general(4,1),P_general(5,1)
            read(4,'(A)',iostat=ios) line
            1339 continue
               !call ReadValues(inlop,line,size,isEnergy,fchk1,fchk2,substr1,substr2,direction,strength,lineP)
                call ReadValues(inlop,line,len_trim(line),isEnergy,".fchk-",".fchk:",substr1,substr2,field_direction,field_strength,postLine)
                read(postLine,*,iostat=ios) P_general(1,2),P_general(2,2),P_general(3,2),P_general(4,2),P_general(5,2)
            read(4,'(A)',iostat=ios,end=1332) line
            1332 continue
        end if

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
    if (o_derivative.le.2) then
        allocate(FF(mF-1))
    else if (o_derivative.eq.3) then
        allocate(FF(mF-2))
    else if (o_derivative.eq.4) then
        allocate(FF(mF-3))
    end if
    if (inlop.eq.0) energyDerivative=-1.0d0 

        !-- String variables to parse into Romberg for printing
    dipoleComponents=(/"X","Y","Z"/)
    alphaComponents=(/"XX","XY","YY","XZ","YZ","ZZ"/)
    betaComponents=(/"XXX","XXY","XYY","YYY","XXZ","XYZ","YZZ","XZZ","YZZ","ZZZ"/)

    a=F(totalFields)/F(totalFields-1)

        !-- General formulas for up to the fourth derivative
    if (o_derivative.eq.1) then !compute finite-field's first derivative
        do i=1,mF-1
            !write(*,*) i,P(mF+i),P(mF-i)
            FF(i)=energyDerivative*(P(mF+i)-P(mF-i))/(2.0d0*F(mF+i))
        end do
    else if (o_derivative.eq.2) then !compute finite-field's second derivative
        do i=1,mF-1
            !write(*,*) i,P(mF+i),P(mF),P(mF-i)
            FF(i)=energyDerivative*(P(mF+i)-2.0d0*P(mF)+P(mF-i))/F(mF+i)**2.0d0
        end do
    else if (o_derivative.eq.3) then !compute finite-field's third derivative
        !do i=1,totalFields
        !    write(*,*) P(i)
        !end do
        do i=1,mF-1
            if (i.eq.mF-1) cycle
            write(*,*) i,P(mF+i+1),P(mF+i),P(mF-i),P(mF-i-1)
            write(*,*) P(mF+i+1)-a*P(mF+i)+a*P(mF-i)-P(mF-i-1),a*(a**2.0d0-1)*F(mF+i)**3.0d0
            FF(i)=3.0d0*energyDerivative*(P(mF+i+1)-a*P(mF+i)+a*P(mF-i)-P(mF-i-1))/(a*(a**2.0d0-1)*F(mF+i)**3.0d0)
        end do
        !do i=1,mF-2
        !    write(*,*) FF(i)
        !end do
    else if (o_derivative.eq.4) then !compute finite-field's fourth derivative
        do i=1,mF-1
            if (i.eq.mF-2.or.i.eq.mF-1) cycle
            !write(*,*) i,P(mF+i+1)-(a**2.0d0)*P(mF+i)+(2.0d0*(a**2.0d0-1)*P(mF))-(a**2.0d0)*P(mF-i)+P(mF-i-1)
            !write(*,*) P(mF+i+1),P(mF+i),P(mF),P(mF-i),P(mF-i-1)
            FF(i)=12.0d0*energyDerivative*(P(mF+i+1)-(a**2.0d0)*P(mF+i)+(2.0d0*(a**2.0d0-1)*P(mF))-(a**2.0d0)*P(mF-i)+P(mF-i-1))/(a**2.0d0*(a**2.0d0-1)*F(mF+i)**4.0d0)
        end do
    end if

    write(*,*) "Field, properties and derivatives"
    do i=1,totalFields
        write(*,*) F(i),P(i),FF(i)
    end do
    write(*,*)

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
    integer :: iterRR,length
    integer :: i,j

    !-- Reference: M. Medved et al. / Journal of Molecular Structure: THEOCHEM 847 (2007) 39–46

    if (derivative_order.le.2) length=mF-1
    if (derivative_order.eq.3) length=mF-2
    if (derivative_order.eq.4) length=mF-3
    if (allocated(RombergT)) deallocate(RombergT); allocate(RombergT(length,length))
    if (allocated(errRombergT)) deallocate(errRombergT); allocate(errRombergT(length-1,length-1))
    if (allocated(errRomberg)) deallocate(errRomberg);  allocate(errRomberg(length-1,length-1))

        !-- Define Romberg matrix
    RombergT(1:length,1)=P
    RombergT(1:length,2:length)=9.9d99
    errRomberg=9.9d99

        !-- Compute & print the Romberg triangle
    do i=2,length   !-- Compute from the second column, the first is the "zero-th" Romberg iteration
        iterRR=i-1  !-- Compute the Romberg iteration
        do j=1,length-i+1 !-- Compute the Romberg for each row
            RombergT(j,i)=(a**(2*iterRR)*RombergT(j,i-1)-RombergT(j+1,i-1))/(a**(2*iterRR)-1.0d0)
        end do
    end do
    
    write(*,*) ("                            ROMBERG TRIANGLE",i=1,3)
    write(*,*) "Iteration:",(j-1,"              ",j=1,length)
    do i=1,length
        do j=1,length
            if (RombergT(i,j).eq.9.9d99) then
                write(*,'(A)',advance='NO') "   "
            else
                write(*,'(xxxxxxxx1pe22.15)',advance="no") RombergT(i,j)
            end if
        end do
        write(*,*)
    end do

        !-- Compute the Romberg absolute error matrix (errRomberg,intent(out)) and save:
        !       # the lowest (errRomberg1) error
        !       # the second-to-lowest (errRomberg2) error 
        !       # its respective value within 'errRomberg'
    write(*,*)
    do i=1,length-1
        do j=1,length-1
            if (RombergT(j,i).ge.9.9d99.or.Rombergt(j+1,i).ge.9.9d99) cycle
            errRomberg(j,i)=abs(RombergT(j,i)-RombergT(j+1,i))
        end do
    end do
    errRombergT=errRomberg

    write(*,*) ("                        ROMBERG ERROR TRIANGLE",i=1,3)
    write(*,*) "Iteration:",(j-1,"              ",j=1,mF-1)
    do i=1,length-1
        do j=1,length-1
            if (errRombergT(i,j).eq.9.9d99) then
                write(*,'(A)',advance='NO') "   "
            else
                write(*,'(xxxxxxxx1pe22.15)',advance="no") errRomberg(i,j)
            end if
        end do
        write(*,*)
    end do

        !-- Get the position of the value with the lowest and second-to-lowest errors
    abs_errRomberg1=minval(errRomberg)                          !-- Get the lowest value of the error triangle
        !err_minloc(1,1) == row, err_minloc(2,1) == column 
    err_minloc(1:2,1)=minloc(errRomberg)                        !-- Get its location
        errRomberg(err_minloc(1,1),err_minloc(2,1))=9.9d99      !-- Swap lowest error for 9.9d99
        err_minloc(2,1)=err_minloc(2,1)+1                       !-- Move it to the right
    abs_errRomberg2=minval(errRomberg)                          !-- Get the SECOND lowest value of the error triangle
    err_minloc(1:2,2)=minloc(errRomberg)                        !-- Get its location
        err_minloc(2,2)=err_minloc(2,2)+1                       !-- Move it to the right
    deallocate(errRomberg)

    write(*,*)
    write(*,'(" RomberG - Minimum value:",xF30.10,xx"and second best:",xF30.10)') RombergT(err_minloc(1,1),err_minloc(2,1)),RombergT(err_minloc(1,2),err_minloc(2,2))
    write(*,'(" RomberG - Minimum absolute (Romberg) errors for the properties:",xF30.10,xx"and:",xF30.10)') abs_errRomberg1,abs_errRomberg2
    !write(*,'(" RomberG - Romberg errors:",xF20.10,"% &&",xF20.10,"%")') &
    !& 1.0d2*abs_errRomberg1/RombergT(err_minloc(1,1),err_minloc(2,1)),1.0d2*abs_errRomberg2/RombergT(err_minloc(1,2),err_minloc(2,2))
    write(*,*)

    End subroutine RombergProcedure

    Subroutine Reps(isEnergy,printProperties,inlop,onlop,components,derivative_order,mF,totalFields,TensorP,secRombergP,mainRombergP)
    implicit none
    character*1,dimension(6) :: derivative_substr
    character :: dipoleComponents(3)*1, alphaComponents(6)*2,betaComponents(10)*3
    logical, intent(inout) :: isEnergy
    logical, intent(in) :: printProperties
    double precision, allocatable, intent(inout), dimension(:,:,:) :: mainRombergP,secRombergP,TensorP
    integer, intent(in) :: inlop,onlop,components,derivative_order,mF,totalFields
    double precision, dimension(3) :: secP,mainP
    integer :: i,j,k,n_dim
    
    derivative_substr=(/"x","y","z","X","Y","Z"/)
    dipoleComponents=(/"X  ","Y  ","Z  "/)
    alphaComponents=(/"XX ","XY ","YY ","XZ ","YZ ","ZZ "/)
    betaComponents=(/"XXX","XXY","XYY","YYY","XXZ","XYZ","YZZ","XZZ","YZZ","ZZZ"/)

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
    write(*,*)
    write(*,'(" RomberG - COMPUTED PROPERTIES FROM THE ROMBERG PROCEDURE")')
    write(*,*)
    do k=1,3
        write(*,*)
        if (inlop.eq.1) write(*,'(" RomberG - Components of μ:",1X,A1,26X,A1,26X,A1)') dipoleComponents(1:components-1) 
        if (inlop.eq.2) write(*,'(" RomberG - Components of α:",1X,A2,26X,A2,26X,A2,25X,A2,23X,A2,22X,A2)') alphaComponents(1:components-1) 
        if (inlop.eq.3) write(*,'(" RomberG - Components of β:",1X,A3,16X,A3,17X,A3,18X,A3,19X,A3,16X,A3,17X,A3,18X,A3,19X,A3,19X,A3)') betaComponents(1:components-1)

        if (inlop.ne.3) write(*,*) "d",trim(derivative_substr(k+3)),"**",derivative_order," -    Value=",(mainRombergP(j,1,k),j=1,components-1)
        if (inlop.ne.3) write(*,*) "d",trim(derivative_substr(k+3)),"**",derivative_order," - R. Error=",(mainRombergP(j,2,k),j=1,components-1)

        if (inlop.eq.3) write(*,'("d",A1,"**",I1," -    Value=",1X,1pe20.10,1X,1pe20.10,1X,1pe20.10,1X,1pe20.10,1X,1pe20.10,1X,1pe20.10,1X,1pe20.10,1X,1pe20.10,1X,1pe20.10,1X,1pe20.10)') &
        & derivative_substr(k+3),derivative_order,mainRombergP(1:components-1,1,k)
        if (inlop.eq.3) write(*,'("d",A1,"**",I1," - R. Error=",1X,1pe20.10,1X,1pe20.10,1X,1pe20.10,1X,1pe20.10,1X,1pe20.10,1X,1pe20.10,1X,1pe20.10,1X,1pe20.10,1X,1pe20.10,1X,1pe20.10)') &
        & derivative_substr(k+3),derivative_order,mainRombergP(1:components-1,2,k)
        !write(*,*) derivative_substr(k+3),"  RomberG error(%)=",(mainRombergP(j,3,k),j=1,components-1)
    end do
    if (printProperties.eqv..TRUE.) then
        do k=1,3
            write(*,*)
            if (inlop.eq.1) write(*,*) ("                         ",dipoleComponents(i),i=1,components-1) 
            if (inlop.eq.2) write(*,*) ("                         ",alphaComponents(i),i=1,components-1) 
            if (inlop.eq.3) write(*,*) ("                         ",betaComponents(i),i=1,components-1) 
            write(*,*) derivative_substr(k+3),"             Value=",(secRombergP(j,1,k),j=1,components-1)
            write(*,*) derivative_substr(k+3),"    Absolute error=",(secRombergP(j,2,k),j=1,components-1)
            !write(*,*) derivative_substr(k+3),"  RomberG error(%)=",(secRombergP(j,3,k),j=1,components-1)
        end do
    end if
    write(*,*)

    End subroutine
    
    Subroutine complex2simple(TensRomP)
    implicit none
    double precision, intent(inout), allocatable, dimension(:,:,:) :: TensRomP
    double precision :: threshold,marking
    integer :: i,j,rows,columns

    columns=3
    rows=size(TensRomP)/columns/3
    threshold=1.0d1

    do i=1,columns
        do j=1,rows
            marking=1.0d2*abs(TensRomP(j,1,i)/TensRomp(j,2,i))
            if (marking.lt.threshold) then
                write(*,'(" RomberG - WARNING! The output tensor components (",xI1,",",xI1,") have significant errors.")') j,i
            else if (marking.gt.threshold.and.abs(TensRomP(j,1,i)).le.1.0e-6) then
                write(*,'(" RomberG - WARNING! The output tensor component (",xI1,",",xI1,") is simplified to zero. ")') j,i
                write(*,'(" RomberG - To avoid this, execute again RomberG without the -S flag")')
                TensRomP(j,1,i)=0.0d0
            end if
        end do
    end do
    End subroutine complex2simple


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
logical :: simpleOutput = .FALSE.
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
    !-- Best tensors for each property
double precision, dimension (5,3) :: bestGamma
double precision, dimension (5,2) :: bestBeta
double precision, dimension (6) :: bestAlpha
double precision, dimension (3) :: bestDipole,tmpDipole
double precision :: BestProp
    !-- Output values
double precision, dimension (3,3,3,3) :: Gamma
double precision, dimension (3,3,3) :: Beta
double precision, dimension (3,3) :: Alpha
double precision :: dipole_module,tmpDipModulus,alpha_average,beta_vec,beta_4,beta_parallel,gamma_parallel,gamma_average,deltaAlpha,beta_perpendicular
    !-- Dummy matrix to store the properties 
double precision, dimension (5,2) :: P_general
double precision :: field_strength
integer :: indexEnergy,indexDipole,indexAlpha,indexBeta,indexGamma
integer :: Xindex,Yindex,Zindex,indexBase,indexField
integer :: stat1,ios,err,dummy,errorType
integer :: i,j,k

type(option_s) :: opts(9)
!-- opts(i)=option_s(long_name,arguments?,short_name)
opts(1) = option_s("input",.TRUE.,"i")
opts(2) = option_s("output",.TRUE.,"o")
opts(3) = option_s("longitudinal",.FALSE.,"L")
opts(4) = option_s("isotropic",.FALSE.,"I")
opts(5) = option_s("totalfields",.TRUE.,"F")
opts(6) = option_s("printP",.FALSE.,"p")
opts(7) = option_s("help",.FALSE.,"h")
opts(8) = option_s("total",.FALSE.,"T")
opts(9) = option_s("simple",.FALSE.,"S")

derivative_substr=(/"x","y","z","X","Y","Z"/)
dipoleComponents=(/"X  ","Y  ","Z  "/)
alphaComponents=(/"XX ","XY ","YY ","XZ ","YZ ","ZZ "/)
betaComponents=(/"XXX","XXY","XYY","YYY","XXZ","XYZ","YZZ","XZZ","YZZ","ZZZ"/)

if (command_argument_count().le.1) then
    write(*,*)
    write(*,*) "./ROMBERG.exe {options} {name_of_the_molecule}"
    write(*,*)
    write(*,*) " To display the help message execute:"
    write(*,*) "./ROMBERG.exe --help {name_of_the_molecule}"
    write(*,*)
    stop
end if

    !-- Get command line arguments
do
    arg=getopt("i:o:F:SLITph",opts)
        !-- Implementar "Simple"
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

        case("F","totalfields")
            read(optarg,*) totalFields
            positiveFields=(totalFields-1)/2
            negativeFields=positiveFields
            mF=positiveFields+1

        case("S","simple")
            simpleOutput=.TRUE.

        case("L","longitudinal")
            doLongitudinal=.TRUE.

        case("I","isotropic")
            doIsotropic=.TRUE.

        case("T","total")
            doLongitudinal=.TRUE.
            doIsotropic=.TRUE.

        case("p","printP")
            printProperties=.TRUE.

        case("h","help")
            write(*,*)
            write(*,*) " Several optionalities are available for RomberG.exe:"
            write(*,*) "    -i, --input          :   Input property - Energy/energy/E ; Dipole/dipole/M ; Alpha/alpha/A ; Beta/beta/B"
            write(*,*) "    -o, --output         :   Output property - Dipole/dipole/M ; Alpha/alpha/A ; Beta/beta/B ; Gamma/gamma/G"
            write(*,*) "    -F, --totalfields    :   Number of total fields probed. Integer required."
            write(*,*) "    -l, --longitudinal   :   Whether to compute longitudinal properties or not."
            write(*,*) "    -I, --isotropic      :   Whether to compute isotropic properties or not."
            write(*,*) "    -S, --simple         :   Print the output matrices without losing precision"
            write(*,*) "    -p, --printP         :   Print the derivatives when computing the Romberg triangle for each component."
            write(*,*) "    -h, --help           :   Displays this message."
            write(*,*)
            write(*,*) " REMEMBER: The execution of RomberG.exe requires the name of the molecule, even for displaying this message."
            write(*,*) 
            stop
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
    call system("cd $(pwd)/"//mol_name//" ; grep -m 1 -A 2 'Polarizability' *.fchk > isAlpha.nlop; sed -i '/Polarizability/d' isAlpha.nlop; cp isAlpha.nlop ../")
    open (unit=4,file="isAlpha.nlop",status="old")
else if (inlop.eq.3) then
        !-- Remember that for beta the sign is changed with respect to the regular convetion!
    isEnergy=.FALSE.
    call system("rm isBeta.nlop")
    call system("cd $(pwd)/"//mol_name//" ; grep -A 2 'HyperPolarizability' *.fchk > isBeta.nlop; sed -i '/HyperPolarizability/d' isBeta.nlop; cp isBeta.nlop ../")
    open (unit=4,file="isBeta.nlop",status="old")
end if

    !-- Read the zero-field name and NLOPs
                read(4,*,iostat=ios) line,P_general(1,1),P_general(2,1),P_general(3,1),P_general(4,1),P_general(5,1)
if (inlop.ge.2) read(4,*,iostat=ios) line,P_general(1,2),P_general(2,2),P_general(3,2),P_general(4,2),P_general(5,2)
if (inlop.eq.3) read(4,*,iostat=ios) line

    !-- Declare the zero-field properties according to Gaussian's .fchk format
        !-- He d'intentar millorar aquesta part del codi
dummy=1
if (inlop.eq.3) then
        !-- Beta is changed to negative later 
    allocate(P_beta(totalFields,11,3))
    P_beta(mF,1,1:3)=0.0d0
    do j=1,2; do i=1,5
            dummy=dummy+1
            P_beta(mF,dummy,1:3)=P_general(i,j)
    end do; end do
    !do i=1,3
    !    write(*,*) (P_beta(mF,j,i),j=1,dummy)
    !end do 
    !stop
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
write(*,*) "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
write(*,*)

!#########################################################
!BestGamma     BestBeta     BestAlpha     BestDipole
!xxxx          xxx          xx            x
!xxyx          xxy          yx            y
!xxyy          yxy          yy            z
!yxyy          yyy          zx
!yyyy          xxz          zy

!xxzx          yxz          zz
!xxzy          yyz
!yxzy          zxz
!yyzy          zyz
!xxzz          zzz

!yxzz
!yyzz
!zxzz
!zyzz
!zzzz
!#########################################################

    !-- (components+1) because the Romberg data is stored from the second column onwards
errorType=1 !-- errorType.eq.2 uses RombergError, errorType.eq.1 uses AbsoluteError
if (inlop.eq.3) then !-- Compute the derivatives of beta
    if (allocated(mainRombergBeta)) deallocate(mainRombergBeta)
    if (allocated(secRombergBeta)) deallocate(secRombergBeta)
    allocate(mainRombergBeta(10,3,3)); allocate(secRombergBeta(10,3,3))
   !call Reps(isEnergy,printProperties,inlop,onlop,components,derivative_order,mF,totalFields,TensorP,secRombergP,mainRombergP)
    call Reps(isEnergy,printProperties,inlop,onlop,11,derivative_order,mF,totalFields,P_beta,secRombergBeta,mainRombergBeta)

    if (simpleOutput.eqv..TRUE.) call complex2simple(mainRombergBeta)

        !-- Get the best property tensor from the Romberg values and errors (errorType=1)
    if (onlop.eq.4) then
        BestGamma=0.0d0
        BestGamma(1,1)=mainRombergBeta(1,1,1) ! \gamma_xxxx
        BestGamma(2,1)=BestProp(mainRombergBeta(2,1,1),mainRombergBeta(1,1,2),mainRombergBeta(2,errorType+1,1),mainRombergBeta(1,errorType+1,2)) ! \gamma_xxyx
        BestGamma(3,1)=BestProp(mainRombergBeta(3,1,1),mainRombergBeta(2,1,2),mainRombergBeta(3,errorType+1,1),mainRombergBeta(2,errorType+1,2)) ! \gamma_xxyy
        BestGamma(4,1)=BestProp(mainRombergBeta(4,1,1),mainRombergBeta(2,1,3),mainRombergBeta(4,errorType+1,1),mainRombergBeta(2,errorType+1,3)) ! \gamma_yxyy
        BestGamma(5,1)=mainRombergBeta(4,1,2) ! \gamma_yyyy
        BestGamma(1,2)=BestProp(mainRombergBeta(5,1,1),mainRombergBeta(1,1,3),mainRombergBeta(5,errorType+1,1),mainRombergBeta(1,errorType+1,3)) ! \gamma_xxzx
        BestGamma(2,2)=mainRombergBeta(6,1,1) ! \gamma_xxzy; choose the \gamma_xxzy with the lowest Romberg error
            if(abs(mainRombergBeta(6,errorType+1,1)).gt.abs(mainRombergBeta(5,errorType+1,2))) BestGamma(2,2)=mainRombergBeta(5,1,2)
            if(abs(mainRombergBeta(6,errorType+1,1)).gt.abs(mainRombergBeta(2,errorType+1,3))) BestGamma(2,2)=mainRombergBeta(2,1,3)
            if(abs(mainRombergBeta(5,errorType+1,2)).gt.abs(mainRombergBeta(2,errorType+1,3))) BestGamma(2,2)=mainRombergBeta(2,1,3)
        BestGamma(3,2)=mainRombergBeta(7,1,1) ! \gamma_yxzy; choose the \gamma_yxzy with the lowest Romberg error
            if(abs(mainRombergBeta(7,errorType+1,1)).gt.abs(mainRombergBeta(6,errorType+1,2))) BestGamma(3,2)=mainRombergBeta(6,1,2)
            if(abs(mainRombergBeta(7,errorType+1,1)).gt.abs(mainRombergBeta(3,errorType+1,3))) BestGamma(3,2)=mainRombergBeta(3,1,3)
            if(abs(mainRombergBeta(6,errorType+1,2)).gt.abs(mainRombergBeta(3,errorType+1,3))) BestGamma(3,2)=mainRombergBeta(3,1,3)
        BestGamma(4,2)=BestProp(mainRombergBeta(7,1,2),mainRombergBeta(4,1,3),mainRombergBeta(2,errorType+1,1),mainRombergBeta(1,errorType+1,2)) ! \gamma_yyzy
        BestGamma(5,2)=BestProp(mainRombergBeta(8,1,1),mainRombergBeta(5,1,3),mainRombergBeta(2,errorType+1,1),mainRombergBeta(1,errorType+1,2)) ! \gamma_xxzz
        BestGamma(1,3)=mainRombergBeta(6,1,3) ! \gamma_yxzz; choose the \gamma_yxzz with the lowest Romberg error
            if(abs(mainRombergBeta(6,errorType+1,3)).gt.abs(mainRombergBeta(8,errorType+1,2))) BestGamma(1,3)=mainRombergBeta(8,1,2)
            if(abs(mainRombergBeta(6,errorType+1,3)).gt.abs(mainRombergBeta(9,errorType+1,1))) BestGamma(1,3)=mainRombergBeta(9,1,1)
            if(abs(mainRombergBeta(8,errorType+1,2)).gt.abs(mainRombergBeta(9,errorType+1,1))) BestGamma(1,3)=mainRombergBeta(9,1,1)
        BestGamma(2,3)=BestProp(mainRombergBeta(9,1,2),mainRombergBeta(7,1,3),mainRombergBeta(2,errorType+1,1),mainRombergBeta(1,errorType+1,2)) ! \gamma_yyzz
        BestGamma(3,3)=BestProp(mainRombergBeta(8,1,3),mainRombergBeta(10,1,1),mainRombergBeta(2,errorType+1,1),mainRombergBeta(1,errorType+1,2)) ! \gamma_zxzz
        BestGamma(4,3)=BestProp(mainRombergBeta(10,1,2),mainRombergBeta(9,1,3),mainRombergBeta(2,errorType+1,1),mainRombergBeta(1,errorType+1,2)) ! \gamma_zyzz
        BestGamma(5,3)=mainRombergBeta(10,1,3) ! \gamma_zzzz
    end if

            !-- Transform the BestBeta matrix into the first hyperpolarizability tensor
    call Gamma2Gamma(BestGamma,Gamma)

else if (inlop.eq.2) then !-- Compute the derivatives of alpha

    if (allocated(mainRombergAlpha)) deallocate(mainRombergAlpha)
    if (allocated(secRombergAlpha)) deallocate(secRombergAlpha)
    allocate(mainRombergAlpha(6,3,3)); allocate(secRombergAlpha(6,3,3))
   !call Reps(isEnergy,printProperties,inlop,onlop,components,derivative_order,mF,totalFields,TensorP,secRombergP,mainRombergP)
    call Reps(isEnergy,printProperties,inlop,onlop,7,derivative_order,mF,totalFields,P_alpha,secRombergAlpha,mainRombergAlpha)

    if (simpleOutput.eqv..TRUE.) call complex2simple(mainRombergAlpha)

        !-- Get the best property tensor from the Romberg values and errors (errorType=2)
    if (onlop.eq.3) then
        BestBeta(1,1)=mainRombergAlpha(1,1,1) ! \beta_xxx
        BestBeta(2,1)=BestProp(mainRombergAlpha(2,1,1),mainRombergAlpha(1,1,2),mainRombergAlpha(2,errorType+1,1),mainRombergAlpha(1,errorType+1,2)) ! \beta_xxy
        BestBeta(3,1)=BestProp(mainRombergAlpha(3,1,1),mainRombergAlpha(2,1,2),mainRombergAlpha(3,errorType+1,1),mainRombergAlpha(2,errorType+1,2)) ! \beta_xyy
        BestBeta(4,1)=mainRombergAlpha(3,1,2) ! \beta_yyy
        BestBeta(5,1)=BestProp(mainRombergAlpha(4,1,1),mainRombergAlpha(1,1,3),mainRombergAlpha(4,errorType+1,1),mainRombergAlpha(1,errorType+1,4)) ! \beta_xxz
        BestBeta(1,2)=mainRombergAlpha(5,1,1) ! \beta_xyz
            if (abs(mainRombergAlpha(5,errorType+1,1)).gt.abs(mainRombergAlpha(4,errorType+1,2))) BestBeta(1,2)=mainRombergAlpha(4,1,2)
            if (abs(mainRombergAlpha(5,errorType+1,1)).gt.abs(mainRombergAlpha(2,errorType+1,3))) BestBeta(1,2)=mainRombergAlpha(2,1,3)
            if (abs(mainRombergAlpha(4,errorType+1,2)).gt.abs(mainRombergAlpha(2,errorType+1,3))) BestBeta(1,2)=mainRombergAlpha(2,1,3)
        BestBeta(2,2)=BestProp(mainRombergAlpha(5,1,2),mainRombergAlpha(3,1,3),mainRombergAlpha(5,errorType+1,2),mainRombergAlpha(3,errorType+1,3)) ! \beta_yyz
        BestBeta(3,2)=BestProp(mainRombergAlpha(6,1,1),mainRombergAlpha(4,1,3),mainRombergAlpha(6,errorType+1,1),mainRombergAlpha(4,errorType+1,3)) ! \beta_xzz
        BestBeta(4,2)=BestProp(mainRombergAlpha(6,1,2),mainRombergAlpha(5,1,3),mainRombergAlpha(6,errorType+1,2),mainRombergAlpha(5,errorType+1,3)) ! \beta_yzz
        BestBeta(5,2)=mainRombergAlpha(6,1,3) ! \beta_zzz

            !-- Transform the BestBeta matrix into the first hyperpolarizability tensor
        call Beta2Beta(BestBeta,Beta)
 
    else if (onlop.eq.4) then
        BestGamma=0.0d0
        BestGamma(1,1)=mainRombergAlpha(1,1,1) ! \gamma_xxxx
        BestGamma(3,1)=BestProp(mainRombergAlpha(1,1,2),mainRombergAlpha(3,1,1),mainRombergAlpha(1,errorType+1,2),mainRombergAlpha(3,errorType+1,1)) ! \gamma_xxyy
        BestGamma(5,1)=mainRombergAlpha(3,1,2) ! \gamma_yyyy
        BestGamma(5,2)=BestProp(mainRombergAlpha(1,1,3),mainRombergAlpha(6,1,1),mainRombergAlpha(1,errorType+1,3),mainRombergAlpha(6,errorType+1,1)) ! \gamma_xxzz
        BestGamma(2,3)=BestProp(mainRombergAlpha(3,1,3),mainRombergAlpha(6,1,2),mainRombergAlpha(3,errorType+1,3),mainRombergAlpha(6,errorType+1,2)) ! \gamma_yyzz
        BestGamma(5,3)=mainRombergAlpha(6,1,3) ! \gamma_zzzz

            !-- Transform the BestGamma matrix into the second hyperpolarizability tensor
        call Gamma2Gamma(BestGamma,Gamma)

    end if

else if (inlop.eq.1) then !-- Compute the derivatives of the dipole moment

        !-- This section has to be revised
    if (allocated(mainRombergDipole)) deallocate(mainRombergDipole)
    if (allocated(secRombergDipole)) deallocate(secRombergDipole)
    allocate (mainRombergDipole(3,3,3)); allocate(secRombergDipole(3,3,3))
   !call Reps(isEnergy,printProperties,inlop,onlop,components,derivative_order,mF,totalFields,TensorP,secRombergP,mainRombergP)
    call Reps(isEnergy,printProperties,inlop,onlop,4,derivative_order,mF,totalFields,P_dipole,secRombergDipole,mainRombergDipole)

    if (simpleOutput.eqv..TRUE.) call complex2simple(mainRombergDipole)

        !-- Get the best property tensor from the Romberg values and errors (errorType=2)
    if (onlop.eq.2) then
        BestAlpha(1)=mainRombergDipole(1,1,1) ! \alpha_xx
        BestAlpha(2)=BestProp(mainRombergDipole(1,1,2),mainRombergDipole(2,1,1),mainRombergDipole(1,errorType+1,2),mainRombergDipole(2,errorType+1,1)) ! \alpha_xy
        BestAlpha(3)=mainRombergDipole(2,1,2) ! \alpha_yy
        BestAlpha(4)=BestProp(mainRombergDipole(1,1,3),mainRombergDipole(3,1,1),mainRombergDipole(1,errorType+1,3),mainRombergDipole(3,errorType+1,1)) ! \alpha_xz
        BestAlpha(5)=BestProp(mainRombergDipole(2,1,3),mainRombergDipole(3,1,2),mainRombergDipole(2,errorType+1,3),mainRombergDipole(3,errorType+1,2)) ! \alpha_yz
        BestAlpha(6)=mainRombergDipole(3,1,3) ! \alpha_zz

            !-- Transform the BestAlpha vector into the symmetric polarizability matrix
        call Alpha2Alpha(BestAlpha,Alpha)

    else if (onlop.eq.3) then
        BestBeta=0.0d0
        BestBeta(1,1)=mainRombergDipole(1,1,1) ! \beta_xxx
        BestBeta(2,1)=mainRombergDipole(2,1,1) ! \beta_xxy
        BestBeta(3,1)=mainRombergDipole(1,1,2) ! \beta_xyy
        BestBeta(4,1)=mainRombergDipole(2,1,2) ! \beta_yyy
        BestBeta(5,1)=mainRombergDipole(3,1,1) ! \beta_xxz
        BestBeta(2,2)=mainRombergDipole(3,1,2) ! \beta_yyz
        BestBeta(3,2)=mainRombergDipole(1,1,3) ! \beta_xzz
        BestBeta(4,2)=mainRombergDipole(2,1,3) ! \beta_yzz
        BestBeta(5,2)=mainRombergDipole(3,1,3) ! \beta_zzz 

            !-- Transform the BestBeta matrix into the symmetric hyperpolarizability tensor
        call Beta2Beta(BestBeta,Beta)

    else if (onlop.eq.4) then
        BestGamma=0.0d0
        BestGamma(1,1)=mainRombergDipole(1,1,1) ! \gamma_xxxx
        BestGamma(2,1)=mainRombergDipole(2,1,1) ! \gamma_xxxy
        BestGamma(4,1)=mainRombergDipole(1,1,2) ! \gamma_xyyy
        BestGamma(5,1)=mainRombergDipole(2,1,2) ! \gamma_yyyy
        BestGamma(1,2)=mainRombergDipole(3,1,1) ! \gamma_xxxz
        BestGamma(4,2)=mainRombergDipole(3,1,2) ! \gamma_yyyz
        BestGamma(3,3)=mainRombergDipole(1,1,3) ! \gamma_xzzz
        BestGamma(4,3)=mainRombergDipole(2,1,3) ! \gamma_yzzz
        BestGamma(5,3)=mainRombergDipole(3,1,3) ! \gamma_zzzz

            !-- Transform the BestGamma matrix into the symmetric hyperpolarizability second order tensor
        call Gamma2Gamma(BestGamma,Gamma)

    end if

else if (inlop.eq.0) then !-- Compute the derivatives of the energy

    if (allocated(mainRombergEnergy)) deallocate(mainRombergEnergy)
    if (allocated(secRombergEnergy)) deallocate(secRombergEnergy)
    allocate (mainRombergEnergy(3,3,3)); allocate(secRombergEnergy(3,3,3))
   !call Reps(isEnergy,printProperties,inlop,onlop,components,derivative_order,mF,totalFields,TensorP,secRombergP,mainRombergP)
    call Reps(isEnergy,printProperties,inlop,onlop,2,derivative_order,mF,totalFields,P_energy,secRombergEnergy,mainRombergEnergy)

    if (simpleOutput.eqv..TRUE.) call complex2simple(mainRombergEnergy)

        !-- Get the best property tensor from the Romberg values and errors (errorType=2)
    if (onlop.eq.1) then ! \mu_x/y/z
        BestDipole(1)=mainRombergEnergy(1,1,1)
        BestDipole(2)=mainRombergEnergy(1,1,2)
        BestDipole(3)=mainRombergEnergy(1,1,3) 

    else if (onlop.eq.2) then ! \alpha_xx/yy/zz
        BestAlpha=0.0d0
        BestAlpha(1)=mainRombergEnergy(1,1,1)
        BestAlpha(3)=mainRombergEnergy(1,1,2)
        BestAlpha(6)=mainRombergEnergy(1,1,3) 

            !-- Transform the BestAlpha vector into the symmetric polarizability matrix
        call Alpha2Alpha(BestAlpha,Alpha)

    else if (onlop.eq.3) then ! \beta_xxx/yyy/zzz
        BestBeta=0.0d0
        BestBeta(1,1)=mainRombergEnergy(1,1,1)
        BestBeta(4,1)=mainRombergEnergy(1,1,2)
        BestBeta(5,2)=mainRombergEnergy(1,1,3) 

            !-- Transform the BestBeta matrix into the symmetric hyperpolarizability tensor
        call Beta2Beta(BestBeta,Beta)

    else if (onlop.eq.4) then ! \gamma_xxxx/yyyy/zzzz
        BestGamma=0.0d0
        BestGamma(1,1)=mainRombergEnergy(1,1,1)
        BestGamma(5,1)=mainRombergEnergy(1,1,2)
        BestGamma(5,3)=mainRombergEnergy(1,1,3) 

            !-- Transform the BestGamma matrix into the symmetric hyperpolarizability second order tensor
        call Gamma2Gamma(BestGamma,Gamma)

    end if

end if

write(*,*)
write(*,*) "######################################################"
write(*,*) "######################################################"
write(*,*) "######################################################"
write(*,*)

    !-- Computing the longitudinal properties and other common values
!dipole_module,tmpDipModulus,alpha_average,beta_vec,beta_4,beta_parallel,gamma_parallel
if (doLongitudinal.eqv..TRUE..or.doIsotropic.eqv..TRUE.) then
    if (onlop.eq.1) then
            !-- Get the longitudinal vector
        dipole_module=dsqrt(BestDipole(1)**2.0d0+BestDipole(2)**2.0d0+BestDipole(3)**2.0d0)
        write(*,'(" RomberG - Dipole moment vector (μX,μY,μZ) = ("xF10.6,",",xF10.6,",",xF10.6,")")') BestDipole(1),BestDipole(2),BestDipole(3)
        write(*,'(" RomberG - The modulus of the dipole moment is",xF10.6" a.u.")') dipole_module

    else if (onlop.eq.2) then
            !-- Get the isotropic polarizability: 10.1021/jp405144f
        alpha_average=(BestAlpha(1)+BestAlpha(3)+BestAlpha(6))/3.0d0
        write(*,'(" RomberG - Elements of the polarizability matrix (αXX,αXY,αYY,αXZ,αYZ) = ("xF10.6,",",xF10.6,",",xF10.6,",",xF10.6,",",xF10.6)') &
        & BestAlpha(1),BestAlpha(2),BestAlpha(3),BestAlpha(4),BestAlpha(5)
        write(*,'(" RomberG - Elements of the polarizability matrix                 (αZZ) =",xF10.6,")")') BestAlpha(6)
        write(*,*)
        write(*,'(" RomberG - Polarizability matrix (αXX,αXY,αXZ):",4X,1pe22.15,4X,1pe22.15,4X,1pe22.15)') Alpha(1,1),Alpha(1,2),Alpha(1,3)
        write(*,'(" RomberG - Polarizability matrix (αYX,αYY,αYZ):",4X,1pe22.15,4X,1pe22.15,4X,1pe22.15)') Alpha(2,1),Alpha(2,2),Alpha(2,3)
        write(*,'(" RomberG - Polarizability matrix (αZX,αZY,αZZ):",4X,1pe22.15,4X,1pe22.15,4X,1pe22.15)') Alpha(3,1),Alpha(3,2),Alpha(3,3)
        write(*,*)
        write(*,'(" RomberG - The average polarizability is ",x1pe22.15" a.u.")') alpha_average

    else if (onlop.eq.3) then
            !-- Get the longitudinal first hyperpolarizability
        write(*,'(" RomberG - Elements of the hyperpolarizability tensor (βXXX,βXXY,βXYY,βYYY,βXXZ) = ("xF10.6,",",xF10.6,",",xF10.6,",",xF10.6,",",xF10.6,")")') &
        & BestBeta(1,1),BestBeta(2,1),BestBeta(3,1),BestBeta(4,1),BestBeta(5,1)
        write(*,'(" RomberG - Elements of the hyperpolarizability tensor (βXYZ,βYYZ,βXZZ,βYZZ,βZZZ) = ("xF10.6,",",xF10.6,",",xF10.6,",",xF10.6,",",xF10.6,")")') &
        & BestBeta(1,2),BestBeta(2,2),BestBeta(3,2),BestBeta(4,2),BestBeta(5,2)

            !-- Print the whole 1st hyperpolarizability tensor
        write(*,*)
        write(*,'(" RomberG - X-layer of the hyperpolarizability tensor (βXXX,βXYX,βXZX):",4X,1pe22.15,4X,1pe22.15,4X,1pe22.15)') Beta(1,1,1),Beta(1,2,1),Beta(1,3,1)
        write(*,'(" RomberG - X-layer of the hyperpolarizability tensor (βYXX,βYYX,βYZX):",4X,1pe22.15,4X,1pe22.15,4X,1pe22.15)') Beta(2,1,1),Beta(2,2,1),Beta(2,3,1)
        write(*,'(" RomberG - X-layer of the hyperpolarizability tensor (βZXX,βZYX,βZZX):",4X,1pe22.15,4X,1pe22.15,4X,1pe22.15)') Beta(3,1,1),Beta(3,2,1),Beta(3,3,1)
        write(*,*)
        write(*,'(" RomberG - Y-layer of the hyperpolarizability tensor (βXXY,βXYY,βXZY):",4X,1pe22.15,4X,1pe22.15,4X,1pe22.15)') Beta(1,1,2),Beta(1,2,2),Beta(1,3,2)
        write(*,'(" RomberG - Y-layer of the hyperpolarizability tensor (βYXY,βYYY,βYZY):",4X,1pe22.15,4X,1pe22.15,4X,1pe22.15)') Beta(2,1,2),Beta(2,2,2),Beta(2,3,2)
        write(*,'(" RomberG - Y-layer of the hyperpolarizability tensor (βZXY,βZYY,βZZY):",4X,1pe22.15,4X,1pe22.15,4X,1pe22.15)') Beta(3,1,2),Beta(3,2,2),Beta(3,3,2)
        write(*,*)
        write(*,'(" RomberG - Z-layer of the hyperpolarizability tensor (βXXZ,βXYZ,βXZZ):",4X,1pe22.15,4X,1pe22.15,4X,1pe22.15)') Beta(1,1,3),Beta(1,2,3),Beta(1,3,3)
        write(*,'(" RomberG - Z-layer of the hyperpolarizability tensor (βYXZ,βYYZ,βYZZ):",4X,1pe22.15,4X,1pe22.15,4X,1pe22.15)') Beta(2,1,3),Beta(2,2,3),Beta(2,3,3)
        write(*,'(" RomberG - Z-layer of the hyperpolarizability tensor (βZXZ,βZYZ,βZZZ):",4X,1pe22.15,4X,1pe22.15,4X,1pe22.15)') Beta(3,1,3),Beta(3,2,3),Beta(3,3,3)
        write(*,*)

    else if (onlop.eq.4) then
            !-- Get the longitudinal second hyperpolarizability
        write(*,'(" RomberG - Elements of the second hyperpolarizability tensor (γXXXX,γXXYX,γXXYY,γYXYY,γYYYY) = ("xF15.6,",",xF15.6,",",xF15.6,",",xF15.6,",",xF15.6,")")') &
        & BestGamma(1,1),BestGamma(2,1),BestGamma(3,1),BestGamma(4,1),BestGamma(5,1)
        write(*,'(" RomberG - Elements of the second hyperpolarizability tensor (γXXZX,γXXZY,γYXZY,γYYZY,γXXZZ) = ("xF15.6,",",xF15.6,",",xF15.6,",",xF15.6,",",xF15.6,")")') &
        & BestGamma(1,2),BestGamma(2,2),BestGamma(3,2),BestGamma(4,2),BestGamma(5,2)
        write(*,'(" RomberG - Elements of the second hyperpolarizability tensor (γYXZZ,γYYZZ,γZXZZ,γZYZZ,γZZZZ) = ("xF15.6,",",xF15.6,",",xF15.6,",",xF15.6,",",xF15.6,")")') &
        & BestGamma(1,3),BestGamma(2,3),BestGamma(3,3),BestGamma(4,3),BestGamma(5,3)

            !-- Print the whole 2nd hyperpolarizability tensor
        write(*,*)
        write(*,'(" RomberG - XX-layer of the second hyperpolarizability tensor (γXXXX,γXYXX,γXZXX):",4X,1pe22.15,4X,1pe22.15,4X,1pe22.15)') Gamma(1,1,1,1),Gamma(1,2,1,1),Gamma(1,3,1,1)
        write(*,'(" RomberG - XX-layer of the second hyperpolarizability tensor (γYXXX,γYYXX,γYZXX):",4X,1pe22.15,4X,1pe22.15,4X,1pe22.15)') Gamma(2,1,1,1),Gamma(2,2,1,1),Gamma(2,3,1,1)
        write(*,'(" RomberG - XX-layer of the second hyperpolarizability tensor (γZXXX,γZYXX,γZZXX):",4X,1pe22.15,4X,1pe22.15,4X,1pe22.15)') Gamma(3,1,1,1),Gamma(3,2,1,1),Gamma(3,3,1,1)
        write(*,*)
        write(*,'(" RomberG - XY-layer of the second hyperpolarizability tensor (γXXXY,γXYXY,γXZXY):",4X,1pe22.15,4X,1pe22.15,4X,1pe22.15)') Gamma(1,1,2,2),Gamma(1,2,2,2),Gamma(1,3,2,2)
        write(*,'(" RomberG - XY-layer of the second hyperpolarizability tensor (γYXXY,γYYXY,γYZXY):",4X,1pe22.15,4X,1pe22.15,4X,1pe22.15)') Gamma(2,1,2,2),Gamma(2,2,2,2),Gamma(2,3,2,2)
        write(*,'(" RomberG - XY-layer of the second hyperpolarizability tensor (γZXXY,γZYXY,γZZXY):",4X,1pe22.15,4X,1pe22.15,4X,1pe22.15)') Gamma(3,1,2,2),Gamma(3,2,2,2),Gamma(3,3,2,2)
        write(*,*)
        write(*,'(" RomberG - XZ-layer of the second hyperpolarizability tensor (γXXXZ,γXYXZ,γXZXZ):",4X,1pe22.15,4X,1pe22.15,4X,1pe22.15)') Gamma(1,1,3,3),Gamma(1,2,3,3),Gamma(1,3,3,3)
        write(*,'(" RomberG - XZ-layer of the second hyperpolarizability tensor (γYXXZ,γYYXZ,γYZXZ):",4X,1pe22.15,4X,1pe22.15,4X,1pe22.15)') Gamma(2,1,3,3),Gamma(2,2,3,3),Gamma(2,3,3,3)
        write(*,'(" RomberG - XZ-layer of the second hyperpolarizability tensor (γZXXZ,γZYXZ,γZZXZ):",4X,1pe22.15,4X,1pe22.15,4X,1pe22.15)') Gamma(3,1,3,3),Gamma(3,2,3,3),Gamma(3,3,3,3)
        write(*,*)
        write(*,'(" RomberG - YY-layer of the second hyperpolarizability tensor (γXXYY,γXYYY,γXZYY):",4X,1pe22.15,4X,1pe22.15,4X,1pe22.15)') Gamma(1,1,1,1),Gamma(1,2,1,1),Gamma(1,3,1,1)
        write(*,'(" RomberG - YY-layer of the second hyperpolarizability tensor (γYXYY,γYYYY,γYZYY):",4X,1pe22.15,4X,1pe22.15,4X,1pe22.15)') Gamma(2,1,1,1),Gamma(2,2,1,1),Gamma(2,3,1,1)
        write(*,'(" RomberG - YY-layer of the second hyperpolarizability tensor (γZXYY,γZYYY,γZZYY):",4X,1pe22.15,4X,1pe22.15,4X,1pe22.15)') Gamma(3,1,1,1),Gamma(3,2,1,1),Gamma(3,3,1,1)
        write(*,*)
        write(*,'(" RomberG - YZ-layer of the second hyperpolarizability tensor (γXXYZ,γXYYZ,γXZYZ):",4X,1pe22.15,4X,1pe22.15,4X,1pe22.15)') Gamma(1,1,2,2),Gamma(1,2,2,2),Gamma(1,3,2,2)
        write(*,'(" RomberG - YZ-layer of the second hyperpolarizability tensor (γYXYZ,γYYYZ,γYZYZ):",4X,1pe22.15,4X,1pe22.15,4X,1pe22.15)') Gamma(2,1,2,2),Gamma(2,2,2,2),Gamma(2,3,2,2)
        write(*,'(" RomberG - YZ-layer of the second hyperpolarizability tensor (γZXYZ,γZYYZ,γZZYZ):",4X,1pe22.15,4X,1pe22.15,4X,1pe22.15)') Gamma(3,1,2,2),Gamma(3,2,2,2),Gamma(3,3,2,2)
        write(*,*)
        write(*,'(" RomberG - ZZ-layer of the second hyperpolarizability tensor (γXXZZ,γXYZZ,γXZZZ):",4X,1pe22.15,4X,1pe22.15,4X,1pe22.15)') Gamma(1,1,3,3),Gamma(1,2,3,3),Gamma(1,3,3,3)
        write(*,'(" RomberG - ZZ-layer of the second hyperpolarizability tensor (γYXZZ,γYYZZ,γYZZZ):",4X,1pe22.15,4X,1pe22.15,4X,1pe22.15)') Gamma(2,1,3,3),Gamma(2,2,3,3),Gamma(2,3,3,3)
        write(*,'(" RomberG - ZZ-layer of the second hyperpolarizability tensor (γZXZZ,γZYZZ,γZZZZ):",4X,1pe22.15,4X,1pe22.15,4X,1pe22.15)') Gamma(3,1,3,3),Gamma(3,2,3,3),Gamma(3,3,3,3)
        write(*,*)

    end if
end if

    !-- Computing the isotropic values
if (doIsotropic.eqv..TRUE.) then
    if (onlop.eq.2) then
            !-- Compute the isotropic polarizability

            !-- Compute deltaAlpha: 10.1007/s42452-025-07291-9
        if (inlop.eq.0) write(*,'(" RomberG - WARNING! Δα requires crossed terms. Do not trust the printed value.")')
        deltaAlpha=dsqrt((Alpha(1,1)-Alpha(2,2))**2.0d0+(Alpha(2,2)-Alpha(3,3))**2.0d0+(Alpha(3,3)-Alpha(1,1))**2.0d0+12.0d0*(Alpha(1,2)**2.0d0+Alpha(2,3)**2.0d0+Alpha(3,1)**2.0d0))
        !deltaAlpha=dsqrt(Alpha(1,1)-Alpha(2,2)**2.0d0+(Alpha(2,2)-Alpha(3,3))**2.0d0+(Alpha(3,3)-Alpha(1,1))**2.0d0+24.0d0*(Alpha(1,2)**2.0d0+Alpha(2,3)**2.0d0+Alpha(3,1)**2.0d0))
        deltaAlpha=deltaAlpha*(dsqrt(2.0d0)/2.0d0)
        write(*,'(" RomberG - Δα = ",1pe22.15)') deltaAlpha

    else if (onlop.eq.3) then
            !-- Compute different values for the isotropic first hyperpolarizability
        
            !-- Get the dipole moment of the zeroth field to compute the Beta4
        call system("rm tmpDipole0.nlop")
        call system("cd $(pwd)/"//mol_name//" ; grep -A 1 'Dipole' *.fchk > tmpDipole0.nlop; sed -i '/Dipole/d' tmpDipole0.nlop; cp tmpDipole0.nlop ../")
        open (unit=44,file="tmpDipole0.nlop",status="old")
        read(44,*,iostat=ios) line,tmpDipole(1),tmpDipole(2),tmpDipole(3)
        tmpDipModulus=dsqrt(tmpDipole(1)**2.0d0+tmpDipole(2)**2.0d0+tmpDipole(3)**2.0d0)

            !-- Vectorial hyperpolarizability
        if (inlop.eq.0) write(*,'(" RomberG - WARNING! The vectorial hyperpolarizability is lacking terms to compute. Do not fully trust the printed value")')
        beta_vec=0.0d0
        do i=1,3
            do j=1,3
                beta_vec=beta_vec+Beta(i,j,j)
            end do
            beta_vec=beta_vec**2.0d0
        end do
        beta_vec=dsqrt(beta_vec)
        write(*,'(" RomberG - Vectorial hyperpolarizability = ",1pe22.15)') beta_vec

            !-- Parallel hyperpolarizability: J. Chem. Phys. 130, 194108 2009
        if (inlop.eq.0) write(*,'(" RomberG - WARNING! The parallel hyperpolarizability is lacking terms to compute. Do not fully trust the printed value")')
        beta_parallel=0.0d0
        do i=1,3
            do j=1,3
                 beta_parallel=beta_parallel+tmpDipole(i)*Beta(i,j,j)
            end do
        end do
        beta_parallel=3.0d0*beta_parallel/5.0d0/tmpDipModulus
        write(*,'(" RomberG - Parallel hyperpolarizability = ",1pe22.15)') beta_parallel

            !-- Perpendicular hyperpolarizability: J. Chem. Phys. 152, 244106 (2020)
        if (inlop.eq.0) write(*,'(" RomberG - WARNING! The perpendicular hyperpolarizability is lacking terms to compute. Do not fully trust the printed value")')
        beta_perpendicular=0.0d0
        do i=1,3
            do j=1,3
                beta_perpendicular=beta_perpendicular+(2.0d0*Beta(i,j,j)-3.0d0*Beta(j,i,j)+2.0d0*Beta(j,j,i))*tmpDipole(i)
            end do
        end do
        beta_perpendicular=beta_perpendicular/5.0d0
        write(*,'(" RomberG - Perpendicular hyperpolarizability = ",1pe22.15)') beta_perpendicular

            !-- "Beta4" from 10.1021/acs.jpca.5c00383; equivalent to parallel hyperpolarizability but without a factor
        beta_4=0.0d0
        do i=1,3
            do j=1,3
               beta_4=beta_4+Beta(i,j,j)*tmpDipole(i) 
            end do
        end do
        beta_4=beta_4/tmpDipModulus


    else if (onlop.eq.4) then

            !-- Compute the parallel second hyperpolarizability
        if (inlop.eq.0.or.inlop.eq.1) write(*,'(" RomberG - WARNING! The parallel second hyperpolarizability is lacking terms to compute. Do not fully trust the printed value")')
        gamma_parallel=0.0d0
        do i=1,3
            do j=1,3
                gamma_parallel=gamma_parallel+Gamma(i,i,j,j)
            end do
        end do    
        gamma_parallel=gamma_parallel/5.0d0
        write(*,'(" RomberG - Parallel second hyperpolarizability  = ",1pe22.15)') gamma_parallel

            !-- Compute the average second hyperpolarizability: 10.1007/s42452-025-07291-9
        if (inlop.eq.0.or.inlop.eq.1) write(*,'(" RomberG - WARNING! The average second hyperpolarizability is lacking terms to compute. Do not fully trust the printed value")')
        gamma_average=BestGamma(1,1)+BestGamma(5,1)+BestGamma(5,3)
        do i=1,2
            do j=2,3
                if (i.eq.j) cycle
                    !\gamma_average=gamma_average+2.0d0(\gamma_xxyy+\gamma_xxzz+\gamma_yyzz), but iijj==ijij=ijji=jjii=... -> \gamma_average=\gamma_average+2.0d0*((6.0d0/2.0d0)*\gamma_xxyy+...)
                gamma_average=gamma_average+6.0d0*Gamma(i,i,j,j) !\gamma_xxyy+\gamma_xxzz+\gamma_yyzz
            end do
        end do
        write(*,'(" RomberG - Average second hyperpolarizability  = ",1pe22.15)') gamma_average

    end if
end if
    !-- Closing units and removing .nlop files, they are still stored on the molecule directory
close(unit=4)
close(unit=44)
write(*,*)
write(*,*) "RomberG - Romberg procedure done!"
write(*,*)

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

Function BestProp(aP,bP,aE,bE)
implicit none
double precision :: BestProp,aP,bP,aE,bE

BestProp=aP
if (abs(aE).gt.abs(bE)) BestProp=bP

End function BestProp

Subroutine Alpha2Alpha(Best,Symmetric)
implicit none
double precision, intent(in), dimension(6) :: Best
double precision, intent(out), dimension(3,3) :: Symmetric

    !-- Diagonal elements
Symmetric(1,1)=Best(1) !\alpha_xx
Symmetric(2,2)=Best(3) !\alpha_yy
Symmetric(3,3)=Best(6) !\alpha_zz

    !-- Off-diagonal elements
Symmetric(1,2)=Best(2) !\alpha_xy
Symmetric(2,1)=Best(2) !\alpha_yx
Symmetric(1,3)=Best(4) !\alpha_xz
Symmetric(3,1)=Best(4) !\alpha_zx
Symmetric(2,3)=Best(5) !\alpha_yz
Symmetric(3,2)=Best(5) !\alpha_yz

End subroutine

Subroutine Beta2Beta(Best,Symmetric)
implicit none
double precision, intent(in), dimension(5,2) :: Best
double precision, intent(out), dimension(3,3,3) :: Symmetric

    !-- Diagonal elements
Symmetric(1,1,1)=Best(1,1) !\beta_xxx
Symmetric(2,2,2)=Best(4,1) !\beta_yyy
Symmetric(3,3,3)=Best(5,2) !\beta_zzz

    !-- Off-diagonal elements
Symmetric(1,1,2)=Best(2,1) !\beta_xxy
Symmetric(1,2,2)=Best(3,1) !\beta_xyy
Symmetric(1,1,3)=Best(5,1) !\beta_xxz
Symmetric(1,2,3)=Best(1,2) !\beta_xyz
Symmetric(2,2,3)=Best(2,2) !\beta_yyz
Symmetric(1,3,3)=Best(3,2) !\beta_xzz
Symmetric(2,3,3)=Best(4,2) !\beta_yzz

    !-- Transposition of the elements
                        !xxy
Symmetric(2,1,1)=Symmetric(1,1,2);Symmetric(1,2,1)=Symmetric(1,1,2)
                        !xyy
Symmetric(2,1,2)=Symmetric(1,2,2);Symmetric(2,2,1)=Symmetric(1,2,2)
                        !xxz
Symmetric(3,1,1)=Symmetric(1,1,3);Symmetric(1,3,1)=Symmetric(1,1,3)
                        !yyz
Symmetric(3,2,2)=Symmetric(2,2,3);Symmetric(2,3,2)=Symmetric(2,2,3)
                        !xzz
Symmetric(3,1,3)=Symmetric(1,3,3);Symmetric(3,3,1)=Symmetric(1,3,3)
                        !yzz
Symmetric(3,3,2)=Symmetric(2,3,3);Symmetric(3,2,3)=Symmetric(2,3,3)
                        !xyz
Symmetric(1,3,2)=Symmetric(1,2,3);Symmetric(2,1,3)=Symmetric(1,2,3)
Symmetric(2,3,1)=Symmetric(1,2,3);Symmetric(3,1,2)=Symmetric(1,2,3)
Symmetric(3,2,1)=Symmetric(1,2,3)

End subroutine

Subroutine Gamma2Gamma(Best,Symmetric)
implicit none
double precision, intent(in), dimension(5,3) :: Best
double precision, intent(out), dimension(3,3,3,3) :: Symmetric

    !-- Diagonal elements
Symmetric(1,1,1,1)=Best(1,1) !\gamma_xxxx
Symmetric(2,2,2,2)=Best(5,1) !\gamma_yyyy
Symmetric(3,3,3,3)=Best(5,3) !\gamma_zzzz

    !-- Off-diagonal elements
Symmetric(1,1,2,1)=Best(2,1) !\gamma_xxyx
Symmetric(1,1,2,2)=Best(3,1) !\gamma_xxyy
Symmetric(2,1,2,2)=Best(4,1) !\gamma_yxyy
Symmetric(1,1,3,1)=Best(1,2) !\gamma_xxzx
Symmetric(1,1,3,2)=Best(2,2) !\gamma_xxzy
Symmetric(2,1,3,2)=Best(3,2) !\gamma_yxzy
Symmetric(2,2,3,2)=Best(4,2) !\gamma_yyzy
Symmetric(1,1,3,3)=Best(5,2) !\gamma_xxzz
Symmetric(2,1,3,3)=Best(1,3) !\gamma_yxzz
Symmetric(2,2,3,3)=Best(2,3) !\gamma_yyzz
Symmetric(3,1,3,3)=Best(3,3) !\gamma_zxzz
Symmetric(3,2,3,3)=Best(4,3) !\gamma_zyzz

    !-- Transposition of the elements
                        !xxyx
Symmetric(2,1,1,1)=Symmetric(1,1,2,1);Symmetric(1,2,1,1)=Symmetric(1,1,2,1);Symmetric(1,1,1,2)=Symmetric(1,1,2,1)
                        !yxyy
Symmetric(1,2,2,2)=Symmetric(2,1,2,2);Symmetric(2,2,1,2)=Symmetric(2,1,2,2);Symmetric(2,2,2,1)=Symmetric(2,2,1,2)
                        !xxzx
Symmetric(3,1,1,1)=Symmetric(1,1,3,1);Symmetric(1,3,1,1)=Symmetric(1,1,3,1);Symmetric(1,1,1,3)=Symmetric(1,3,1,1)
                        !yyzy
Symmetric(3,2,2,2)=Symmetric(2,2,3,2);Symmetric(2,3,2,2)=Symmetric(2,2,3,2);Symmetric(2,2,2,3)=Symmetric(2,2,3,2)
                        !zxzz
Symmetric(1,3,3,3)=Symmetric(3,1,3,3);Symmetric(3,3,1,3)=Symmetric(3,1,3,3);Symmetric(3,3,3,1)=Symmetric(3,3,1,3)
                        !zyzz
Symmetric(2,3,3,3)=Symmetric(3,2,3,3);Symmetric(3,3,2,3)=Symmetric(3,2,3,3);Symmetric(3,3,3,2)=Symmetric(3,3,2,3)
                        !xxyy
Symmetric(1,2,1,2)=Symmetric(1,1,2,2);Symmetric(1,2,2,1)=Symmetric(1,1,2,2);Symmetric(2,1,1,2)=Symmetric(1,1,2,2)
Symmetric(2,1,2,1)=Symmetric(1,1,2,2);Symmetric(2,2,1,1)=Symmetric(1,1,2,2)
                        !xxzz
Symmetric(1,3,1,3)=Symmetric(1,1,3,3);Symmetric(1,3,3,1)=Symmetric(1,1,3,3);Symmetric(3,1,1,3)=Symmetric(1,1,3,3)
Symmetric(3,1,3,1)=Symmetric(1,1,3,3);Symmetric(3,3,1,1)=Symmetric(1,1,3,3)
                        !yyzz
Symmetric(2,3,2,3)=Symmetric(2,2,3,3);Symmetric(2,3,3,2)=Symmetric(2,2,3,3);Symmetric(3,2,2,3)=Symmetric(2,2,3,3)
Symmetric(3,1,3,1)=Symmetric(2,2,3,3);Symmetric(3,3,2,2)=Symmetric(2,2,3,3)
                        !xxzy
Symmetric(1,1,2,3)=Symmetric(1,1,3,2);Symmetric(1,2,1,3)=Symmetric(1,1,3,2);Symmetric(1,2,3,1)=Symmetric(1,1,3,2)
Symmetric(1,3,1,2)=Symmetric(1,1,3,2);Symmetric(1,3,2,1)=Symmetric(1,1,3,2);Symmetric(2,1,1,3)=Symmetric(1,1,3,2)
Symmetric(2,1,3,1)=Symmetric(1,1,3,2);Symmetric(3,1,1,2)=Symmetric(1,1,3,2);Symmetric(3,1,2,1)=Symmetric(1,1,3,2)
Symmetric(3,2,1,1)=Symmetric(1,1,3,2);Symmetric(2,3,1,1)=Symmetric(1,1,3,2)
                        !yxzy
Symmetric(1,2,2,3)=Symmetric(2,1,3,2);Symmetric(1,2,3,2)=Symmetric(2,1,3,2);Symmetric(1,3,2,2)=Symmetric(2,1,3,2)
Symmetric(2,2,3,1)=Symmetric(2,1,3,2);Symmetric(2,3,1,2)=Symmetric(2,1,3,2);Symmetric(2,3,2,1)=Symmetric(2,1,3,2)
Symmetric(3,1,2,2)=Symmetric(2,1,3,2);Symmetric(3,2,1,2)=Symmetric(2,1,3,2);Symmetric(3,2,2,1)=Symmetric(2,1,3,2)
Symmetric(2,1,2,3)=Symmetric(2,1,3,2);Symmetric(2,2,1,3)=Symmetric(2,1,3,2)
                        !yxzz
Symmetric(1,2,3,3)=Symmetric(2,1,3,3);Symmetric(1,3,2,3)=Symmetric(2,1,3,3);Symmetric(1,3,3,2)=Symmetric(2,1,3,3)
Symmetric(2,3,1,3)=Symmetric(2,1,3,3);Symmetric(2,3,3,1)=Symmetric(2,1,3,3);Symmetric(3,1,2,3)=Symmetric(2,1,3,3)
Symmetric(3,1,3,2)=Symmetric(2,1,3,3);Symmetric(3,2,1,3)=Symmetric(2,1,3,3);Symmetric(3,2,3,1)=Symmetric(2,1,3,3)
Symmetric(3,3,1,2)=Symmetric(2,1,3,3);Symmetric(3,3,2,1)=Symmetric(2,1,3,3)

End subroutine
