program variables
    implicit none ! all variables will be explicitly declared
    integer :: amount
    real :: pi
    logical :: isOkay
    integer :: age
    real :: radius
    real :: height
    real :: area
    real :: volume
    real(kind=8) num_double
    real(kind=4) num_float
    num_double = 1/3.
    num_float = 1/3.
    print *, 'double', num_double
    print *, 'float', num_float

    amount = 10
    pi = 3.1415926
    isOkay = .false.
    print *, 'amount: ', amount
    print *, 'pi: ', pi
    print *, 'logical: ', isOkay
    
    ! print *, 'Please enter your age: '
    ! read(*,*) age
    ! print *, 'Your age is: ', age

    radius = 1
    height = 10
    area = pi * radius**2
    volume = area * height
    print *, 'area: ', area
    print *, 'vol: ', volume
end program variables