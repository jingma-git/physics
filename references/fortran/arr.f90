program arr
    implicit none
    integer :: arr1(4)
    real :: arr2(3, 3)
    integer, allocatable :: arr3(:)
    real, allocatable :: arr4(:, :)

    character(len=4) :: first_name
    character(len=5) :: last_name
    character(10) :: full_name
    
    arr1 = [1, 2, 3, 4]
    print *, arr1
    arr1(1:2) = 0
    arr1(3:) = 1
    print *, arr1
    arr2(:, :) = 0
    print *, arr2

    allocate(arr3(10))
    allocate(arr4(10, 10))

    deallocate(arr3)
    deallocate(arr4)


    first_name = 'John'
    last_name = 'Smith'
    full_name = first_name//' '//last_name
    print *, full_name
end program arr