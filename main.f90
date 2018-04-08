program main
use tools
use types
use functions

    implicit none
    real (kind=rk), dimension (:,:), allocatable :: a
    real (kind=rk), dimension (:), allocatable :: x
    real (kind=rk), dimension (:), allocatable :: b

    integer, parameter :: n=4

    allocate ( a(n,n) )
    allocate ( x(n) )
    allocate ( b(n) )

    a = transpose(reshape((/1.0,  6.0, -5.0, 2.0, &
                            2.0, -3.0,  12.0, 4.0, &
                            2.0, -6.0,  1.0,  7.0, &
                            4.0,  4.0, -4.0,  4.0  /), &
                            (/4,4/)))

    b = (/1.0,  5.0,  1.0,  5.0/)

    call gauss(a,b,x)

    print*,"Solution set for Ax=b is:",x

end program main
