PROGRAM test

    use healpix_modules
    implicit none
    real(kind=8),dimension(3) :: vector
    integer :: ipring
    print*,' pi = ',PI
    print*,' number of pixels in a Nside=8 map:',nside2npix(8)

    vector(1) = 13*3; vector(2) = 4*3; vector(3) = -5*3
    call vec2pix_ring(8,vector,ipring)
    print*,' vector=',vector
    print*,' ipring=',ipring
END PROGRAM test