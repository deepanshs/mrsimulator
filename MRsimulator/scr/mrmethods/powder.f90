

subroutine lineshape(spec, m, nt, fstart, fwidth, aniso, eta, iso, &
                     n, xr, yr, zr, powamp)!, npros)

!use omp_lib
implicit none

integer*4, intent(in) :: m, nt, n!, npros
!f2py integer*4, intent(in) :: m, nt, n!, npros

double precision, intent(in) :: fstart, fwidth
!f2py double precision, intent(in) :: fstart, fwidth

double precision, dimension(0:nt, 0:2*nt), intent(in) :: xr, yr, zr, powamp
!f2py double precision, intent(in) :: xr, yr, zr, powamp

double precision, intent(in), dimension(0:n-1) :: aniso, eta, iso
!f2py double precision, intent(in) :: aniso, eta, iso

double precision, intent(out), dimension(0:m-1, 0:n-1) :: spec
!f2py double precision, intent(out) :: spec

integer*4 :: i,j, arraySize
double precision, dimension(0:nt, 0:2*nt) :: powfreq !xr, yr, zr, powamp,
double precision, dimension(0:m-1) :: temp1, temp2!, temp3, temp4
double precision :: finc, sxx, syy, szz, tempAniso, tempEta, tempIso, amp1, amp2



!call SetupPowder2(nt, xr, yr, zr, powamp)

finc = fwidth/m
!spec = 0.0

do arraySize = 0, n-1
    tempAniso = aniso(arraySize)
    tempEta = eta(arraySize)
    tempIso = iso(arraySize)

!! old
    sxx = tempIso - 0.5*tempAniso* (1. + tempEta)
    syy = tempIso - 0.5*tempAniso* (1. - tempEta)
    szz = tempIso + tempAniso

    powfreq(:,:) = sxx*xr(:,:) + syy*yr(:,:) + szz*zr(:,:)

    temp1 = 0.0
    temp2 = 0.0


    !/* Interpolate between frequecies by setting up tents */

    !$omp parallel sections num_threads(2) shared(nt, powamp, powfreq, &
    !$omp          m, fstart, finc, temp1, temp2) &
    !$omp private(amp1, amp2, i, j)

    !$omp section

    !$omp parallel num_threads(2) shared(nt, powamp, powfreq, &
    !$omp          m, fstart, finc, temp1, temp2) &
    !$omp private(amp1, amp2, i, j)
    !$omp do
    do i=0, nt-1
        do j=0, nt-1
            amp1 = powamp(i+1,j)+powamp(i,j+1)+powamp(i,j)
            amp2 = powamp(i+1,j)+powamp(i,j+1)+powamp(i+1,j+1)
            call tent(powfreq(i+1,j), powfreq(i,j+1), powfreq(i,j), &
                        amp1, temp1, m, fstart, finc)
            call tent(powfreq(i+1,j), powfreq(i,j+1), powfreq(i+1,j+1), &
                        amp2, temp1, m, fstart, finc)
        enddo
    enddo
    !$omp end do
    !$omp end parallel

    !$omp section

    !$omp parallel num_threads(2) shared(nt, powamp, powfreq, &
    !$omp          m, fstart, finc, temp1, temp2) &
    !$omp private(amp1, amp2, i, j)
    !$omp do
    do i=0, nt-1
        do j=nt, 2*nt-1
            amp1 = powamp(i,j)+powamp(i+1,j+1)+powamp(i+1,j)
            amp2 = powamp(i,j)+powamp(i+1,j+1)+powamp(i,j+1)
            call tent(powfreq(i,j), powfreq(i+1,j+1), powfreq(i+1,j), &
                        amp1, temp2, m, fstart, finc);
            call tent(powfreq(i,j), powfreq(i+1,j+1), powfreq(i,j+1), &
                        amp2, temp2, m, fstart, finc);
        enddo
    enddo
    !$omp end do
    !$omp end parallel

    !$omp end parallel sections
    spec(:, arraySize) = temp1 + temp2! + temp3 + temp4
enddo


end subroutine lineshape








subroutine SetupPowder2(nt, xr, yr, zr, rrr)
implicit none

integer*4, intent(in) :: nt
double precision, dimension(0:nt, 0:2*nt), intent(out) :: xr, yr, zr, rrr

integer*4 :: i,j
double precision :: x,y,z, x2,y2,z2,r2

xr = 0.0
yr = 0.0
zr = 0.0
rrr = 0.0

!/* Do the (x + y + z = nt) face of the octahedron
!z -> 0 to nt-1
!y -> 0 to nt-z
!x -> nt - y - z
!*/

do j=0, nt-1
    do i=0, nt-j
        x = nt-i-j
        y = i
        z = j
        x2 = x*x
        y2 = y*y
        z2 = z*z
        r2 = x2 + y2 + z2
        xr(i,j) = x2/r2
        yr(i,j) = y2/r2
        zr(i,j) = z2/r2
        rrr(i,j) = 1.0/(r2*sqrt(r2))
    enddo
enddo

!/* Do the (-x + y + z = nt) face of the octahedron
!z -> 0 to nt-1
!y -> 0 to nt-z
!x -> -nt + y + z
!*/

do j=0, nt-1
    do i=nt-j+1, nt
        x = nt-i-j
        y = nt-j
        z = nt-i
        x2 = x*x
        y2 = y*y
        z2 = z*z
        r2 = x2 + y2 + z2
        xr(i,j) = x2/r2
        yr(i,j) = y2/r2
        zr(i,j) = z2/r2
        rrr(i,j) = 1.0/(r2*sqrt(r2))
    enddo
enddo

!/* Do the (-x - y + z = nt) face of the octahedron
!*/

do j=nt, 2*nt-1
    do i=j-nt+1, nt-1
        x = -nt-i+j
        y = nt-j
        z = nt-i
        x2 = x*x
        y2 = y*y
        z2 = z*z
        r2 = x2 + y2 + z2
        xr(i,j) = x2/r2
        yr(i,j) = y2/r2
        zr(i,j) = z2/r2
        rrr(i,j) = 1.0/(r2*sqrt(r2))
    enddo
enddo

!/* Do the (x - y + z = nt) face of the octahedron
!*/

do j=nt, 2*nt-1
    do i=1, j-nt
        x = -nt-i+j
        y = -i
        z = 2*nt-j
        x2 = x*x
        y2 = y*y
        z2 = z*z
        r2 = x2 + y2 + z2
        xr(i,j) = x2/r2
        yr(i,j) = y2/r2
        zr(i,j) = z2/r2
        rrr(i,j) = 1.0/(r2*sqrt(r2))
    enddo
enddo

xr(0,nt) = 0.0
yr(0,nt) = 0.0
zr(0,nt) = 1.0
r2 = nt
rrr(0,nt) = 1.0/(r2*r2*r2)

do j=0, nt-1
    i = 2*nt-j
    xr(0,i) = xr(0,j)
    yr(0,i) = yr(0,j)
    zr(0,i) = zr(0,j)
    rrr(0,i) = rrr(0,j)
enddo

do i=0, nt
    xr(nt,nt+i) = xr(i,0)
    yr(nt,nt+i) = yr(i,0)
    zr(nt,nt+i) = zr(i,0)
    rrr(nt,nt+i) = rrr(i,0)
enddo


i = 2*nt
do j=1, nt-1
    xr(nt-j,i) = xr(nt,j)
    yr(nt-j,i) = yr(nt,j)
    zr(nt-j,i) = zr(nt,j)
    rrr(nt-j,i) = rrr(nt,j)
enddo

end subroutine SetupPowder2



subroutine tent(freq1, freq2, freq3, amp, spec, points, fstart, finc)

implicit none
integer*4, intent(in) :: points
double precision, intent(in) :: freq1, freq2, freq3, amp, fstart, finc
double precision, dimension(0:points-1), intent(inout) :: spec

double precision :: df1,df2,f1,f2,top, t
double precision, dimension(0:2) :: f

integer*4 :: p,pmid,pmax, i,j

f(0) = freq1
f(1) = freq2
f(2) = freq3
do j=1, 2
    t = f(j);
    i=j-1;
    do while(i .ge. 0 .and. f(i) .gt. t)
        f(i+1) = f(i)
        i = i-1
    end do
    f(i+1)=t
enddo

top = amp*2./(f(2)-f(0))
p = floor((f(0)-fstart)/finc)
pmid = floor((f(1)-fstart)/finc)
pmax = floor((f(2)-fstart)/finc)
df1 = 0.
df2 = 0.
if((f(1)-f(0)) .ne. 0.) df1 = top/(2*(f(1)-f(0)))
if((f(2)-f(1)) .ne. 0.) df2 = top/(2*(f(2)-f(1)))
if((pmax .lt. points) .and. (p .ge. 0)) then
    if(p .ne. pmid) then
        f2 = finc * ( dble(p) +1.) + fstart;
        spec(p) = spec(p) + (f2-f(0))*(f2-f(0))*df1
        p = p+1
        f1 = f2
        do while(p .ne. pmid)
            f2 = finc * ( dble(p) + 1.) + fstart
            spec(p) = spec(p) + finc * ( (f2-f(0)) + (f1-f(0)) ) * df1
            p = p+1
            f1 = f2
        enddo
        spec(p) = spec(p) + (f(1)-f1)*((f(1)-f(0))+(f1-f(0)))*df1
    else
        spec(p) = spec(p) + (f(1)-f(0))*top/2
    endif

    if(p .ne. pmax) then
        f2 = finc * ( dble(pmid) + 1.) + fstart
        spec(p) = spec(p) + (f2-f(1))*((f(2)-f2)+(f(2)-f(1)))*df2
        p = p+1
        f1 = f2
        do while(p .ne. pmax)
            f2 = finc*( dble(p) + 1.)+fstart
            spec(p) = spec(p) + finc*((f(2)-f1)+(f(2)-f2))*df2
            p = p+1
            f1 = f2
        enddo
        spec(p) = spec(p) + (f(2)-f1)*(f(2)-f1)*df2
    else
        spec(p) = spec(p) + (f(2)-f(1))*top/2
    endif
endif

end subroutine tent
