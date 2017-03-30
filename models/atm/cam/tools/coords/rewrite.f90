real a(49)
open(1,file='hya.48')
open(2,file='hyb.48')
do i=1,49
read(1,*)a(i)
end do
write(*,'(6f15.10)') a

print*

do i=1,49
read(2,*)a(i)
end do
write(*,'(6f15.10)') a

end
