	program calr
	use foxy
	implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	config file variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=50) :: configFile,tmpFld
    character(len=:), allocatable :: parsed         !< String containing the parsed XML data.
    character(len=:), allocatable :: tagcontent
    type(xml_file)                :: xfile          !< XML file handler.
    integer :: fileid,ntags,i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	real*8 :: x0,y0,x1,y1,x2,y2,x3,y3,r,t,a
	integer :: inc,ref,tras

    configFile="config.xml"
    call xfile%parse(filename=configFile)

    tmpFld="r"//trim(adjustl(xfile%content("inc")))//"Syin.dat"
    open(unit=1,file=tmpFld, status="old",action="read")
    tmpFld="r"//trim(adjustl(xfile%content("ref")))//"Syr.dat"
    open(unit=2,file=tmpFld, status="old",action="read")
    tmpFld="r"//trim(adjustl(xfile%content("tras")))//"Syt.dat"
    open(unit=3,file=tmpFld, status="old",action="read")

    open(unit=5,file="RTA_2Dfdtd.data", status="unknown")
	do
	    read(1,*,end=12)x0,y0
	    read(2,*)x1,y1
	    read(3,*) x2,y2
	    r=-y1/y0
	    t=y2/y0
	    a=1.d0-r-t
	    write(5,*)x0,r,t,a
	end do
 12	print *,"end calculating RTA"
	end program calr


