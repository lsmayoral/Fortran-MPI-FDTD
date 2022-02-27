! post-processig fortran program to calculate the reflectance
!config file contains the number of node which contains the corresponding data,
!allowing the program to be compiled only once.
!When making any change which produce that the incident/reflected field move to differet node, 
!just change the config.xlm accordingly.

! in the following config.xlm example, the node 0 contains the data of the incident field
! and the node 2 contains the data of the reflected field.
!<inc> 0 </inc> 
!<ref> 2 </ref>
!
!
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
	real*8 :: x0,y0,I0,x1,y1,I1,r,rI
	integer :: inc,ref,tras

    configFile="config.xml"
    call xfile%parse(filename=configFile)

    tmpFld="r"//trim(adjustl(xfile%content("inc")))//"pump_Syin.dat"
    open(unit=1,file=tmpFld, status="old",action="read")
    tmpFld="r"//trim(adjustl(xfile%content("ref")))//"R_Syr.dat"
    open(unit=2,file=tmpFld, status="old",action="read")

    open(unit=5,file="R_1Dfdtd.data", status="unknown")
	do
	    read(1,*,end=12)x0,y0
	    read(2,*)x1,y1
	    r=-y1/y0
	    write(5,*)x0,r
	end do
 12	print *,"end calculating R"
	end program calr

