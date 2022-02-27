!< StringiFor `unique` test.
program unique
!-----------------------------------------------------------------------------------------------------------------------------------
!< StringiFor `unique` test.
!-----------------------------------------------------------------------------------------------------------------------------------
use, intrinsic :: iso_fortran_env, only : stdout => output_unit
use stringifor, only : string
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
type(string) :: astring        !< A string.
logical      :: test_passed(1) !< List of passed tests.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
test_passed = .false.

astring = '+++ab-++cre-++cre-ab+++++'
write(stdout, "(A)") 'Original: "'//astring//'"'
test_passed(1) = astring%unique(substring='+')//''=='+ab-+cre-+cre-ab+'
write(stdout, "(A,L1)") 'Unique:   "'//astring%unique(substring='+')//'", is correct? ', test_passed(1)

write(stdout, "(A,L1)") new_line('a')//'Are all tests passed? ', all(test_passed)
stop
!-----------------------------------------------------------------------------------------------------------------------------------
endprogram unique
