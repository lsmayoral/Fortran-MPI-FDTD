!< StringiFor `split` test.
program split
!< StringiFor `split` test.
use, intrinsic :: iso_fortran_env, only : stdout => output_unit
use stringifor, only : string

implicit none
type(string)              :: astring         !< A string.
type(string), allocatable :: strings(:)      !< A set of strings.
logical                   :: test_passed(12) !< List of passed tests.
integer                   :: s               !< Counter.

test_passed = .false.

astring = '+ab-++cre-++cre-ab+'
write(stdout, "(A)") 'Split "'//astring//'" at "+"'
call astring%split(tokens=strings, sep='+')
test_passed(1) = (strings(1)//''=='ab-'.and.strings(2)//''=='cre-'.and.strings(3)//''=='cre-ab')
do s=1, size(strings)
  write(stdout, "(A)") '+ "'//strings(s)//'"'
enddo

astring = 'ab-++cre-++cre-ab+'
write(stdout, "(A)") 'Split "'//astring//'" at "+"'
call astring%split(tokens=strings, sep='+')
test_passed(2) = (strings(1)//''=='ab-'.and.strings(2)//''=='cre-'.and.strings(3)//''=='cre-ab')
do s=1, size(strings)
  write(stdout, "(A)") '+ "'//strings(s)//'"'
enddo

astring = 'ab-++cre-++cre-ab'
write(stdout, "(A)") 'Split "'//astring//'" at "+"'
call astring%split(tokens=strings, sep='+')
test_passed(3) = (strings(1)//''=='ab-'.and.strings(2)//''=='cre-'.and.strings(3)//''=='cre-ab')
do s=1, size(strings)
  write(stdout, "(A)") '+ "'//strings(s)//'"'
enddo

astring = 'Hello '//new_line('a')//'World!'
write(stdout, "(A)") 'Split "'//astring//'" at "new_line"'
call astring%split(tokens=strings, sep=new_line('a'))
test_passed(4) = (strings(1)//''=='Hello '.and.strings(2)//''=='World!')
do s=1, size(strings)
  write(stdout, "(A)") '+ "'//strings(s)//'"'
enddo

astring = 'Hello World!'
write(stdout, "(A)") 'Split "'//astring//'" at "default" (namely space)'
call astring%split(tokens=strings)
test_passed(5) = (strings(1)//''=='Hello'.and.strings(2)//''=='World!')
do s=1, size(strings)
  write(stdout, "(A)") '+ "'//strings(s)//'"'
enddo

astring = '+ab-'
write(stdout, "(A)") 'Split "'//astring//'" at "+"'
call astring%split(tokens=strings, sep='+')
test_passed(6) = (strings(1)//''=='ab-')
do s=1, size(strings)
  write(stdout, "(A)") '+ "'//strings(s)//'"'
enddo

astring = '+ab-'
write(stdout, "(A)") 'Split "'//astring//'" at "-"'
call astring%split(tokens=strings, sep='-')
test_passed(7) = (strings(1)//''=='+ab')
do s=1, size(strings)
  write(stdout, "(A)") '+ "'//strings(s)//'"'
enddo

astring = '+ab-+cd-'
write(stdout, "(A)") 'Split "'//astring//'" at "+"'
call astring%split(tokens=strings, sep='+')
test_passed(8) = (strings(1)//''=='ab-'.and.strings(2)//''=='cd-')
do s=1, size(strings)
  write(stdout, "(A)") '+ "'//strings(s)//'"'
enddo

astring = 'ab-+cd-+'
write(stdout, "(A)") 'Split "'//astring//'" at "+"'
call astring%split(tokens=strings, sep='+')
test_passed(9) = (strings(1)//''=='ab-'.and.strings(2)//''=='cd-')
do s=1, size(strings)
  write(stdout, "(A)") '+ "'//strings(s)//'"'
enddo

astring = '+ab-+cd-+'
write(stdout, "(A)") 'Split "'//astring//'" at "+"'
call astring%split(tokens=strings, sep='+')
test_passed(10) = (strings(1)//''=='ab-'.and.strings(2)//''=='cd-')
do s=1, size(strings)
  write(stdout, "(A)") '+ "'//strings(s)//'"'
enddo

astring = '1-2-3-4-5-6-7-8'
write(stdout, "(A)") 'Split "'//astring//'" at "-" in max 3 or 4 tokens'
call astring%split(tokens=strings, sep='-', max_tokens=3)
test_passed(11) = (strings(1)//''=='1'.and.strings(2)//''=='2'.and.strings(3)//''=='3'.and.strings(4)//''=='4-5-6-7-8')
do s=1, size(strings)
  write(stdout, "(A)") '+ "'//strings(s)//'"'
enddo

astring = '-1-2-3-4-5-6-7-8-'
write(stdout, "(A)") 'Split "'//astring//'" at "-" in chunks of 3'
call astring%split_chunked(tokens=strings, sep='-', chunks=3)
test_passed(12) = (strings(1)//''=='1'.and.strings(2)//''=='2'.and.strings(3)//''=='3'.and.strings(4)//''=='4'.and. &
                   strings(5)//''=='5'.and.strings(6)//''=='6'.and.strings(7)//''=='7'.and.strings(8)//''=='8')
do s=1, size(strings)
  write(stdout, "(A)") '+ "'//strings(s)//'"'
enddo

write(stdout, "(A,L1)") new_line('a')//'Are all tests passed? ', all(test_passed)
endprogram split
