program volatile_doctest
use befor64
 use befor64
 use penf
 character(len=:), allocatable :: code64
 call b64_encode_up(up=1._R8P, code=code64)
 print "(A)", code64
endprogram volatile_doctest