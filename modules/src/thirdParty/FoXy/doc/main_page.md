project: FoXy
project_dir: ./src
output_dir: ./doc/html/publish/
exclude_dir: third_party
project_github: https://github.com/Fortran-FOSS-Programmers/FoXy
summary: Fortran XML parser for poor people
author: Fortran-FOSS-Programmers Group
github: http://fortran-foss-programmers.github.io/
website: http://fortran-foss-programmers.github.io/
! md_extensions: markdown.extensions.toc(anchorlink=True)
md_extensions: markdown.extensions.toc()
               markdown.extensions.smarty(smart_quotes=False)
               markdown.extensions.extra
               markdown_checklist.extension
docmark: <
display: public
         protected
         private
source: true
warn: true
graph: true
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html

{!README-FoXy.md!}
