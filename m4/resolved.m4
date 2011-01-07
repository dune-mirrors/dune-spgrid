AC_DEFUN([DUNE_RESOLVED_ABS_TOPSRCDIR],[
  ac_resolved_top_srcdir="$(pwd -P)"
  AC_SUBST([resolved_top_srcdir],[$ac_resolved_top_srcdir])
])
