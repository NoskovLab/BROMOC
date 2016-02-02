!    BROMOC  -  CG-GCMC-BD
!    Electrodiffusion, Gran Canonical Monte Carlo, Brownian,Dynamics 
!    and Coarse Grain Model DNA Simulation Program.
!    Copyright (C) 2014 Pablo M. De Biase (pablodebiase@gmail.com)
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.

logical*1 function getopt (args, largs, nargs, flgs, lflgs, nflgs, opts, nopts, optarg)
! 
!.....getopt - parses the command line for arguments, flags and options.
!
!     The routine returns .false. if an error occurs, .true. otherwise.
!
!     Argument is any token in the commmand line that doesn't begin with
!     '-'. A standalone '-' will be considered a special argument, used
!     to skip an arg without filling it.
!
!     Flag is any char of a token that starts with a '-'.
!
!     Option is a flag that accepts an argument, which must be the token
!     following the option.
!
!     A '--' in the input line signals the end of options and flags (in
!     this way, arguments beginning with '-' are allowed).
!
!..Input parameters:
!
!  nargs ......... maximum number of allowed arguments.
!  nflgs ......... number of valid flags.
!  flgs .......... list of valid flags.
!  nopts ......... number of valid options.
!  opts .......... list of valid options.
!
!..Output parameters:
!
!  args() ........ returns the values of the specified arguments.
!  largs() ....... returns .true. on each specified argument (that is,
!                  '-' will increase the count of arguments but with
!                  a .false. in lflags).
!  lflgs() ....... returns a .true. in the position of the flags
!                  actually found in the command line.
!  lopts() ....... returns .true. in the position of used options.
!  optarg() ...... contains the arguments associated to the options
!                  actually found.
!
!..Example:
!
!  Given a program that expects for the next entrance:
!
!       myprog [input] [output] [-abcde] [-o arg1] [-s arg2]
!
!  where square braces [] contain optional arguments, the main program
!  may use the next call to getopt to parse the command line:
!
!       nargs = 2
!       nflgs = 5
!       flgs = 'abcde'
!       nopts = 2
!       opts = 'os'
!       if (.not. getopt (args, largs, nargs, flgs, lflgs, nflgs,
!      &            opts, lopts, nopts, optarg)) then
!          --- Something went wrong ---
!       endif
!
! PARAMETERS:
!
implicit none
integer         nargs, nflgs, nopts
logical*1           largs(*), lflgs(*)
character*(*)     args(*), optarg(*)
character*(*)     flgs, opts

! LOCAL VARIABLES:

integer         iargc
integer         iargs, iflgs, iopts
integer         ntokens, itoken, ltoken, lp
character*(128)   token
logical*1           endopts

!.....get the number of tokens in the command line:

iargs = 0
ntokens = iargc ()
if (ntokens.lt.0) then
  getopt = .false.
  return
endif

!.....loop over tokens:

endopts = .false.
itoken = 1
do while (itoken.le.ntokens)
  call getarg (itoken, token)
  itoken = itoken + 1
  ltoken = len_trim(token)

!.....is it an arg? (or we found the end of options)

  if (token(1:1).ne.'-' .or. endopts) then
    if (iargs.lt.nargs) then
      iargs = iargs + 1
      args(iargs) = token(1:ltoken)
      largs(iargs) = .true.
    else
      write (0,*) 'getopt: too many arguments'
      getopt = .false.
      return
    endif

!.....dummy argument:

  else if (ltoken.eq.1) then
    if (iargs.lt.nargs) then
      iargs = iargs + 1
      args(iargs) = token(1:ltoken)
      largs(iargs) = .false.
    else
      write (0,*) 'getopt: too many arguments'
      getopt = .false.
      return
    endif

!.....end of options:

  else if (token(1:ltoken).eq.'--') then
    endopts = .true.

!.....everything else are flags and options:

  else

!.....loop for flags:

    do lp = 2, ltoken-1
      iflgs = 1
      do while (iflgs.lt.nflgs .and. token(lp:lp).ne.flgs(iflgs:iflgs))
        iflgs = iflgs + 1
      enddo
      if (iflgs.eq.nflgs .and. token(lp:lp).ne.flgs(iflgs:iflgs)) then
        write (0,*) 'getopt: -'//token(lp:lp)//' is not a valid flag'
        getopt = .false.
        return
      endif
      lflgs(iflgs) = .true.
    enddo

!.....last one can be a flag or an option:

    lp = ltoken
    iflgs = 1
    do while (iflgs.lt.nflgs .and. token(lp:lp).ne.flgs(iflgs:iflgs))
      iflgs = iflgs + 1
    enddo

!.....it is a flag:

    if (iflgs.lt.nflgs .or. token(lp:lp).eq.flgs(iflgs:iflgs)) then
      lflgs(iflgs) = .true.

!.....not a flag:

    else
      iopts = 1
      do while (iopts.lt.nopts .and. token(lp:lp).ne.opts(iopts:iopts))
        iopts = iopts + 1
      enddo

!.....it is an option:

      if (iopts.ne.nopts .or. token(lp:lp).eq.opts(iopts:iopts)) then
        if (itoken.le.ntokens) then
          call getarg (itoken, token)
          itoken = itoken + 1
          ltoken = len_trim(token)
          optarg(iopts) = token(1:ltoken)
        else
          write (0,*) 'getopt: missing argument in option -'//token(lp:lp)
          getopt = .false.
          return
        endif

!.....not an option:

      else
        write (0,*) 'getopt: -'//token(lp:lp)//' is not a valid flag'
        getopt = .false.
        return
      endif
    endif
  endif
enddo

!.....end of parsing ok:

getopt = .true.
return
end function
