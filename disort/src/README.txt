FORTRANtoIDL -- Version 1.0 August 2001

This library of functions gives some examples of dynamicly linked modules
(DLMs) which call Fortran functions/subroutines.  While there are likely other
approaches to do this, I have adoped writing function 'wrapping' routines in C
which do type checking and communicate with IDL.  The Fortran routines (mainly
F77 in these examples) are then called from within the C wrapper.

I wrote most of these functions early in 2001 when I too was trying to use
Fortran functions from within IDL.  I experiemented with the call_external
routine but found it lacked the flexibility that I needed to use many of the
Fortran routines I have inhereted over the years.  This library is not intended
as an introduction to DLM programming.  If you are new to DLM programming, I
suggest you find a copy of Ronn Kling's book, "IDL and C: Using DLM's to extend
your IDL code" which is available at http://www.kilvarock.com.  Likewise, this
is not a really designed as a tutorial in writing mixed C-Fortran functions.
Calling Fortran functions from C is discussed in Chapter 11 of the "Fortran
Programmer's Guide" which is part of the Sun Workshop documentation (available
at http://docs.sun.com).  I have included this chapter (as a PDF file) in the
distribution because finding it on the reorganized Sun documentation site has
proved rather difficult.

Be forwarned that communicating between IDL and C is actually the easy part;
communicating between C and Fortran can lead to some realy nasty bugs which are
potentially OS dependent and may hang/crash some systems! There are undoubtably
inconsistancies, errors and ommissions in these routines.  I have tried to
ensure the code is portable, simple, clean and bug-free.  If you do find bugs
or feel the documentation is lacking, write me and I will make the necessary
changes.  I have tested this code on a variety of unix systems and with a
variety of compilers including:

   IDL 5.3/5.4 Linux GNU compiler Collection
   IDL 5.3/5.4 Linux GNU gcc + NAGWare Fortran 95 Compiler
   IDL 5.3/5.4 Linux GNU gcc + Compaq/Fujitsu Fortran 95 Compiler
   IDL 5.3     DEC/Compaq Development Suite
   IDL 5.2     Sun Workshop Compiler Suite

I do not have access to the Windows/MacOS version of IDL (or the required
compilers) so I have not tested these routines on either of these platforms.
There should not be any real difficulties in compiling for either of these
systems.  If you do have problems with either of these platforms, please let me
know.

I am freely distributing this code in the hope that it will be of use to those
who, like me, learn by example.  If you do find this code helpful and/or useful
drop me a thank you; better still, consider contributing something back to the
IDL community yourself.

A description of each IDL function is given in Tutorial.txt

Good luck,

Randall Skelton
skelton@atm.ox.ac.uk
