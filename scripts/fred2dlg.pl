#!/usr/bin/env perl

use Fortran::F90Namelist;
my $nl = Fortran::F90Namelist->new() or die "Couldn't get object\n";

$nl->parse(file     => 'aorsa2d.in.fred',
           all      => 1,
           namelist => 'nlist');

print $nl->output(format => 'idl', name => 'fred');




