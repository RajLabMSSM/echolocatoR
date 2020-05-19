================================================================
========   bigWigToBedGraph   ====================================
================================================================
### kent source version 393 ###
bigWigToBedGraph - Convert from bigWig to bedGraph format.
usage:
   bigWigToBedGraph in.bigWig out.bedGraph
options:
   -chrom=chr1 - if set restrict output to given chromosome
   -start=N - if set, restrict output to only that over start
   -end=N - if set, restict output to only that under end
   -udcDir=/dir/to/cache - place to put cache for remote bigBed/bigWigs





This file is from:

http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/README.txt

This directory contains applications for stand-alone use, 
built on a Mac OSX 10.14.6 intel machine.  (Mojave)
Darwin Kernel Version 18.7.0, gcc version:
Apple clang version 11.0.0 (clang-1100.0.33.8)
Target: x86_64-apple-darwin18.7.0
Thread model: posix

kent source tree v392 January 2020.

For help on the bigBed and bigWig applications see:
http://genome.ucsc.edu/goldenPath/help/bigBed.html
http://genome.ucsc.edu/goldenPath/help/bigWig.html

View the file 'FOOTER.txt' to see the usage statement for 
each of the applications.

The shared libraries used by these binaries are:  (from: otool -L <binary>)

/usr/lib/libc++.1.dylib (compatibility version 1.0.0, current version 400.9.4)
/usr/lib/libSystem.B.dylib (compatibility version 1.0.0, current version 1252.250.1)

##############################################################################
Thank you to Bob Harris for permission to distribute a binary
version of the lastz and lastz_D programs, from:

   https://github.com/lastz/lastz

Version 1.04.00 as of April 2018:

-rwxrwxr-x 1 514360 Apr  6 11:24 lastz-1.04.00
-rwxrwxr-x 1 514360 Apr  6 11:24 lastz_D-1.04.00

$ gmd5sum lastz*
4aa388bf0743e48d45704dc9194d5888  lastz-1.04.00
bae26692e8a35313a1f149e37c8c7e6e  lastz_D-1.04.00

##############################################################################
This entire directory can by copied with the rsync command
into the local directory ./

rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/macOSX.x86_64/ ./

Or from our mirror site:

rsync -aP rsync://hgdownload-sd.soe.ucsc.edu/genome/admin/exe/macOSX.x86_64/ ./

Individual programs can by copied by adding their name, for example:

rsync -aP \
   rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/macOSX.x86_64/faSize ./
