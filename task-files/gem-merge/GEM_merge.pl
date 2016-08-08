#####################################################################
#                                                                   #
#  MatrixStitcher                                                   #
#                                                                   #
#  Alex Feltus <ffeltus@clemson.edu>                                #
#  First Write: 04-10-2016                                          #
#  Last Revisit: 04-11-2016                                         #
#####################################################################

use strict;
use warnings;

### Usage
if (@ARGV < 1) {
   print "usage: perl MatrixSticher.pl <FULL_PATHNAME_TO_DATASETS>\n";
   die;
}

### Read all files in directory with a bunch of two-column files (WITH HEADER) with an ID in the first column and expression value in the next column seperated by whitespace
my $dirpath = $ARGV[0];
my %MatrixValues = ();
my @files = ();
my @OneFileData =();
my $ValueCounter = 0;

opendir (DIR, $dirpath) or die $!;
while (my $file = readdir(DIR)) {

   # Use a regular expression to ignore files beginning with a period or ends in '.pl'
   next if ($file =~ m/^\./);
   next if ($file =~ m/\.pl/);
   
   # Keep track of file names for column headers.  One could parse the file names here, I reckon.
   push (@files, $file);

   # Build the hash of arrays
   my $fullpath =  $dirpath."/".$file;
   @OneFileData = get_file_data ($fullpath);

   shift @OneFileData; #Remove the header
   foreach my $line (@OneFileData) {
        chomp $line;
        my @elements = split (/\s+/, $line);
        my $FirstColName = $elements[0];
        my $SecondColValue = $elements[1];
        push(@{$MatrixValues{$FirstColName}}, $SecondColValue);
   }
   @OneFileData = ();
   
   #Find the length of the arrays and see if any new keys were added or a key didn't exist.  Add 'NA' to all previous files if a label was added and an NA to the current file if a label was missed.

   $ValueCounter++;
   foreach my $RowID ( sort keys %MatrixValues ) {
        my $ValueCount = scalar @{$MatrixValues{$RowID}};
        #print "$ValueCount\t$ValueCounter\n";

   # Label wasn't in dataset but seen previously so add 'NA
   if ($ValueCount == $ValueCounter - 1 ) {
      push(@{$MatrixValues{$RowID}}, "NA");
      $ValueCount++;
      }

   # Label was seen for the first time so add 'NA' to all previous datasets and value to the new one
   if ($ValueCount < $ValueCounter  && $ValueCount == 1 ) {
                    my $StoreValue = shift @{$MatrixValues{$RowID}};
                    for (my $i=1; $i < $ValueCounter; $i++) {
                        push (@{$MatrixValues{$RowID}}, "NA");
                        }
                    push(@{$MatrixValues{$RowID}}, $StoreValue);
                    }
   }
}

closedir(DIR);


### Print the stitched matrix
# Print column headers
print "RowID\t";
foreach my $columnheader (@files) {
        print "$columnheader\t";
}
print "\n";

# Print stitched matrix values
foreach my $RowID ( sort keys %MatrixValues ) {
        print "$RowID\t";
        foreach my $ColumnValue (@{$MatrixValues{$RowID}}) {
                print "$ColumnValue\t";
        }
        print "\n";
}

exit;
###########################################################################
#                             SUBROUTINES                                 #
###########################################################################

# subroutine that opens a file
sub get_file_data {

    my ($filename) = @_;
    my @filedata = ();

    unless (open (GET_FILE_DATA, $filename)){
           print STDERR "Cannot open file \"$filename\"\n\n";
           exit;
    }

    @filedata = <GET_FILE_DATA>;

    close GET_FILE_DATA;

    return @filedata;

}


