# This file goes through all .cc- files in apps/atint/src 
# It checks the dbg-status of each file (i.e. which line of donotlog, dolog and dotrace is not commented out)
# and comments out debug statements accordingly

use File::Slurp;

my $dirname = "./apps/atint/src";
opendir(DIR, $dirname);

my @files = readdir(DIR);

for my $f (@files) {
  
  #We are only interested in files ending in .cc
  if($f =~ /.*\.cc$/) {
    #Open file
    my $file = read_file($dirname."/".$f);

  
    #Check the debug status of the file
    my $dotrace = 1;
    my $dolog = 1;
    if($file =~ /\/\/.*using namespace atintlog\:\:dotrace.*/) { 
      $dotrace = 0;
    }
    if($file =~ /\/\/.*using namespace atintlog\:\:dolog.*/) { 
      $dolog = $dotrace;
    }

    print "$f: ".($dotrace? "Tracing" : ($dolog? "Logging" : "No debugging"))."\n";

    #Comment out trace statements if necessary
    if(!$dotrace)  {
      $file =~ s/^(\s*)dbgtrace(.*)/$1\/\/dbgtrace$2/mg;
    }
    #Otherwise uncomment trace statements
    else {
      $file =~ s/^(\s*)\/\/(\s*dbgtrace.*)/$1$2/mg;
    }
    #Same for log statements
    if(!$dolog) {
      $file =~ s/^(\s*)dbglog(.*)/$1\/\/dbgtrace$2/mg;
    }
    else {
      $file =~ s/^(\s*)\/\/(\s*dbglog.*)/$1$2/mg;
    }

    #Write changes to file
    write_file($dirname."/".$f,$file);
    
  }



}

closedir(DIR);