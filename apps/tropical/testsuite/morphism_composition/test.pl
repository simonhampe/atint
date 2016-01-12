for my $i (1 .. 8) {
   my $f = load($i."f");
   my $g = load($i."g");
   compare_object($i, morphism_composition($f,$g));
}
