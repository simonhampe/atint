for my $ID (qw(empty linear torus l32)) {
   check_boolean($ID, load($ID)->is_fan());
}

for my $ID (qw(ptcoll c42 planarcurve)) {
   check_boolean($ID, !load($ID)->is_fan());
}
