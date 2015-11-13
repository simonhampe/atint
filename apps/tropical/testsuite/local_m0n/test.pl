my $c1a = new RationalCurve(N_LEAVES=>7, INPUT_STRING=>"(1,2) + (3,4)");
my $c1b = new RationalCurve(N_LEAVES=>7, INPUT_STRING=>"(1,3) + (1,3,7) + (2,5,6)");
my $c2 = new RationalCurve(N_LEAVES=>11, INPUT_STRING=>"(1,2) + (1,2,3,4) + (1,2,3,4,5) + (1,2,3,4,5,6) + (10,11) + (9,10,11) + (8,9,10,11)");

compare_object("1", local_m0n<Max>($c1a, $c1b));

compare_object("2", local_m0n<Min>( $c2));
