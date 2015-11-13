my $c1 = new RationalCurve(N_LEAVES=>5,INPUT_STRING=>"(1,2) + (4,5)");
my $c2 = new RationalCurve(N_LEAVES=>4, INPUT_STRING=>"(1,2)");

compare_object( "1", rational_curve_immersion<Max>([[1,0],[1,0],[-4,0],[1,0],[1,0]], $c1));

compare_object( "2", rational_curve_immersion<Min>([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]],$c2));
