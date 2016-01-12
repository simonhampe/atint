my $r1 = new RationalCurve(N_LEAVES=>3,INPUT_STRING=>"");
my $r2 = new RationalCurve(N_LEAVES=>4,INPUT_STRING=>"(1,2)");

my $s1a = new RationalCurve(N_LEAVES=>5, INPUT_STRING=>"(1,2,3)");
my $s1b = new RationalCurve(N_LEAVES=>5, INPUT_STRING=>"(1,2)");
my $s1c = new RationalCurve(N_LEAVES=>5, INPUT_STRING=>"(1,2) + (4,5)");

compare_object('1', 3*$r1);

compare_object('2', 3*$r2);

compare_object('3', 2*$s1a + 4*$s1b);

compare_object('4', $s1c - $s1a);
