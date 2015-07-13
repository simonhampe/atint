my $r1 = new RationalCurve(N_LEAVES=>4, INPUT_STRING=>"(1,2)");
my $s1 = insert_leaves($r1, new Vector<Int>([0,1,0,0]));

my $r2 = new RationalCurve(N_LEAVES=>6, INPUT_STRING=>"(1,2) + (3,4) + (5,6)");
my $s2 = insert_leaves($r2,[0,2,3,1,3,2,0]);

compare_object("1", $s1)
	and
compare_object("2", $s2);
