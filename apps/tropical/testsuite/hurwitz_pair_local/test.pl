my @h1 = hurwitz_pair_local<Max>(1, (new Vector<Int>(1,1,1,1,-4)), new RationalCurve(N_LEAVES=>5, INPUT_STRING=>"(1,2)"), Verbose=>0);
my @h2 = hurwitz_pair_local<Min>(1, (new Vector<Int>(1,2,3,4,-10)), new RationalCurve(N_LEAVES=>5, INPUT_STRING=>"(1,2)"), Verbose=>0);

compare_object("1subdiv", $h1[0])
	and 
compare_object("1cycle", $h1[1])
	and
compare_object("2subdiv", $h2[0])
	and
compare_object("2cycle", $h2[1]);
