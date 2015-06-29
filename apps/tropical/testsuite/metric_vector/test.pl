my $r1 = new RationalCurve(N_LEAVES=>3,INPUT_STRING=>"");
my $r2 = new RationalCurve(INPUT_STRING=>"(1,2) + 2*(4,5)");
my $r3 = new RationalCurve(INPUT_STRING=>"7*(1,2,3)",N_LEAVES=>5);

compare_data('1', $r1->metric_vector()) 
		and
compare_data('2', $r2->metric_vector())
		and
compare_data('3', $r3->metric_vector());
