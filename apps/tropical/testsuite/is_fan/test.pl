compare_output {
	for my $i (qw(empty ptcoll linear torus l32 c42 planarcurve)) {
		my $x = load("$i");
		print $x->name(), ": ",$x->is_fan(),"\n";
	}
} 'out'
