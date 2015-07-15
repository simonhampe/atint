my $c1 = cross_variety<Max>(3,2);
my $c2 = uniform_linear_space<Min>(4,3);

compare_object("1", triangulate_cycle($c1))
	and
compare_object("2", triangulate_cycle($c2));
