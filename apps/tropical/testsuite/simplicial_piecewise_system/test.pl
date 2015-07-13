my $l1 = uniform_linear_space<Max>(2,1);
my $l2 = uniform_linear_space<Min>(3,2);

compare_data("1", simplicial_piecewise_system($l1))
	and
compare_data("2", simplicial_piecewise_system($l2));
