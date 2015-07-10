my $p1 = uniform_linear_space<Max>(3,2);

compare_data("1", cutting_functions($p1, new Vector<Integer>(1,1,1,1)))
	and
compare_data("2", cutting_functions(load("p2"), new Vector<Integer>(1,0,1,1,0,1)));
