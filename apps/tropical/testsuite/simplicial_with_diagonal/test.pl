my $l1 = uniform_linear_space<Max>(2,1);
my $l2 = uniform_linear_space<Min>(3,2);

compare_object("1", simplicial_with_diagonal($l1))
	and
compare_object("2", simplicial_with_diagonal($l2));
