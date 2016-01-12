my $c1 = uniform_linear_space<Max>(3,1);
my $c2 = cross_variety<Min>(3,2);
my $c3 = new Cycle<Max>(VERTICES=>[[1,0,0,0],[1,0,1,0],[0,0,0,1],[0,0,-1,-1],[0,0,1,-1]],MAXIMAL_POLYTOPES=>[[0,1],[0,2],[0,3],[1,2],[1,4]], WEIGHTS=>[1,1,1,1,1]); 
my $c4 = empty_cycle<Min>(5);

compare_object("1", recession_fan($c1));

compare_object("2", recession_fan($c2));

compare_object("3", recession_fan($c3));

compare_object("4", recession_fan($c4));
