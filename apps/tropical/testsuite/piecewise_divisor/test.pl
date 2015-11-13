my $F1 = new Cycle<Max>(VERTICES=>[[1,0,0,0],[0,0,1,1],[0,0,1,-1],[0,0,-1,0]],MAXIMAL_POLYTOPES=>[[0,1],[0,2],[0,3]],WEIGHTS=>[1,1,2]);
my $F2 = uniform_linear_space<Min>(3,2);

compare_object("1", piecewise_divisor($F1,[[0,3]],[3]));

compare_object("2", piecewise_divisor($F2, [[0,1],[0,2]],[1,1]));
