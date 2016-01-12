my $t = cross_variety<Max>(3,2);
my $w1 = new Cycle<Min>(VERTICES=>[[1,0,0]],MAXIMAL_POLYTOPES=>[[0]],LINEALITY_SPACE=>[[0,1,0]],WEIGHTS=>[1]);
my $w2 = new Cycle<Min>(VERTICES=>[[1,2,3]],MAXIMAL_POLYTOPES=>[[0]],LINEALITY_SPACE=>[[0,1,0]],WEIGHTS=>[2]);

check_boolean("1", check_cycle_equality($t,$t));

check_boolean("2", check_cycle_equality($w1,$w2,0));
