check_if_configured("bundled:cdd") or return;
prefer_now 'polytope::cdd';

my $a = cross_variety<Max>(2,1,0);
my $b = new Cycle<Min>(VERTICES=>thomog([[1,0,0],[0,1,0],[0,1,1],[0,0,1],[0,-1,0],[0,-1,-1],[0,0,-1]]),MAXIMAL_POLYTOPES=>[[0,1],[0,2],[0,3],[0,4],[0,5],[0,6]],WEIGHTS=>[1,1,1,1,1,1]); 
my $c = new Cycle<Min>(VERTICES=>thomog([[1,0,0],[0,1,0],[0,1,1],[0,0,1],[0,-1,0],[0,-1,-1],[0,0,-1]]),MAXIMAL_POLYTOPES=>[[0,1],[0,2],[0,3],[0,4],[0,5],[0,6]],WEIGHTS=>[1,3,2,1,3,2]);

compare_object("1", decomposition_polytope($a));

compare_object("2", decomposition_polytope($b));

compare_object("3", decomposition_polytope($c));
