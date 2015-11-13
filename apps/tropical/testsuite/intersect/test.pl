my $a = uniform_linear_space<Max>(2,1);
my $b = new Cycle<Max>(VERTICES=>[[1,0,0,0],[1,0,1,0],[0,0,0,1],[0,0,-1,-1],[0,0,1,-1]],MAXIMAL_POLYTOPES=>[[0,1],[0,2],[0,3],[1,2],[1,4]],WEIGHTS=>[1,1,1,1,1]);
my $c = uniform_linear_space<Min>(2,1);
my $d = empty_cycle<Max>(2);
my $e = uniform_linear_space<Max>(3,2);
my $f = new Hypersurface<Max>(POLYNOMIAL=>toTropicalPolynomial("max(-49/1000+3*x0 , -523/250+2*x0+x1 , -501/250+2*x0+x2 , -2+2*x0+x3 , -801/100+x0+2*x1 , -499/100+x0+x1+x2 , -5997/1000+x0+x1+x3 , -7919/1000+x0+2*x2 , -5981/1000+x0+x2+x3 , -799/100+x0+2*x3 , -18023/1000+3*x1 , -6021/500+2*x1+x2 , -1759/125+2*x1+x3 , -12091/1000+x1+2*x2 , -2761/250+x1+x2+x3 , -6983/500+x1+2*x3 , -17933/1000+3*x2 , -6979/500+2*x2+x3 , -7023/500+x2+2*x3 , -2249/125+3*x3)"));
my $g = affine_linear_space<Min>([[0,1,1,0,0],[0,1,0,1,0]]);
my $h = affine_linear_space<Min>([[0,1,-1,0,0],[0,0,0,0,1]]);

compare_object("1", intersect($a,$b));

compare_object("2", intersect($c,$c));

compare_object("3", intersect($a,$d));

compare_object("4", intersect($e,$f));

compare_object("5", intersect($g,$h));
