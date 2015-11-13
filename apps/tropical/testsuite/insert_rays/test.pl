my $c1 = uniform_linear_space<Max>(3,2);
my $c2 = uniform_linear_space<Min>(3,2);
my $c3 = cross_variety<Max>(3,2);
my $c4 = affine_linear_space<Min>([[0,1,0,1],[0,0,1,1]]);

my $r1 = new Matrix<Rational>([[0,0,-1,-1,0],[0,0,-1,1,1]]);
my $r2 = new Matrix<Rational>([[1,0,1,1,0]]);
my $r3 = new Matrix<Rational>([[0,0,1,-1,0],[1,0,2,-2,1]]);
my $r4 = new Matrix<Rational>([[1,0,1,1,2]]);

compare_object("1", insert_rays($c1,$r1));

compare_object("2", insert_rays($c2,$r2));

compare_object("3", insert_rays($c3,$r3));

compare_object("4", insert_rays($c4,$r4));
