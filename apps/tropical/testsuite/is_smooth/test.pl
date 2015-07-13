my $a = uniform_linear_space<Max>(3,2);
my $b = uniform_linear_space<Min>(2,1);
my $c = cartesian_product( (uniform_linear_space<Max>(2,1)), (projective_torus<Max>(1)));
my $d = affine_linear_space<Min>( [[0,2,3,1,0],[0,1,-2,0,1]]);
my $e = point_collection<Max>([[0,0,0,0]], [1]);
my $f = cross_variety<Min>(3,2,0);

my @s1 = is_smooth($a);
my @s2 = is_smooth($b);
my @s3 = is_smooth($c);
my @s4 = is_smooth($d);
my @s5 = is_smooth($e);
my @s6 = is_smooth($f);

compare_values("1b",1,$s1[0])
	and
compare_values("2b",1,$s2[0])
	and
compare_values("3b",1,$s3[0])
	and
compare_values("4b",1,$s4[0])
	and
compare_values("5b",1,$s5[0])
	and
compare_values("6b",0,$s6[0])
	and
compare_object("1mat",$s1[1])
	and
compare_object("2mat",$s2[1])
	and
compare_object("3mat",$s3[1])
	and
compare_object("4mat",$s4[1])
	and
compare_object("5mat",$s5[1])
	and
compare_object("1morph",$s1[2])
	and
compare_object("2morph",$s2[2])
	and
compare_object("3morph",$s3[2])
	and
compare_object("4morph",$s4[2])
	and
compare_object("5morph",$s5[2]);

