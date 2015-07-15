
my $r1 = new RationalCurve(N_LEAVES=>5, INPUT_STRING=>"(1,2)");
my $r2 = new RationalCurve(N_LEAVES=>6, INPUT_STRING=>"(1,2) + (3,4) + (5,6)");
my $r3 = new RationalCurve(N_LEAVES=>4, INPUT_STRING=>"");

my $v1 = new Vector<Rational>(0,-1,0,0,0,0,0);
my $v2 = new Vector<Rational>(0,2,1,1,0,1,1,0,2,0,0);
my $v3 = new Vector<Rational>(0,0,0,0);

compare_object("1", rational_curve_from_metric(new Vector<Rational>(0,1,1,1,1,1,1,0,0,0)))
	and
compare_object("2", rational_curve_from_metric(new Vector<Rational>(0,1,2,2,1,2,2,1,1,0)))
	and
compare_object("3", rational_curve_from_metric(new Vector<Rational>(0,0,0,0,0,0)))
	and
compare_object("1a", rational_curve_from_matroid_coordinates<Max>($v1))
	and
compare_object("2a", rational_curve_from_matroid_coordinates<Min>($v2))
	and
compare_object("3a", rational_curve_from_matroid_coordinates<Max>($v3))
	and
compare_data("1b", matroid_coordinates_from_curve<Max>($r1))
	and
compare_data("2b", matroid_coordinates_from_curve<Min>($r2))
	and
compare_data("3b", matroid_coordinates_from_curve<Max>($r3));
