my $a = uniform_linear_space<Max>(2,1);
my $b = shift_cycle($a, [0,1,0]);
my $c = uniform_linear_space<Min>(3,2);
my $d = cross_variety<Min>(3,2);
my $e = shift_cycle($a, [0,1,-1]);

my @i1 = intersect_check_transversality($a,$b);
my @i2 = intersect_check_transversality($c,$d);
my @i3 = intersect_check_transversality($a,$e);

compare_values("1b",0,$i1[1]);

compare_values("2b",0,$i2[1]);

compare_values("3b",1,$i3[1]);

compare_object("1", $i1[0]);

compare_object("2", $i2[0]);

compare_object("3",$i3[0]);
