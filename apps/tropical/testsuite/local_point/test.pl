check_if_configured("bundled:cdd") or return;
prefer_now 'polytope::cdd';

my $a = uniform_linear_space<Max>(2,1);
my $b = cross_variety<Min>(3,2);

compare_object("1", local_point($a,[1,0,0,-1]));

compare_object("2", local_point($b,[1,0,2,2,2]));
