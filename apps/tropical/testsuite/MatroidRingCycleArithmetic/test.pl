my $m1 = new matroid::Matroid(N_ELEMENTS=>4,BASES=>[[0,1],[0,2],[1,3],[2,3]]);
my $r1 = matroid_ring_cycle<Max>($m1);
my $r2 = matroid_ring_cycle<Max>( matroid::uniform_matroid(3,4));
my $result1 = 3*$r1 - 2*($r2 * $r2);
my $result2 = $result1 * $r2;

my $r3 = matroid_ring_cycle<Min>(matroid::matroid_from_graph(complete(5)));
my $r4 = matroid_ring_cycle<Min>( matroid::direct_sum(matroid::uniform_matroid(4,5), matroid::uniform_matroid(4,5)));
my $result3 = - (3*$r3) * $r4;

my $r5 = zero_in_matroid_ring<Max>(4);
my $result4 = $result1 + 17*$r5;

compare_object("1", $result1);
compare_object("2", $result2);
compare_object("3", $result3);
compare_object("4", $result4);
