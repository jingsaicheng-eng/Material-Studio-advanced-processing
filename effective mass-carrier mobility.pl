#!perl
use strict;
use warnings;
use MaterialsScript qw(:all);

# -------------------------------------------------------------------------
# Effective mass from a CASTEP .bands file.
#
# CASTEP .bands contains:
#   - k-points in fractional reciprocal-lattice coordinates
#   - eigenvalues in atomic units (Hartree)
#
# This script:
#   1. Reads the .bands text file.
#   2. Converts eigenvalues to eV.
#   3. Converts fractional k-points to cumulative k distance in A^-1.
#   4. Finds VBM and CBM automatically, or uses a manually selected band edge.
#   5. Fits E(k) = a*(k-k0)^2 + b*(k-k0) + c near the band edge.
#   6. Calculates m*/m_e = 7.619964231 / (d2E/dk2), where d2E/dk2 = 2a.
#
# Notes:
#   - The result is the effective mass along the sampled band-path direction,
#     not the full effective-mass tensor.
#   - At a high-symmetry vertex, left and right fits represent different
#     reciprocal-space directions. Use the side matching the desired direction.
# -------------------------------------------------------------------------

# ================= User inputs =================

# Use a full Windows path. Forward slashes are safest in Materials Studio Perl.
my $BandsFile = "D:/msproject/bb_Files/Documents/GaAs CASTEP Energy/GaAs_BandStr.bands";

# CASTEP .bands cell vectors are usually written in atomic units.
# Use "bohr" or "A". If your effective mass looks scaled incorrectly, check this.
my $CellVectorUnit = "bohr";

# "AutoBandEdges" finds VBM and CBM from the Fermi energy.
# "Manual" analyzes one selected spin/band/k-point.
my $Mode = "AutoBandEdges";

# Manual mode settings. Band, spin, and k-point indices are 1-based.
my $ManualSpin = 1;
my $ManualBand = 1;
my $ManualKPoint = 1;
my $ManualCarrier = "electron";  # electron or hole

# Number of points used for one-sided fits, including the band-edge point.
# Needs at least 3. 4 or 5 is usually reasonable.
my $OneSideFitPoints = 5;

# Number of points on each side for the centered fit.
# 2 gives 5 points total if both sides are available.
my $CenteredFitPointsEachSide = 2;

# Energy tolerance for classifying occupied/unoccupied states around Fermi.
my $FermiTolerance_eV = 1.0e-6;

my $OutputStdName = "EffectiveMass.std";

# ================= Constants =================
my $HARTREE_TO_EV = 27.211386245988;
my $BOHR_TO_A = 0.529177210903;
my $HBAR2_OVER_ME_EV_A2 = 7.619964231073;  # hbar^2 / m_e
my $PI = 3.14159265358979323846;

# ================= Main =================
my $data = read_castep_bands($BandsFile);
my @cell_A = cell_vectors_in_angstrom($data->{cell}, $CellVectorUnit);
my @recip = reciprocal_vectors_2pi(\@cell_A);
my @k_cart = map { frac_to_cart($_, \@recip) } @{$data->{kfrac}};
my @kdist = cumulative_k_distance(\@k_cart);

my @jobs;
if ($Mode =~ /^auto/i) {
    my ($vbm, $cbm) = find_band_edges($data);
    die "Could not find VBM below Fermi energy.\n" unless defined $vbm;
    die "Could not find CBM above Fermi energy. The system may be metallic, "
        . "or the Fermi energy is not between valence and conduction bands.\n"
        unless defined $cbm;

    push @jobs, make_job("hole", $vbm);
    push @jobs, make_job("electron", $cbm);
} elsif ($Mode =~ /^manual/i) {
    push @jobs, {
        carrier => lc($ManualCarrier),
        spin => $ManualSpin - 1,
        band => $ManualBand - 1,
        kindex => $ManualKPoint - 1,
        energy_eV => $data->{eig}[$ManualSpin - 1][$ManualKPoint - 1][$ManualBand - 1]
    };
} else {
    die "Mode must be AutoBandEdges or Manual.\n";
}

my @results;
foreach my $job (@jobs) {
    foreach my $side ("center", "left", "right") {
        my $fit = fit_effective_mass($data, \@kdist, $job, $side);
        push @results, $fit if defined $fit;
    }
}

write_results_std($data, \@cell_A, \@recip, \@kdist, \@jobs, \@results);
print_results($data, \@jobs, \@results);

# ================= Reading .bands =================
sub read_castep_bands {
    my ($file) = @_;
    open my $fh, "<", $file or die "Could not open '$file': $!\n";

    my $nk = last_number(next_nonempty_line($fh));
    my $nspin = last_number(next_nonempty_line($fh));
    my $nelec = last_number(next_nonempty_line($fh));
    my $nband = last_number(next_nonempty_line($fh));
    my $fermi_ha = last_number(next_nonempty_line($fh));

    my $cell_header = next_nonempty_line($fh);
    die "Expected 'Unit cell vectors' line, got: $cell_header\n"
        unless $cell_header =~ /unit\s+cell\s+vectors/i;

    my @cell;
    for (my $i = 0; $i < 3; ++$i) {
        my $line = next_nonempty_line($fh);
        my @v = numbers_from_line($line);
        die "Could not read unit-cell vector from line: $line\n" unless @v >= 3;
        push @cell, [@v[0..2]];
    }

    my (@kfrac, @weight, @eig);
    for (my $ik_read = 0; $ik_read < $nk; ++$ik_read) {
        my $line = next_nonempty_line($fh);
        die "Expected K-point line, got: $line\n" unless $line =~ /k-?point/i;

        my @tok = split(/\s+/, trim($line));
        die "Malformed K-point line: $line\n" unless @tok >= 6;
        my $ik = int($tok[1]) - 1;
        $kfrac[$ik] = [$tok[2] + 0.0, $tok[3] + 0.0, $tok[4] + 0.0];
        $weight[$ik] = $tok[5] + 0.0;

        for (my $is = 0; $is < $nspin; ++$is) {
            my $spin_line = next_nonempty_line($fh);
            die "Expected Spin component line, got: $spin_line\n"
                unless $spin_line =~ /spin\s+component/i;

            for (my $ib = 0; $ib < $nband; ++$ib) {
                my $e_line = next_nonempty_line($fh);
                my @nums = numbers_from_line($e_line);
                die "Could not read eigenvalue from line: $e_line\n" unless @nums;
                $eig[$is][$ik][$ib] = $nums[-1] * $HARTREE_TO_EV;
            }
        }
    }

    close $fh;
    return {
        file => $file,
        nk => int($nk),
        nspin => int($nspin),
        nelec => $nelec + 0.0,
        nband => int($nband),
        fermi_eV => $fermi_ha * $HARTREE_TO_EV,
        cell => \@cell,
        kfrac => \@kfrac,
        weight => \@weight,
        eig => \@eig
    };
}

sub next_nonempty_line {
    my ($fh) = @_;
    while (defined(my $line = <$fh>)) {
        $line =~ s/\r?\n$//;
        next if $line =~ /^\s*$/;
        return $line;
    }
    die "Unexpected end of file.\n";
}

sub last_number {
    my ($line) = @_;
    my @nums = numbers_from_line($line);
    die "No number found in line: $line\n" unless @nums;
    return $nums[-1];
}

sub numbers_from_line {
    my ($line) = @_;
    my @nums = ($line =~ /[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][-+]?\d+)?/g);
    foreach my $v (@nums) {
        $v =~ tr/Dd/Ee/;
    }
    return @nums;
}

sub trim {
    my ($s) = @_;
    $s =~ s/^\s+//;
    $s =~ s/\s+$//;
    return $s;
}

# ================= Band edge selection =================
sub find_band_edges {
    my ($data) = @_;
    my ($vbm, $cbm);
    my $ef = $data->{fermi_eV};

    for (my $is = 0; $is < $data->{nspin}; ++$is) {
        for (my $ik = 0; $ik < $data->{nk}; ++$ik) {
            for (my $ib = 0; $ib < $data->{nband}; ++$ib) {
                my $e = $data->{eig}[$is][$ik][$ib];
                if ($e <= $ef + $FermiTolerance_eV) {
                    if (!defined($vbm) || $e > $vbm->{energy_eV}) {
                        $vbm = {
                            spin => $is, kindex => $ik, band => $ib, energy_eV => $e
                        };
                    }
                }
                if ($e >= $ef - $FermiTolerance_eV) {
                    if (!defined($cbm) || $e < $cbm->{energy_eV}) {
                        $cbm = {
                            spin => $is, kindex => $ik, band => $ib, energy_eV => $e
                        };
                    }
                }
            }
        }
    }
    return ($vbm, $cbm);
}

sub make_job {
    my ($carrier, $edge) = @_;
    return {
        carrier => $carrier,
        spin => $edge->{spin},
        band => $edge->{band},
        kindex => $edge->{kindex},
        energy_eV => $edge->{energy_eV}
    };
}

# ================= Effective mass fit =================
sub fit_effective_mass {
    my ($data, $kdist, $job, $side) = @_;
    my $ik0 = $job->{kindex};
    my @idx = fit_indices($data->{nk}, $ik0, $side);
    return undef unless @idx >= 3;

    my (@x, @y, @used);
    my %seen_x;
    foreach my $ik (@idx) {
        my $xx = $kdist->[$ik] - $kdist->[$ik0];
        my $key = sprintf("%.12g", $xx);
        next if $seen_x{$key}++;

        push @x, $xx;
        push @y, $data->{eig}[$job->{spin}][$ik][$job->{band}];
        push @used, $ik;
    }
    return undef unless @x >= 3;

    my ($a, $b, $c) = quadratic_fit(\@x, \@y);
    return undef unless defined $a;

    my $curvature = 2.0 * $a;  # eV A^2
    return undef if abs($curvature) < 1.0e-14;

    my $m_signed = $HBAR2_OVER_ME_EV_A2 / $curvature;
    my $m_carrier = abs($m_signed);
    my $r2 = fit_r2(\@x, \@y, $a, $b, $c);

    return {
        carrier => $job->{carrier},
        side => $side,
        spin => $job->{spin},
        band => $job->{band},
        kindex => $job->{kindex},
        energy_eV => $job->{energy_eV},
        npoints => scalar(@x),
        used_kpoints => join(" ", map { $_ + 1 } @used),
        kmin_Ainv => min(@x),
        kmax_Ainv => max(@x),
        a => $a,
        b => $b,
        c => $c,
        curvature_eV_A2 => $curvature,
        m_signed_me => $m_signed,
        m_carrier_me => $m_carrier,
        r2 => $r2
    };
}

sub fit_indices {
    my ($nk, $ik0, $side) = @_;
    my @idx;

    if ($side eq "left") {
        my $start = $ik0 - $OneSideFitPoints + 1;
        $start = 0 if $start < 0;
        @idx = ($start .. $ik0);
    } elsif ($side eq "right") {
        my $end = $ik0 + $OneSideFitPoints - 1;
        $end = $nk - 1 if $end >= $nk;
        @idx = ($ik0 .. $end);
    } elsif ($side eq "center") {
        my $start = $ik0 - $CenteredFitPointsEachSide;
        my $end = $ik0 + $CenteredFitPointsEachSide;
        $start = 0 if $start < 0;
        $end = $nk - 1 if $end >= $nk;
        @idx = ($start .. $end);
    } else {
        die "Unknown fit side '$side'.\n";
    }
    return @idx;
}

sub quadratic_fit {
    my ($x, $y) = @_;
    my $n = scalar(@$x);
    return undef unless $n >= 3;

    my ($s0, $s1, $s2, $s3, $s4) = ($n, 0, 0, 0, 0);
    my ($t0, $t1, $t2) = (0, 0, 0);
    for (my $i = 0; $i < $n; ++$i) {
        my $xi = $x->[$i];
        my $yi = $y->[$i];
        my $x2 = $xi * $xi;
        $s1 += $xi;
        $s2 += $x2;
        $s3 += $x2 * $xi;
        $s4 += $x2 * $x2;
        $t0 += $yi;
        $t1 += $xi * $yi;
        $t2 += $x2 * $yi;
    }

    my @A = (
        [$s4, $s3, $s2, $t2],
        [$s3, $s2, $s1, $t1],
        [$s2, $s1, $s0, $t0]
    );
    my @sol = solve_3x3(\@A);
    return @sol;
}

sub solve_3x3 {
    my ($A) = @_;
    for (my $col = 0; $col < 3; ++$col) {
        my $pivot = $col;
        for (my $r = $col + 1; $r < 3; ++$r) {
            $pivot = $r if abs($A->[$r][$col]) > abs($A->[$pivot][$col]);
        }
        return () if abs($A->[$pivot][$col]) < 1.0e-20;
        if ($pivot != $col) {
            my $tmp = $A->[$col];
            $A->[$col] = $A->[$pivot];
            $A->[$pivot] = $tmp;
        }
        my $p = $A->[$col][$col];
        for (my $j = $col; $j < 4; ++$j) {
            $A->[$col][$j] /= $p;
        }
        for (my $r = 0; $r < 3; ++$r) {
            next if $r == $col;
            my $f = $A->[$r][$col];
            for (my $j = $col; $j < 4; ++$j) {
                $A->[$r][$j] -= $f * $A->[$col][$j];
            }
        }
    }
    return ($A->[0][3], $A->[1][3], $A->[2][3]);
}

sub fit_r2 {
    my ($x, $y, $a, $b, $c) = @_;
    my $n = scalar(@$y);
    return "" if $n < 2;

    my $mean = 0.0;
    foreach my $yy (@$y) { $mean += $yy; }
    $mean /= $n;

    my ($ss_res, $ss_tot) = (0.0, 0.0);
    for (my $i = 0; $i < $n; ++$i) {
        my $pred = $a * $x->[$i] * $x->[$i] + $b * $x->[$i] + $c;
        $ss_res += ($y->[$i] - $pred) * ($y->[$i] - $pred);
        $ss_tot += ($y->[$i] - $mean) * ($y->[$i] - $mean);
    }
    return $ss_tot > 0.0 ? 1.0 - $ss_res / $ss_tot : "";
}

# ================= k-space geometry =================
sub cell_vectors_in_angstrom {
    my ($cell, $unit) = @_;
    my $scale = 1.0;
    $scale = $BOHR_TO_A if $unit =~ /^bohr$/i || $unit =~ /^au$/i;
    my @out;
    foreach my $v (@$cell) {
        push @out, [$v->[0] * $scale, $v->[1] * $scale, $v->[2] * $scale];
    }
    return @out;
}

sub reciprocal_vectors_2pi {
    my ($cell) = @_;
    my ($a, $b, $c) = @$cell;
    my $vol = dot($a, cross($b, $c));
    die "Unit cell volume is zero.\n" if abs($vol) < 1.0e-20;

    my $factor = 2.0 * $PI / $vol;
    my $b1 = scale_vec(cross($b, $c), $factor);
    my $b2 = scale_vec(cross($c, $a), $factor);
    my $b3 = scale_vec(cross($a, $b), $factor);
    return ($b1, $b2, $b3);
}

sub frac_to_cart {
    my ($kfrac, $recip) = @_;
    my ($b1, $b2, $b3) = @$recip;
    return [
        $kfrac->[0] * $b1->[0] + $kfrac->[1] * $b2->[0] + $kfrac->[2] * $b3->[0],
        $kfrac->[0] * $b1->[1] + $kfrac->[1] * $b2->[1] + $kfrac->[2] * $b3->[1],
        $kfrac->[0] * $b1->[2] + $kfrac->[1] * $b2->[2] + $kfrac->[2] * $b3->[2]
    ];
}

sub cumulative_k_distance {
    my ($kcart) = @_;
    my @dist = (0.0);
    for (my $i = 1; $i < @$kcart; ++$i) {
        $dist[$i] = $dist[$i - 1] + norm(sub_vec($kcart->[$i], $kcart->[$i - 1]));
    }
    return @dist;
}

sub cross {
    my ($u, $v) = @_;
    return [
        $u->[1] * $v->[2] - $u->[2] * $v->[1],
        $u->[2] * $v->[0] - $u->[0] * $v->[2],
        $u->[0] * $v->[1] - $u->[1] * $v->[0]
    ];
}

sub dot {
    my ($u, $v) = @_;
    return $u->[0] * $v->[0] + $u->[1] * $v->[1] + $u->[2] * $v->[2];
}

sub scale_vec {
    my ($v, $s) = @_;
    return [$v->[0] * $s, $v->[1] * $s, $v->[2] * $s];
}

sub sub_vec {
    my ($u, $v) = @_;
    return [$u->[0] - $v->[0], $u->[1] - $v->[1], $u->[2] - $v->[2]];
}

sub norm {
    my ($v) = @_;
    return sqrt(dot($v, $v));
}

# ================= Output =================
sub write_results_std {
    my ($data, $cell_A, $recip, $kdist, $jobs, $results) = @_;
    my $std = Documents->New($OutputStdName);

    $std->Cell(0, 0) = "CASTEP effective mass from .bands";
    $std->Cell(1, 0) = "Bands file";
    $std->Cell(1, 1) = $data->{file};
    $std->Cell(2, 0) = "Mode";
    $std->Cell(2, 1) = $Mode;
    $std->Cell(3, 0) = "Cell vector unit input";
    $std->Cell(3, 1) = $CellVectorUnit;
    $std->Cell(4, 0) = "k-points";
    $std->Cell(4, 1) = $data->{nk};
    $std->Cell(5, 0) = "spin components";
    $std->Cell(5, 1) = $data->{nspin};
    $std->Cell(6, 0) = "bands";
    $std->Cell(6, 1) = $data->{nband};
    $std->Cell(7, 0) = "electrons";
    $std->Cell(7, 1) = $data->{nelec};
    $std->Cell(8, 0) = "Fermi energy eV";
    $std->Cell(8, 1) = $data->{fermi_eV};
    $std->Cell(9, 0) = "mass formula";
    $std->Cell(9, 1) = "m*/me = 7.619964231 / curvature(eV A^2)";

    my $row = 11;
    $std->Cell($row, 0) = "Selected band edges";
    $row++;
    $std->Cell($row, 0) = "carrier";
    $std->Cell($row, 1) = "spin";
    $std->Cell($row, 2) = "band";
    $std->Cell($row, 3) = "k index";
    $std->Cell($row, 4) = "kfrac x";
    $std->Cell($row, 5) = "kfrac y";
    $std->Cell($row, 6) = "kfrac z";
    $std->Cell($row, 7) = "k distance A^-1";
    $std->Cell($row, 8) = "energy eV";
    $std->Cell($row, 9) = "E-Ef eV";
    $row++;
    foreach my $job (@$jobs) {
        my $ik = $job->{kindex};
        my $kf = $data->{kfrac}[$ik];
        $std->Cell($row, 0) = $job->{carrier};
        $std->Cell($row, 1) = $job->{spin} + 1;
        $std->Cell($row, 2) = $job->{band} + 1;
        $std->Cell($row, 3) = $ik + 1;
        $std->Cell($row, 4) = $kf->[0];
        $std->Cell($row, 5) = $kf->[1];
        $std->Cell($row, 6) = $kf->[2];
        $std->Cell($row, 7) = $kdist->[$ik];
        $std->Cell($row, 8) = $job->{energy_eV};
        $std->Cell($row, 9) = $job->{energy_eV} - $data->{fermi_eV};
        $row++;
    }

    $row += 2;
    $std->Cell($row, 0) = "Fit results";
    $row++;
    my @headers = (
        "carrier", "fit side", "spin", "band", "k index", "points",
        "used k-points", "k min A^-1", "k max A^-1", "edge energy eV",
        "curvature eV A^2", "m signed / me", "carrier mass / me", "R2",
        "a", "b", "c"
    );
    for (my $i = 0; $i < @headers; ++$i) {
        $std->Cell($row, $i) = $headers[$i];
    }
    $row++;
    foreach my $r (@$results) {
        $std->Cell($row, 0) = $r->{carrier};
        $std->Cell($row, 1) = $r->{side};
        $std->Cell($row, 2) = $r->{spin} + 1;
        $std->Cell($row, 3) = $r->{band} + 1;
        $std->Cell($row, 4) = $r->{kindex} + 1;
        $std->Cell($row, 5) = $r->{npoints};
        $std->Cell($row, 6) = $r->{used_kpoints};
        $std->Cell($row, 7) = $r->{kmin_Ainv};
        $std->Cell($row, 8) = $r->{kmax_Ainv};
        $std->Cell($row, 9) = $r->{energy_eV};
        $std->Cell($row, 10) = $r->{curvature_eV_A2};
        $std->Cell($row, 11) = $r->{m_signed_me};
        $std->Cell($row, 12) = $r->{m_carrier_me};
        $std->Cell($row, 13) = $r->{r2};
        $std->Cell($row, 14) = $r->{a};
        $std->Cell($row, 15) = $r->{b};
        $std->Cell($row, 16) = $r->{c};
        $row++;
    }
}

sub print_results {
    my ($data, $jobs, $results) = @_;
    print "\nEffective mass analysis\n";
    print "Bands file: $data->{file}\n";
    print "Fermi energy: " . fmt($data->{fermi_eV}) . " eV\n";
    print "Output STD: $OutputStdName\n";
    foreach my $job (@$jobs) {
        print "\n" . uc($job->{carrier}) . " edge: spin "
            . ($job->{spin} + 1) . ", band " . ($job->{band} + 1)
            . ", k-point " . ($job->{kindex} + 1)
            . ", E = " . fmt($job->{energy_eV}) . " eV\n";
        foreach my $r (@$results) {
            next unless $r->{carrier} eq $job->{carrier}
                && $r->{spin} == $job->{spin}
                && $r->{band} == $job->{band}
                && $r->{kindex} == $job->{kindex};
            print "  $r->{side}: m*/me = "
                . fmt($r->{m_carrier_me})
                . "  curvature = " . fmt($r->{curvature_eV_A2})
                . " eV A^2, R2 = " . fmt($r->{r2}) . "\n";
        }
    }
}

sub min {
    my @v = @_;
    my $m = $v[0];
    foreach my $x (@v) { $m = $x if $x < $m; }
    return $m;
}

sub max {
    my @v = @_;
    my $m = $v[0];
    foreach my $x (@v) { $m = $x if $x > $m; }
    return $m;
}

sub fmt {
    my ($v) = @_;
    return "" unless defined $v;
    return $v unless $v =~ /^[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?$/;
    return sprintf("%.10g", $v);
}
