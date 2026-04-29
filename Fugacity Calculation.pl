#!perl

use strict;
use warnings;
use MaterialsScript qw(:all);

# -------------------------------------------------------------------------
# Materials Studio standalone fugacity converter.
# Author Jingsai Cheng
# E-mail jingsaicheng@gmail.com
# Scope:
#   Only converts pressure/composition/temperature to fugacity.
#   It does not run Sorption, FixedPressure, or any other MS calculation task.
#
# MS Sorption note:
#   FixedPressure needs fugacity in kPa:
#       $task->Fugacity($sorbate) = <Fugacity_kPa>;
#
# Method:
#   Peng-Robinson EOS with classical mixture rules.
#   For a mixture:
#       f_i = y_i * phi_i * P_total
#   where P is in kPa and f_i is returned in kPa.
# -------------------------------------------------------------------------

# ================= User input =================

my $Temperature_K = 300.0;

# Pressure input. Unit can be "MPa" or "kPa".
my $PressureInputUnit = "MPa";
my @TotalPressure_List = (100.0);

# Mode 1: calculate one specified gas or mixture.
my $CalculateAllSupportedPureGases = "No";
my @Components = (
    { Gas => "N2", Ratio => 1.0 },
);

# Mode 2: if CalculateAllSupportedPureGases = Yes, calculate each gas below
# as a pure gas at every pressure in TotalPressure_List.
my @PureGasList = qw(
    N2 O2 CO2 CO CH4 C2H2 C2H4 C2H6 C3H8
    H2 H2O NO NO2 N2O SO2 H2S NH3 Ar He Ne Kr Xe
);

# CH3 is not a normal stable gas. If entered, treat it as CH4.
my $TreatCH3AsCH4 = "Yes";

# Optional binary interaction parameters for PR mixing.
# Default is zero for all gas pairs. Add values here when you have MS/project
# calibrated values, for example:
#   my %BinaryInteraction = ( "CO2-CH4" => 0.10 );
my %BinaryInteraction = ();

my $OutputStudyTableName = "FugConv.std";

# ================= Gas database =================

# Tc in K, Pc in kPa, omega is acentric factor.
my %Critical = (
    N2     => { Tc => 126.192, Pc => 3395.8,  omega => 0.0372, note => "" },
    O2     => { Tc => 154.581, Pc => 5043.0,  omega => 0.0210, note => "" },
    CO2    => { Tc => 304.128, Pc => 7377.3,  omega => 0.2239, note => "" },
    CO     => { Tc => 132.86,  Pc => 3499.0,  omega => 0.0490, note => "" },
    CH4    => { Tc => 190.564, Pc => 4599.2,  omega => 0.0114, note => "methane" },
    C2H2   => { Tc => 308.32,  Pc => 6138.0,  omega => 0.1870, note => "acetylene" },
    C2H4   => { Tc => 282.35,  Pc => 5041.8,  omega => 0.0865, note => "ethylene" },
    C2H6   => { Tc => 305.32,  Pc => 4872.0,  omega => 0.0995, note => "ethane" },
    C3H8   => { Tc => 369.83,  Pc => 4248.0,  omega => 0.1520, note => "propane" },
    H2     => { Tc => 33.145,  Pc => 1296.4,  omega => -0.2190, note => "hydrogen; PR estimate" },
    H2O    => { Tc => 647.096, Pc => 22064.0, omega => 0.3440, note => "water vapor; use vapor partial pressure for humidity" },
    NO     => { Tc => 180.0,   Pc => 6485.0,  omega => 0.58,   note => "approximate NO constants" },
    NO2    => { Tc => 431.4,   Pc => 10100.0, omega => 0.85,   note => "NO2 may dimerize; EOS estimatec" },
    N2O    => { Tc => 309.52,  Pc => 7245.0,  omega => 0.162,  note => "nitrous oxide" },
    SO2    => { Tc => 430.64,  Pc => 7884.0,  omega => 0.251,  note => "" },
    H2S    => { Tc => 373.1,   Pc => 8963.0,  omega => 0.100,  note => "hydrogen sulfide" },
    NH3    => { Tc => 405.4,   Pc => 11333.0, omega => 0.256,  note => "ammonia; associating gas estimatej" },
    AR     => { Tc => 150.687, Pc => 4863.0,  omega => -0.002, note => "argons" },
    HE     => { Tc => 5.195,   Pc => 227.5,   omega => -0.385, note => "helium; PR is rough for quantum gas" },
    NE     => { Tc => 44.492,  Pc => 2678.6,  omega => -0.039, note => "neon" },
    KR     => { Tc => 209.48,  Pc => 5525.0,  omega => 0.000,  note => "krypton" },
    XE     => { Tc => 289.73,  Pc => 5842.0,  omega => 0.004,  note => "xenon" },
);

my %Alias = (
    ARGON => "AR",
    HELIUM => "HE",
    NEON => "NE",
    KRYPTON => "KR",
    XENON => "XE",
    METHANE => "CH4",
    ACETYLENE => "C2H2",
    ETHYLENE => "C2H4",
    ETHANE => "C2H6",
    PROPANE => "C3H8",
);

# ================= Main =================

die "Temperature_K must be positive.\n" unless $Temperature_K > 0.0;
die "Pressure list is empty.\n" unless @TotalPressure_List;
my $pressure_unit = normalize_pressure_unit($PressureInputUnit);

my @cases = build_cases();
my $std = Documents->New(unique_document_name($OutputStudyTableName, "std"));
write_headings($std);

my $row = 1;
foreach my $pressure_input (@TotalPressure_List) {
    die "Total pressure must be positive.\n" unless $pressure_input > 0.0;
    my $pressure_kpa = pressure_to_kpa($pressure_input, $pressure_unit);
    my $pressure_mpa = $pressure_kpa / 1000.0;

    foreach my $case (@cases) {
        my $result = mixture_fugacity_pr(
            $Temperature_K, $pressure_kpa, $case->{Components}
        );

        for (my $i = 0; $i < @{$case->{Components}}; ++$i) {
            my $component = $case->{Components}->[$i];
            my $yi = $component->{MoleFraction};
            my $partial_pressure_kpa = $yi * $pressure_kpa;
            my $phi = $result->{Phi}->[$i];
            my $ln_phi = $result->{LnPhi}->[$i];
            my $fugacity_kpa = $phi * $partial_pressure_kpa;

            write_row(
                $std, $row, $case->{Name}, $component, $Temperature_K,
                $pressure_mpa, $pressure_kpa, $partial_pressure_kpa,
                $phi, $ln_phi, $fugacity_kpa, $result
            );
            ++$row;
        }
    }
}

write_notes($std);
$std->Save;

print "\nFugacity conversion finished.\n";
print "Temperature: $Temperature_K K\n";
print "Pressure input unit: $pressure_unit\n";
print "Fugacity output unit for MS Sorption FixedPressure: kPa\n";
print "Output table: " . $std->Name . "\n";

# ================= Peng-Robinson EOS =================

sub mixture_fugacity_pr {
    my ($temperature, $pressure_kpa, $components) = @_;

    my $r = 8.314462618;          # kPa L mol-1 K-1
    my $sqrt2 = sqrt(2.0);

    my @a;
    my @b;
    my @y;

    for (my $i = 0; $i < @$components; ++$i) {
        my $gas = $components->[$i]->{EOSGas};
        my $data = $Critical{$gas};
        my ($ai, $bi) = pure_pr_ab(
            $temperature, $data->{Tc}, $data->{Pc}, $data->{omega}, $r
        );
        push @a, $ai;
        push @b, $bi;
        push @y, $components->[$i]->{MoleFraction};
    }

    my $amix = 0.0;
    my @aij;
    for (my $i = 0; $i < @a; ++$i) {
        for (my $j = 0; $j < @a; ++$j) {
            my $kij = binary_kij($components->[$i]->{EOSGas},
                                 $components->[$j]->{EOSGas});
            $aij[$i][$j] = sqrt($a[$i] * $a[$j]) * (1.0 - $kij);
            $amix += $y[$i] * $y[$j] * $aij[$i][$j];
        }
    }

    my $bmix = 0.0;
    for (my $i = 0; $i < @b; ++$i) {
        $bmix += $y[$i] * $b[$i];
    }

    my $a_red = $amix * $pressure_kpa
              / ($r * $r * $temperature * $temperature);
    my $b_red = $bmix * $pressure_kpa / ($r * $temperature);

    my @roots = pr_z_roots($a_red, $b_red);
    my $z = largest_gas_root(\@roots, $b_red);

    my $arg1 = $z - $b_red;
    my $arg2 = ($z + (1.0 + $sqrt2) * $b_red)
             / ($z + (1.0 - $sqrt2) * $b_red);
    die "Invalid PR logarithm argument. Check pressure, temperature, or gas constants.\n"
        if $arg1 <= 0.0 || $arg2 <= 0.0;

    my @ln_phi;
    my @phi;
    for (my $i = 0; $i < @a; ++$i) {
        my $sum_y_aij = 0.0;
        for (my $j = 0; $j < @a; ++$j) {
            $sum_y_aij += $y[$j] * $aij[$i][$j];
        }

        my $bi_over_b = $b[$i] / $bmix;
        my $term = 2.0 * $sum_y_aij / $amix - $bi_over_b;
        my $ln_phi_i =
            $bi_over_b * ($z - 1.0)
            - log($arg1)
            - $a_red / (2.0 * $sqrt2 * $b_red) * $term * log($arg2);

        push @ln_phi, $ln_phi_i;
        push @phi, exp($ln_phi_i);
    }

    return {
        Z => $z,
        A => $a_red,
        B => $b_red,
        Phi => \@phi,
        LnPhi => \@ln_phi,
    };
}

sub pure_pr_ab {
    my ($temperature, $tc, $pc, $omega, $r) = @_;

    my $tr = $temperature / $tc;
    my $kappa = 0.37464 + 1.54226 * $omega - 0.26992 * $omega * $omega;
    my $alpha = (1.0 + $kappa * (1.0 - sqrt($tr))) ** 2;

    my $a = 0.45724 * $r * $r * $tc * $tc / $pc * $alpha;
    my $b = 0.07780 * $r * $tc / $pc;

    return ($a, $b);
}

sub binary_kij {
    my ($gas_a, $gas_b) = @_;

    my $key1 = $gas_a . "-" . $gas_b;
    my $key2 = $gas_b . "-" . $gas_a;
    return $BinaryInteraction{$key1} if exists $BinaryInteraction{$key1};
    return $BinaryInteraction{$key2} if exists $BinaryInteraction{$key2};
    return 0.0;
}

sub pr_z_roots {
    my ($a_red, $b_red) = @_;

    my $c2 = -1.0 + $b_red;
    my $c1 = $a_red - 3.0 * $b_red * $b_red - 2.0 * $b_red;
    my $c0 = -($a_red * $b_red - $b_red * $b_red - $b_red ** 3);

    return cubic_real_roots($c2, $c1, $c0);
}

sub cubic_real_roots {
    my ($a, $b, $c) = @_;

    my $p = $b - $a * $a / 3.0;
    my $q = 2.0 * $a ** 3 / 27.0 - $a * $b / 3.0 + $c;
    my $disc = ($q / 2.0) ** 2 + ($p / 3.0) ** 3;
    my @roots;

    if ($disc >= -1.0e-14) {
        $disc = 0.0 if $disc < 0.0;
        my $sqrt_disc = sqrt($disc);
        my $u = cbrt(-$q / 2.0 + $sqrt_disc);
        my $v = cbrt(-$q / 2.0 - $sqrt_disc);
        push @roots, $u + $v - $a / 3.0;
    } else {
        my $radius = 2.0 * sqrt(-$p / 3.0);
        my $cos_arg = (3.0 * $q / (2.0 * $p)) * sqrt(-3.0 / $p);
        $cos_arg = 1.0 if $cos_arg > 1.0;
        $cos_arg = -1.0 if $cos_arg < -1.0;
        my $theta = acos($cos_arg);
        for (my $k = 0; $k < 3; ++$k) {
            push @roots,
                $radius * cos(
                    ($theta - 2.0 * 3.141592653589793 * $k) / 3.0
                ) - $a / 3.0;
        }
    }

    @roots = sort { $a <=> $b } @roots;
    return @roots;
}

sub largest_gas_root {
    my ($roots, $b_red) = @_;

    my @valid = grep { defined($_) && $_ > $b_red } @$roots;
    die "No valid gas-like PR Z root found.\n" unless @valid;
    @valid = sort { $a <=> $b } @valid;
    return $valid[-1];
}

sub cbrt {
    my ($x) = @_;
    return 0.0 if abs($x) < 1.0e-30;
    return $x > 0.0 ? exp(log($x) / 3.0) : -exp(log(-$x) / 3.0);
}

sub acos {
    my ($x) = @_;
    $x = 1.0 if $x > 1.0;
    $x = -1.0 if $x < -1.0;
    return atan2(sqrt(1.0 - $x * $x), $x);
}

# ================= Input/case handling =================

sub build_cases {
    my @cases;

    if ($CalculateAllSupportedPureGases =~ /^y/i) {
        foreach my $gas (@PureGasList) {
            my @pure = normalize_components([{ Gas => $gas, Ratio => 1.0 }]);
            push @cases, {
                Name => "Pure_" . canonical_gas_name($gas),
                Components => \@pure,
            };
        }
    } else {
        my @mix = normalize_components(\@Components);
        push @cases, {
            Name => component_case_name(\@mix),
            Components => \@mix,
        };
    }

    return @cases;
}

sub normalize_pressure_unit {
    my ($unit) = @_;
    $unit = uc($unit || "MPA");
    $unit =~ s/\s+//g;

    return "MPa" if $unit eq "MPA";
    return "kPa" if $unit eq "KPA";

    die "PressureInputUnit must be MPa or kPa.\n";
}

sub pressure_to_kpa {
    my ($pressure, $unit) = @_;

    return $pressure * 1000.0 if $unit eq "MPa";
    return $pressure if $unit eq "kPa";

    die "PressureInputUnit must be MPa or kPa.\n";
}

sub normalize_components {
    my ($input_components) = @_;

    my @normalized;
    my $sum_ratio = 0.0;

    foreach my $component (@$input_components) {
        my $input_gas = canonical_gas_name($component->{Gas});
        my $ratio = $component->{Ratio};
        die "Ratio for $input_gas must be positive.\n"
            unless defined($ratio) && $ratio > 0.0;

        my ($eos_gas, $alias_note) = resolve_eos_gas($input_gas);
        my $data = $Critical{$eos_gas};
        die "No gas data for '$input_gas'. Add it to Critical.\n"
            unless defined $data;

        push @normalized, {
            InputGas => $input_gas,
            EOSGas => $eos_gas,
            Ratio => $ratio,
            MoleFraction => 0.0,
            Note => join_notes($alias_note, $data->{note}),
        };
        $sum_ratio += $ratio;
    }

    die "Sum of component ratios must be positive.\n" unless $sum_ratio > 0.0;

    foreach my $component (@normalized) {
        $component->{MoleFraction} = $component->{Ratio} / $sum_ratio;
    }

    return @normalized;
}

sub canonical_gas_name {
    my ($gas) = @_;

    $gas = uc($gas);
    $gas =~ s/\s+//g;
    return $Alias{$gas} if exists $Alias{$gas};
    return $gas;
}

sub resolve_eos_gas {
    my ($gas) = @_;

    if ($gas eq "CH3") {
        die "CH3 is not a normal stable gas. Use CH4 or set TreatCH3AsCH4=Yes.\n"
            unless $TreatCH3AsCH4 =~ /^y/i;
        return ("CH4", "Input CH3 treated as CH4 methane");
    }

    return ($gas, "");
}

sub component_case_name {
    my ($components) = @_;
    my @items;
    foreach my $component (@$components) {
        push @items, $component->{InputGas} . $component->{Ratio};
    }
    return join("_", @items);
}

# ================= Output =================

sub write_headings {
    my ($table) = @_;

    $table->ColumnHeading(0) = "Case";
    $table->ColumnHeading(1) = "InputGas";
    $table->ColumnHeading(2) = "EOSGas";
    $table->ColumnHeading(3) = "Temperature_K";
    $table->ColumnHeading(4) = "TotalPressure_MPa";
    $table->ColumnHeading(5) = "TotalPressure_kPa";
    $table->ColumnHeading(6) = "InputRatio";
    $table->ColumnHeading(7) = "MoleFraction";
    $table->ColumnHeading(8) = "PartialPressure_kPa";
    $table->ColumnHeading(9) = "Phi";
    $table->ColumnHeading(10) = "LnPhi";
    $table->ColumnHeading(11) = "Fugacity_kPa";
    $table->ColumnHeading(12) = "Fugacity_MPa";
    $table->ColumnHeading(13) = "Mixture_Z";
    $table->ColumnHeading(14) = "Mixture_A";
    $table->ColumnHeading(15) = "Mixture_B";
    $table->ColumnHeading(16) = "Sorption_Input_Value_kPa";
    $table->ColumnHeading(17) = "Ready_Line";
    $table->ColumnHeading(18) = "Note";
}

sub write_row {
    my ($table, $row, $case_name, $component, $temperature,
        $pressure_mpa, $pressure_kpa, $partial_pressure_kpa,
        $phi, $ln_phi, $fugacity_kpa, $result) = @_;

    $table->Cell($row, 0) = $case_name;
    $table->Cell($row, 1) = $component->{InputGas};
    $table->Cell($row, 2) = $component->{EOSGas};
    $table->Cell($row, 3) = $temperature;
    $table->Cell($row, 4) = $pressure_mpa;
    $table->Cell($row, 5) = $pressure_kpa;
    $table->Cell($row, 6) = $component->{Ratio};
    $table->Cell($row, 7) = $component->{MoleFraction};
    $table->Cell($row, 8) = $partial_pressure_kpa;
    $table->Cell($row, 9) = $phi;
    $table->Cell($row, 10) = $ln_phi;
    $table->Cell($row, 11) = $fugacity_kpa;
    $table->Cell($row, 12) = $fugacity_kpa / 1000.0;
    $table->Cell($row, 13) = $result->{Z};
    $table->Cell($row, 14) = $result->{A};
    $table->Cell($row, 15) = $result->{B};
    $table->Cell($row, 16) = $fugacity_kpa;
    $table->Cell($row, 17) =
        '$task->Fugacity($sorbate) = ' . $fugacity_kpa . '; # kPa';
    $table->Cell($row, 18) = $component->{Note};
}

sub write_notes {
    my ($table) = @_;

    $table->Cell(0, 20) = "Meaning";
    $table->Cell(0, 21) =
        "Use Fugacity_kPa as the value for Sorption FixedPressure Fugacity(component).";
    $table->Cell(1, 20) = "EOS";
    $table->Cell(1, 21) =
        "Peng-Robinson mixture EOS. Binary interaction kij defaults to 0.";
    $table->Cell(2, 20) = "No MS calculation";
    $table->Cell(2, 21) =
        "This script only converts pressure/composition to fugacity.";
    $table->Cell(3, 20) = "PressureInputUnit";
    $table->Cell(3, 21) = $pressure_unit;
}

sub join_notes {
    my ($a, $b) = @_;
    return $b if (!defined($a) || $a eq "") && defined($b);
    return $a if (!defined($b) || $b eq "") && defined($a);
    return "" if (!defined($a) || $a eq "") && (!defined($b) || $b eq "");
    return "$a; $b";
}

sub unique_document_name {
    my ($name, $extension) = @_;
    $name =~ s/\.$extension$//i;
    $name =~ s/[^A-Za-z0-9]+//g;
    $name = "FugConv" if $name eq "";

    for (my $attempt = 0; $attempt < 1000; ++$attempt) {
        my $suffix = $attempt == 0 ? "" : $attempt;
        my $base = substr($name . $suffix, 0, 24);
        my $candidate = "$base.$extension";
        my $existing = eval { $Documents{$candidate}; };
        return $candidate unless defined $existing;
    }

    die "Could not create a unique .$extension document name.\n";
}
