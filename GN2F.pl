#!perl
use strict;
use warnings;
use Getopt::Long;

# -------------------------------------------------------------------------
# Forcite ExternalForceStrength converter
# Author Jingsai Cheng
# E-mail jingsaicheng@gmail.com
# Supports two input modes:
#   1. mode = "GPa": input target pressure P(GPa), area side lengths A/B(Angstrom).
#   2. mode = "N"  : input total force F(N), atom count n.
#
# Forcite ExternalForceStrength unit:
#   kcal/mol/Angstrom per atom
#
# Conversion:
#   1 kcal/mol/Angstrom = 6.9477e-11 N
#   P(GPa) = 6.9477 * n * f / S
#   f = P(GPa) * S(Angstrom^2) / (6.9477 * n)
# -------------------------------------------------------------------------

# ================= User input, edit here =================

my $Mode = "GPa";      # "GPa" or "N"
my $Value = 0.25; # If Mode=GPa: pressure in GPa. If Mode=N: total force in N.
my $AtomCount = 55;   # n, atoms in the ExternalForceSet

# Area S = LengthA * LengthB, unit Angstrom^2.
# For a square face, set LengthA = LengthB.
my $LengthA_A = 10.0;
my $LengthB_A = 50.0;

# Optional: directly set area. Leave 0 to use LengthA_A * LengthB_A.
my $Area_A2 = 0.0;

# ================= Command-line override =================
# Example:
#   perl external_force_converter.pl --mode GPa --value 0.00101325 --n 100 --a 50 --b 50
#   perl external_force_converter.pl --mode N   --value 1e-9       --n 100

GetOptions(
    "mode=s"  => \$Mode,
    "value=f" => \$Value,
    "n=i"     => \$AtomCount,
    "a=f"     => \$LengthA_A,
    "b=f"     => \$LengthB_A,
    "s=f"     => \$Area_A2,
) or die "Bad input options.\n";

# ================= Constants =================

my $N_PER_KCAL_MOL_A = 6.9477e-11;
my $ATM_TO_GPA = 0.000101325;

# ================= Calculation =================

die "AtomCount n must be > 0.\n" unless defined($AtomCount) && $AtomCount > 0;

my $S_A2 = $Area_A2 > 0.0 ? $Area_A2 : $LengthA_A * $LengthB_A;
die "Area S must be > 0 Angstrom^2.\n" unless defined($S_A2) && $S_A2 > 0.0;

my ($P_GPa, $F_total_N, $F_total_kcal_mol_A, $F_per_atom_N, $ExternalForceStrength);

if ($Mode =~ /^GPa$/i) {
    $P_GPa = $Value;
    $F_total_N = $P_GPa * $S_A2 * 1.0e-11;
    $F_total_kcal_mol_A = $F_total_N / $N_PER_KCAL_MOL_A;
    $ExternalForceStrength = $F_total_kcal_mol_A / $AtomCount;
    $F_per_atom_N = $F_total_N / $AtomCount;
} elsif ($Mode =~ /^N$/i) {
    $F_total_N = $Value;
    $P_GPa = $F_total_N / ($S_A2 * 1.0e-11);
    $F_total_kcal_mol_A = $F_total_N / $N_PER_KCAL_MOL_A;
    $ExternalForceStrength = $F_total_kcal_mol_A / $AtomCount;
    $F_per_atom_N = $F_total_N / $AtomCount;
} else {
    die "Mode must be 'GPa' or 'N'.\n";
}

my $P_atm = $P_GPa / $ATM_TO_GPA;

print "\nExternal force conversion result\n";
print "--------------------------------\n";
printf "Mode                         : %s\n", $Mode;
printf "Atom count n                 : %d\n", $AtomCount;
printf "Area S                       : %.3f Angstrom^2\n", $S_A2;
printf "Pressure                     : %.6g GPa = %.3f atm\n", $P_GPa, $P_atm;
printf "Total force                  : %.6g N\n", $F_total_N;
printf "Total force                  : %.3f kcal/mol/Angstrom\n", $F_total_kcal_mol_A;
printf "Force per atom               : %.6g N\n", $F_per_atom_N;
printf "ExternalForceStrength input  : %.3f kcal/mol/Angstrom per atom\n", $ExternalForceStrength;
print "\nForcite setting example:\n";
printf "ExternalForceStrength => %.3f,\n", $ExternalForceStrength;
print  "ExternalForceX => 0,\n";
print  "ExternalForceY => 0,\n";
print  "ExternalForceZ => -1,\n";
print  "CounterExternalForce => \"No\",\n";

